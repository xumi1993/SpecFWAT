! GPU preparation uses the SPECFEM FK conventions (GPL-3.0-or-later).
! The common reflection solve below is adapted from couple_with_injection.f90.
module fk_gpu_backend
  use iso_c_binding
  use config, only: noderank
  use specfem_par, only: CUSTOM_REAL, myrank, PI, TINYVAL, IMAIN, ACOUSTIC_SIMULATION, &
    NGLLSQUARE, num_abs_boundary_faces, abs_boundary_ispec, abs_boundary_ijk, &
    abs_boundary_normal, ibool, xstore, ystore, zstore, kappastore, mustore, deltat
  use specfem_par_coupling, only: vp=>alpha_FK, vs=>beta_FK, rho=>rho_FK, H=>h_FK, nlayer, &
    kpsv=>type_kpsv_fk, amplitude_fk, phi_FK, xx0, yy0, zz0, tt0, Z_REF_for_FK, &
    npt, NF_FOR_STORING, NF_FOR_FFT, NP_RESAMP, Veloc_FK, Tract_FK, ipt_table
  use specfem_par_elastic, only: ispec_is_elastic
  implicit none
  private
  public :: compute_fk_gpu
  interface
    integer(c_int) function fk_cuda_compute(nl, np, nf, local_rank, layers, eta_a, eta_b, &
        gamma, emat, bottom, cache, geom, point_layer, acoustic, opts, vel, tract) bind(C)
      import
      integer(c_int), value :: nl, np, nf, local_rank
      real(c_double), intent(in) :: layers(*), geom(*), opts(*)
      complex(c_double_complex), intent(in) :: eta_a(*), eta_b(*), gamma(*), emat(*), bottom(*), cache(*)
      integer(c_int), intent(in) :: point_layer(*), acoustic(*)
      real(c_double), intent(out) :: vel(*), tract(*)
    end function
  end interface
contains
  subroutine compute_fk_gpu(ray_p, Tg, DF_FK)
    real(kind=CUSTOM_REAL), intent(in) :: ray_p, Tg, DF_FK
    integer, parameter :: CUSTOM_CMPLX=8
    real(kind=CUSTOM_REAL), parameter :: THRESHOLD_VS=1.e-6
    integer :: nf2, ii, i, j, ilayer, ilayer_ac, nn, iface, igll, ispec, ipt, iglob, k, status
    real(kind=CUSTOM_REAL) :: om, C_1, two_mul, xi, mu, kap, height
    logical :: have_fluid_layer
    real(kind=CUSTOM_REAL), allocatable :: fvec(:)
    complex(kind=8) :: eta_alpha(nlayer), eta_beta(nlayer), gamma0(nlayer), gamma1(nlayer)
    complex(kind=8) :: E_mat(4,4), N_mat(4,4), Pmat(4,4), Qmat(2,2), Qmat_I(2,2), N1_mat(2,2)
    complex(kind=8) :: C_3, eta_p, eta_s, a,b,c,d,delta_mat,MM, vec(4)
    complex(kind=8), allocatable :: coeff(:,:)
    complex(c_double_complex), allocatable :: bottom(:,:), cache(:,:,:)
    real(c_double), allocatable :: layers(:,:), geom(:,:), velocity(:,:,:), traction(:,:,:)
    integer(c_int), allocatable :: point_layer(:), acoustic(:)
    real(c_double) :: opts(12)

    nf2 = NF_FOR_STORING + 1
    allocate(fvec(nf2), coeff(2,nf2), bottom(4,nf2), cache(4,nlayer,nf2), &
      layers(3,nlayer), geom(10,npt), point_layer(npt), acoustic(npt), &
      velocity(3,npt,NF_FOR_STORING), traction(3,npt,NF_FOR_STORING))
    do ii=1,nf2
      fvec(ii) = (ii-1)*DF_FK
    enddo
    nn = int(-tt0/deltat)
    ilayer_ac = 0
    C_1 = 0
    C_3 = 0
    ! check if top layers are in fluid material
    have_fluid_layer = .false.
    do j = nlayer,1,-1
      if (vs(j) < THRESHOLD_VS) then
        ilayer_ac = j
        have_fluid_layer = .true.
        exit
      endif
    enddo

    ! check no elastic layers for ilayer < ilayer_ac
    if (have_fluid_layer .and. (ilayer_ac > 1)) then
      if (sum(vs(1:ilayer_ac)) / ilayer_ac > THRESHOLD_VS) then
        if (myrank == 0) print *,'Also check the FK model, make sure fluid layers are on the top of elastic layers'
        call exit_MPI(myrank, 'Invalid FK fluid layer ordering')
      endif
    endif

    ! make sure ACOUSTIC_SIMULATION is enabled
    if (have_fluid_layer .and. (.not. ACOUSTIC_SIMULATION)) then
      if (myrank == 0) then
        print *,'FK model contains fluid layers, but ACOUSTIC_SIMULATION is not enabled'
      endif
      call exit_MPI(myrank, 'FK fluid layers require acoustic mesh')
    endif

    if (have_fluid_layer) then
      if (myrank == 0) write(IMAIN,*) 'FK simulation: acoustic + elastic'
    else
      if (myrank == 0) write(IMAIN,*) 'FK simulation: elastic'
    endif

    ! compute temporary variables
    do i = 1,nlayer
      eta_alpha(i) = -cmplx(0,1) * sqrt( 1.0 / vp(i)**2 - ray_p**2 )
      gamma0(i) = 2.0 * vs(i)**2 * ray_p**2
      if (vs(i) < THRESHOLD_VS) then
        gamma1(i) = 0.
        eta_beta(i) = 0.
      else
        gamma1(i) = 1.0 - 1.0/gamma0(i)
        eta_beta(i) = -cmplx(0,1) * sqrt( 1.0 / vs(i)**2 - ray_p**2 )
      endif
    enddo

    if (myrank == 0) then
      write(IMAIN,*) '    starting from ',nn,' points before time 0'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! amplitude in half space
    if (kpsv == 1) then  ! P-SV
      C_3 = amplitude_fk * cmplx(0,1.) * ray_p * vp(nlayer)      ! amp. of incoming P in the bot. layer
      eta_p = sqrt(1.0/vp(nlayer)**2 - ray_p**2)                 ! vertical slowness for lower layer
      if (myrank == 0) write(IMAIN,*) '  Incoming P : C_3,  ray_p, eta = ', C_3, ray_p, eta_p
    else
      ! SV-wave
      ! for C_2 = sin(inc) (u=[cos(inc), sin(inc)])
      C_1 = amplitude_fk * ray_p * vs(nlayer)                   ! amp. of incoming S in the bot. layer
      eta_s = sqrt(1.0/vs(nlayer)**2 - ray_p**2)                ! vertical slowness for lower layer

      eta_p = eta_s ! Explicit incident slowness for SV (legacy code leaves eta_p unset).
      if (myrank == 0 ) write(IMAIN,*) '  Incoming S :  C_1,  ray_p, eta = ', C_1, ray_p, eta_s
    endif

    !E matrix
    ! initializes matrix
    two_mul = 2.0 * rho(nlayer) * vs(nlayer) * vs(nlayer)
    E_mat(:,:) = (1.0,0.0)

    ! Tong et al. (2014), appendix (A10) E_0:
    ! note: E_mat is not omega dependent
    E_mat(1,1) =  eta_beta(nlayer) / ray_p
    E_mat(1,2) = -E_mat(1,1)
    E_mat(2,3) =  eta_alpha(nlayer) / ray_p
    E_mat(2,4) = -E_mat(2,3)
    E_mat(3,1) = two_mul * gamma1(nlayer)
    E_mat(3,2) = E_mat(3,1)
    E_mat(3,3) = two_mul * eta_alpha(nlayer) / ray_p
    E_mat(3,4) = -E_mat(3,3)
    E_mat(4,1) = two_mul * eta_beta(nlayer) / ray_p
    E_mat(4,2) = -E_mat(4,1)
    E_mat(4,3) = E_mat(3,1)
    E_mat(4,4) = E_mat(3,1)

    ! now loop every frequency to determine coefs in half space
    do ii = 1,nf2
      om = 2.0 * PI * fvec(ii)

      ! apply propagation matrix in elastic layers
      N_mat = E_mat
      ilayer = 1
      if (have_fluid_layer) ilayer = ilayer_ac + 1
      do i = nlayer-1,ilayer,-1
        call fk_propagator_psv(om,eta_alpha(i),eta_beta(i),rho(i), &
                              vs(i),H(i),ray_p,gamma1(i),Pmat)
        N_mat = matmul(Pmat,N_mat) * gamma0(i)
      enddo

      ! apply propagation matrix in acoustic layers
      if (have_fluid_layer) then
        Qmat_I(:,:) = 0.0_CUSTOM_CMPLX
        Qmat_I(1,1) = cmplx(1.0,0.0,kind=CUSTOM_CMPLX)
        Qmat_I(2,2) = cmplx(1.0,0.0,kind=CUSTOM_CMPLX)
        do j = ilayer_ac,1,-1
          call fk_propagator_ac(om,eta_alpha(j),rho(j),H(j),ray_p,Qmat(:,:))
          Qmat_I = matmul(Qmat,Qmat_I)
        enddo
        Qmat = Qmat_I
      endif

      !determine coefs in half space
      if (.not. have_fluid_layer) then
        ! inverse matrix
        a = N_mat(3,2); b = N_mat(3,4); c = N_mat(4,2); d = N_mat(4,4)
        delta_mat = a*d - b*c
        if (abs(delta_mat) > TINYVAL) then
          if (kpsv == 1) then
            coeff(1,ii) = -(d*N_mat(3,3) - b*N_mat(4,3)) / delta_mat * C_3
            coeff(2,ii) = -(-c*N_mat(3,3) + a*N_mat(4,3)) / delta_mat * C_3
          else
            coeff(1,ii) = -(d*N_mat(3,1) - b*N_mat(4,1)) / delta_mat * C_1
            coeff(2,ii) = -(-c*N_mat(3,1) + a*N_mat(4,1)) / delta_mat * C_1
          endif
        else
          coeff(1,ii) = (0.d0,0.d0)
          coeff(2,ii) = (0.d0,0.d0)
        endif
      else
        N1_mat(1,1) = N_mat(3,2)
        N1_mat(1,2) = N_mat(3,4)
        if (kpsv == 1) then
          N1_mat(2,1) = Qmat(2,1) * N_mat(2,2) - Qmat(2,2) * N_mat(4,2)
          N1_mat(2,2) = Qmat(2,1) * N_mat(2,4) - Qmat(2,2) * N_mat(4,4)
          MM = Qmat(2,1) * N_mat(2,3) - Qmat(2,2) * N_mat(4,3)
        else
          N1_mat(2,1) = Qmat(2,1) * N_mat(2,2) - Qmat(2,2) * N_mat(4,2)
          N1_mat(2,2) = Qmat(2,1) * N_mat(2,4) - Qmat(2,2) * N_mat(4,4)
          MM = Qmat(2,1) * N_mat(2,1) - Qmat(2,2) * N_mat(4,1)
        endif
        a = N1_mat(1,1); b = N1_mat(1,2); c = N1_mat(2,1); d = N1_mat(2,2)
        delta_mat = a*d - b*c
        if (abs(delta_mat) > TINYVAL) then
          if (kpsv == 1) then
            coeff(1,ii) = -( d * N_mat(3,3) - b * MM) / delta_mat * C_3
            coeff(2,ii) = -(-c * N_mat(3,3) + a * MM) / delta_mat * C_3
          else
            coeff(1,ii) = -( d * N_mat(3,1) - b * MM) / delta_mat * C_1
            coeff(2,ii) = -(-c * N_mat(3,1) + a * MM) / delta_mat * C_1
          endif
        else
          coeff(1,ii) = (0.d0,0.d0)
          coeff(2,ii) = (0.d0,0.d0)
        endif
      endif

        ! if (ii == 10 .and. myrank == 0) then
        !   print *,'debug: Rayleigh coeff ',coeff(1,ii),coeff(2,ii),delta_mat
        !   print *,'N_mat'
        !   print *,N_mat
        ! endif
    enddo


    ! Cache the state at the bottom of each layer, once per frequency.
    ! This removes complete-layer propagation from the point-frequency kernel.
    do ii=1,nf2
      om = 2.0 * PI * fvec(ii)
      vec = 0
      if (kpsv == 1) then
        vec(2)=coeff(1,ii); vec(3)=C_3; vec(4)=coeff(2,ii)
      else
        vec(1)=C_1; vec(3)=coeff(1,ii); vec(4)=coeff(2,ii)
      endif
      bottom(:,ii)=vec
      N_mat=E_mat
      cache(:,nlayer,ii)=matmul(N_mat,vec)
      do j=nlayer-1,ilayer_ac+1,-1
        cache(:,j,ii)=matmul(N_mat,vec)
        call fk_propagator_psv(om,eta_alpha(j),eta_beta(j),rho(j),vs(j),H(j),ray_p,gamma1(j),Pmat)
        N_mat=gamma0(j)*matmul(Pmat,N_mat)
      enddo
      if (have_fluid_layer) then
        vec=matmul(N_mat,vec)
        vec(1)=vec(2)
        vec(2)=-vec(4)
        do j=ilayer_ac,1,-1
          cache(:,j,ii)=vec
          call fk_propagator_ac(om,eta_alpha(j),rho(j),H(j),ray_p,Qmat)
          vec(1:2)=matmul(Qmat,vec(1:2))
        enddo
      endif
    enddo
    layers(1,:)=rho; layers(2,:)=vs; layers(3,:)=H
    ipt=0
    do iface=1,num_abs_boundary_faces
      ispec=abs_boundary_ispec(iface)
      do igll=1,NGLLSQUARE
        ipt=ipt+1
        ipt_table(igll,iface)=ipt
        i=abs_boundary_ijk(1,igll,iface)
        j=abs_boundary_ijk(2,igll,iface)
        k=abs_boundary_ijk(3,igll,iface)
        iglob=ibool(i,j,k,ispec)
        geom(1:3,ipt)=[xstore(iglob),ystore(iglob),zstore(iglob)-Z_REF_for_FK]
        geom(7:9,ipt)=abs_boundary_normal(:,igll,iface)
        ilayer=nlayer
        if (geom(3,ipt) > 0) then
          do j=nlayer-1,1,-1
            if (geom(3,ipt) <= sum(H(j:nlayer-1))) then
              ilayer=j
              exit
            endif
          enddo
        endif
        height=real(geom(3,ipt),CUSTOM_REAL)-sum(H(ilayer+1:nlayer-1))
        geom(10,ipt)=height
        point_layer(ipt)=ilayer-1
        if (ispec_is_elastic(ispec)) then
          if (ilayer <= ilayer_ac) call exit_MPI(myrank, 'Elastic FK point in fluid layer')
          acoustic(ipt)=0
          mu=mustore(i,abs_boundary_ijk(2,igll,iface),k,ispec)
          kap=kappastore(i,abs_boundary_ijk(2,igll,iface),k,ispec)
          xi=mu/(kap+4.0/3.0*mu)
          geom(4,ipt)=1.0-2.0*xi
          geom(5,ipt)=(1.0-xi)*mu
          geom(6,ipt)=(3.0*kap-2.0*mu)/(6.0*kap+2.0*mu)
        else
          if (ilayer > ilayer_ac+1 .or. height < 0 .or. ilayer_ac == 0) &
            call exit_MPI(myrank, 'Acoustic FK point outside fluid layers')
          acoustic(ipt)=1
          geom(4:6,ipt)=[1.d0,0.d0,0.5d0]
        endif
      enddo
    enddo
    opts=[dble(ray_p),dble(phi_FK),dble(xx0),dble(yy0),dble(zz0),dble(Tg),dble(DF_FK), &
      dble(eta_p%re),dble(deltat),dble(nn),dble(ilayer_ac),dble(1._CUSTOM_REAL/(NF_FOR_FFT*deltat))]
    status=fk_cuda_compute(nlayer,npt,NF_FOR_FFT,noderank,layers,eta_alpha,eta_beta, &
      gamma1,E_mat,bottom,cache,geom,point_layer,acoustic,opts,velocity,traction)
    if (status /= 0) call exit_MPI(myrank, 'CUDA FK computation failed; see CUDA diagnostic')
    Veloc_FK(:,:,1:NF_FOR_STORING)=real(velocity, CUSTOM_REAL)
    Tract_FK(:,:,1:NF_FOR_STORING)=real(traction, CUSTOM_REAL)
  end subroutine
end module
