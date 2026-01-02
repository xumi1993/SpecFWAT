module smooth_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

implicit none

contains
subroutine smooth_sem_pde(dat_in, sigma_h, sigma_v, dat, is_sph)
  use constants
  use config, only : distance_min_glob
  use specfem_par
  ! use specfem_par_elastic, only: ispec_is_elastic, &
  use specfem_par_elastic, only: &
      nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic
  ! use specfem_par_acoustic, only: ispec_is_acoustic
  ! use specfem_par_poroelastic, only: ispec_is_poroelastic
  use pml_par, only: is_CPML
  ! use wavefield_discontinuity_par,only: IS_WAVEFIELD_DISCONTINUITY

  implicit none 
  integer, parameter :: PRINT_INFO_PER_STEP = 1000
  logical, parameter :: ZERO_PML = .true.
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(in) :: dat_in
  real(kind=CUSTOM_REAL), intent(in) :: sigma_h, sigma_v
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, intent(out) :: dat
  logical, intent(in) :: is_sph
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat_bak
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: rotate_r
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    dx_elem,dy_elem,dz_elem, stemp1,stemp2,stemp3, snewtemp1,snewtemp2,snewtemp3
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dat_glob, ddat_glob
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rvol
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: rvol_local, jacobian_all
  integer :: i,j,k,l,iglob,ier,ispec,ispec_p,iphase,ispec_irreg, is_sph_int

  ! character(len=MAX_STRING_LEN) :: arg(7)
  ! character(len=MAX_STRING_LEN) :: input_dir, output_dir
  !character(len=MAX_STRING_LEN) :: prname_lp
  ! character(len=MAX_STRING_LEN*2) :: ks_file
  ! character(len=MAX_STRING_LEN*2) :: local_data_file


  !character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  !character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  ! character(len=MAX_STRING_LEN) :: kernel_name
  integer :: num_elements
  ! real t1,t2,tnow,tlast

  real(kind=CUSTOM_REAL) :: ch, cv, cmax
  real(kind=CUSTOM_REAL) :: min_val, max_val, min_val_glob, max_val_glob
  
  ! real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  ! real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  ! real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  ! real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  ! real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  real(kind=CUSTOM_REAL) :: xl,yl,zl,rl,rxl,ryl,rzl,&
    xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl 
  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3
  integer :: ntstep, istep
  double precision :: weight
  !real(kind=CUSTOM_REAL), dimension(MAX_KERNEL_NAMES) :: max_old,max_new,max_old_all,max_new_all
  !real(kind=CUSTOM_REAL), dimension(MAX_KERNEL_NAMES) :: min_old,min_new,min_old_all,min_new_all
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_vector_ext_mesh_smooth
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_recv_vector_ext_mesh_smooth
  ! logical :: BROADCAST_AFTER_READ
 
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (is_sph) then
    is_sph_int = 1
  else
    is_sph_int = 0
  endif

  !! broadcast distance_min_glob to other processors
  call bcast_all_singlecr(distance_min_glob)
  
  allocate(rotate_r(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1013')

  do ispec = 1, NSPEC_AB
    do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
      iglob = ibool(i,j,k,ispec)
      xl = xstore(iglob)
      yl = ystore(iglob)
      zl = zstore(iglob)
      rl = sqrt(xl*xl+yl*yl+zl*zl)
      rotate_r(1,i,j,k,ispec) = xl / rl
      rotate_r(2,i,j,k,ispec) = yl / rl
      rotate_r(3,i,j,k,ispec) = zl / rl
    enddo;enddo;enddo
  enddo

  ! deallocate(xstore,ystore,zstore,kappastore,mustore)
  ! deallocate(ispec_is_acoustic, ispec_is_elastic, ispec_is_poroelastic)

  !! determine ch, cv, ntstep
  !cmax = distance_min_glob ** 2 / 6.0
  cmax = distance_min_glob ** 2 / 9.0
  if (sigma_v >= sigma_h) then
    cv = cmax
    ch = cv * (sigma_h ** 2) / (sigma_v ** 2)
  else
    ch = cmax
    cv = ch * (sigma_v ** 2) / (sigma_h ** 2)
  endif
  ntstep = int(ceiling((max(sigma_h,sigma_v)**2)/(2.0*cmax)))
  
  ! if (myrank == 0) print *, 'cv=', cv, 'ch=', ch, 'ntstep=', ntstep

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! initialize time iteration
   ! set up GLL points, weights and derivation matrices for reference element
   ! (between -1,1)
  call define_derivation_matrices(xigll,yigll,zigll,wxgll,wygll,wzgll, &
                                  hprime_xx,hprime_yy,hprime_zz, &
                                  hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                                  wgllwgll_xy,wgllwgll_xz,wgllwgll_yz)
  ! define transpose of derivation matrix
  do j = 1,NGLLY
    do i = 1,NGLLX
      hprime_xxT(j,i) = hprime_xx(i,j)
      hprime_yyT(j,i) = hprime_yy(i,j)
      hprime_zzT(j,i) = hprime_zz(i,j)

      hprimewgll_xxT(j,i) = hprimewgll_xx(i,j)
    enddo
  enddo

  allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_AB),dat_bak(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(jacobian_all(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(dat_glob(NGLOB_AB))
  allocate(ddat_glob(NGLOB_AB))
  allocate(buffer_send_vector_ext_mesh_smooth( &
              max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  allocate(buffer_recv_vector_ext_mesh_smooth( &
              max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
   ! prepare assemble array
  allocate(rvol(NGLOB_AB)) 
  rvol(:) = 0.0
  allocate(rvol_local(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  do ispec = 1, NSPEC_AB
    ispec_irreg = irregular_element_number(ispec)
    do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
      weight =  wxgll(i)*wygll(j)*wzgll(k)
      if (ispec_irreg /= 0) then
        jacobianl = jacobianstore(i,j,k,ispec_irreg)
      else
        jacobianl = jacobian_regular
      endif
      jacobian_all(i,j,k,ispec) = jacobianl
      rvol_local(i,j,k,ispec) = real(dble(jacobianl)*weight,kind=CUSTOM_REAL)
      iglob = ibool(i,j,k,ispec)
      rvol(iglob) = rvol(iglob) + rvol_local(i,j,k,ispec)
    enddo;enddo;enddo
  enddo
  call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rvol, &
                       num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                       my_neighbors_ext_mesh)
  rvol(:) = 1.0 / rvol(:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read in data to be smoothed
  ! data file
  ! write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'
  ! local_data_file = trim(prname) // trim(kernel_name) // '.bin'

  ! open(unit = IIN,file = trim(local_data_file),status='old',action='read', &
  !      form ='unformatted',iostat=ier)
  ! if (ier /= 0) then
  !   print *,'Error opening data file: ',trim(local_data_file)
  !   stop 'Error opening data file'
  ! endif

  ! read(IIN) dat
  ! close(IIN)

  ! back up original
  dat(:,:,:,:) = dat_in(:,:,:,:)
  if (ZERO_PML) then
    dat_bak(:,:,:,:) = 0._CUSTOM_REAL
  else
    dat_bak(:,:,:,:) = dat_in(:,:,:,:) 
  endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! project
  dat_glob(:) = 0.0
  do ispec = 1, NSPEC_AB
    do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
      iglob = ibool(i,j,k,ispec)
      dat_glob(iglob) = dat_glob(iglob) + dat(i,j,k,ispec) * rvol_local(i,j,k,ispec)
    enddo;enddo;enddo
  enddo
  call assemble_MPI_send_smooth(NPROC,NGLOB_AB, &
          dat_glob, buffer_send_vector_ext_mesh_smooth, &
          buffer_recv_vector_ext_mesh_smooth, &
          num_interfaces_ext_mesh, max_nibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh, ibool_interfaces_ext_mesh, &
          my_neighbors_ext_mesh, &
          request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
  call assemble_MPI_w_ord_smooth(NPROC,NGLOB_AB, &
          dat_glob, buffer_recv_vector_ext_mesh_smooth, num_interfaces_ext_mesh, &
          max_nibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
          request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
          my_neighbors_ext_mesh,myrank)
  ! if (myrank == 0) print *, 'Before smoothing: '

  dat_glob(:) = dat_glob(:) * rvol(:)
  min_val = minval(dat_glob)
  max_val = maxval(dat_glob)
  call min_all_cr(min_val, min_val_glob)
  call max_all_cr(max_val, max_val_glob)
  ! if (myrank == 0) then
  ! !   ! print *, '  '//trim(kernel_name)
  !   print *, '    minval:', min_val_glob
  !   print *, '    maxval:', max_val_glob
  ! !   if (myrank == 0) call cpu_time(tlast)
  ! endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! broadcast glob array back to local array
  do ispec = 1, NSPEC_AB
    do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
      iglob = ibool(i,j,k,ispec)
      dat(i,j,k,ispec) = dat_glob(iglob)
    enddo;enddo;enddo
  enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(.not. GPU_MODE) then 
    do istep = 1, ntstep
      ddat_glob(:) = 0.0
      do iphase = 1,2
        if (iphase == 1) then
          num_elements = nspec_outer_elastic
        else
          num_elements = nspec_inner_elastic
        endif
        do ispec_p = 1,num_elements
          ispec = phase_ispec_inner_elastic(ispec_p,iphase)
          ispec_irreg = irregular_element_number(ispec)
          if (ispec_irreg /= 0) then
            call get_gradient_element(dat(:,:,:,ispec), &
              dx_elem, dy_elem, dz_elem, &
              xixstore(:,:,:,ispec), xiystore(:,:,:,ispec), xizstore(:,:,:,ispec), &
              etaxstore(:,:,:,ispec), etaystore(:,:,:,ispec), etazstore(:,:,:,ispec), &
              gammaxstore(:,:,:,ispec), gammaystore(:,:,:,ispec), gammazstore(:,:,:,ispec))
          else
            call get_gradient_element_regular(dat(:,:,:,ispec), &
              xix_regular,dx_elem, dy_elem, dz_elem)
          endif
          do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
            if (ispec_irreg /= 0) then
              xixl = xixstore(i,j,k,ispec_irreg)
              xiyl = xiystore(i,j,k,ispec_irreg)
              xizl = xizstore(i,j,k,ispec_irreg)
              etaxl = etaxstore(i,j,k,ispec_irreg)
              etayl = etaystore(i,j,k,ispec_irreg)
              etazl = etazstore(i,j,k,ispec_irreg)
              gammaxl = gammaxstore(i,j,k,ispec_irreg)
              gammayl = gammaystore(i,j,k,ispec_irreg)
              gammazl = gammazstore(i,j,k,ispec_irreg)
              jacobianl = jacobianstore(i,j,k,ispec_irreg)
            else
              xixl = xix_regular
              xiyl = xix_regular
              xizl = xix_regular
              etaxl = xix_regular
              etayl = xix_regular
              etazl = xix_regular
              gammaxl = xix_regular
              gammayl = xix_regular
              gammazl = xix_regular
              jacobianl = jacobian_regular
            endif
            if (is_sph) then
              !! spherical coordinate
              rxl = rotate_r(1,i,j,k,ispec)
              ryl = rotate_r(2,i,j,k,ispec)
              rzl = rotate_r(3,i,j,k,ispec)
              stemp1(i,j,k) = ((cv-ch) * (rxl*xixl+ryl*xiyl+rzl*xizl) * &
                (rxl*dx_elem(i,j,k)+ryl*dy_elem(i,j,k)+rzl*dz_elem(i,j,k)) +&
                ch * (xixl*dx_elem(i,j,k)+xiyl*dy_elem(i,j,k)&
                    +xizl*dz_elem(i,j,k))) * jacobianl
              stemp2(i,j,k) = ((cv-ch) * (rxl*etaxl+ryl*etayl+rzl*etazl) * &
                (rxl*dx_elem(i,j,k)+ryl*dy_elem(i,j,k)+rzl*dz_elem(i,j,k)) +&
                ch * (etaxl*dx_elem(i,j,k)+etayl*dy_elem(i,j,k)&
                    +etazl*dz_elem(i,j,k))) * jacobianl
              stemp3(i,j,k) = ((cv-ch) * (rxl*gammaxl+ryl*gammayl+rzl*gammazl) * &
                (rxl*dx_elem(i,j,k)+ryl*dy_elem(i,j,k)+rzl*dz_elem(i,j,k)) +&
                ch * (gammaxl*dx_elem(i,j,k)+gammayl*dy_elem(i,j,k)&
                    +gammazl*dz_elem(i,j,k))) * jacobianl  
            else
              !! Cartesian coordinate
              stemp1(i,j,k) = ((cv-ch) * xizl * dz_elem(i,j,k) + &
                ch * (xixl * dx_elem(i,j,k) + xiyl * dy_elem(i,j,k) + &
                xizl * dz_elem(i,j,k))) * jacobianl
              stemp2(i,j,k) = ((cv-ch) * etazl * dz_elem(i,j,k) + &
                ch * (etaxl * dx_elem(i,j,k) + etayl * dy_elem(i,j,k) + &
                etazl * dz_elem(i,j,k))) * jacobianl
              stemp3(i,j,k) = ((cv-ch) * gammazl * dz_elem(i,j,k) + &
                ch * (gammaxl * dx_elem(i,j,k) + gammayl * dy_elem(i,j,k) + &
                gammazl * dz_elem(i,j,k))) * jacobianl
            endif
          enddo;enddo;enddo
          do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
            snewtemp1(i,j,k) = 0.0
            snewtemp2(i,j,k) = 0.0
            snewtemp3(i,j,k) = 0.0
            do l = 1, NGLLX
              fac1 = hprimewgll_xx(l,i)
              snewtemp1(i,j,k) = snewtemp1(i,j,k) + stemp1(l,j,k) * fac1
              fac2 = hprimewgll_yy(l,j)
              snewtemp2(i,j,k) = snewtemp2(i,j,k) + stemp2(i,l,k) * fac2
              fac3 = hprimewgll_zz(l,k)
              snewtemp3(i,j,k) = snewtemp3(i,j,k) + stemp3(i,j,l) * fac3
            enddo
            fac1 = wgllwgll_yz(j,k)
            fac2 = wgllwgll_xz(i,k)
            fac3 = wgllwgll_xy(i,j)
            iglob = ibool(i,j,k,ispec)
            ddat_glob(iglob) = ddat_glob(iglob) - (fac1*snewtemp1(i,j,k)+&
                            fac2 * snewtemp2(i,j,k) + fac3 * snewtemp3(i,j,k))
          enddo;enddo;enddo 
        enddo  ! ispec_p = 1, num_elements 
        !! assemble MPI
        if (iphase == 1) then
          call assemble_MPI_send_smooth(NPROC,NGLOB_AB,&
            ddat_glob,buffer_send_vector_ext_mesh_smooth,&
            buffer_recv_vector_ext_mesh_smooth,&
            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
            my_neighbors_ext_mesh, &
            request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
        else
          call assemble_MPI_w_ord_smooth(NPROC,NGLOB_AB,&
            ddat_glob,buffer_recv_vector_ext_mesh_smooth,num_interfaces_ext_mesh,&
            max_nibool_interfaces_ext_mesh, &
            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
            request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
            my_neighbors_ext_mesh,myrank)
        endif
        !!!!!!!!!!!!!!!!!
      enddo !iphase = 1,2

      ddat_glob(:) = ddat_glob(:) * rvol(:)
      !! update
      dat_glob(:) = dat_glob(:) + ddat_glob(:)
      !! info
      ! if (mod(PRINT_INFO_PER_STEP, istep) == 0) then
      !   if (myrank == 0) print *, 'Step:', istep
      !   min_val = minval(dat_glob)
      !   max_val = maxval(dat_glob)
      !   call min_all_cr(min_val, min_val_glob)
      !   call max_all_cr(max_val, max_val_glob)
      !   if (myrank == 0) then
      !     print *, '    minval:', min_val_glob
      !     print *, '    maxval:', max_val_glob
      !     call cpu_time(tnow)
      !     print *, 'time since last message:', tnow-tlast
      !     call cpu_time(tlast)
      !   endif
      ! endif
      !!!!!!!!!!!!!
      !! broadcast glob array back to local array
      do ispec = 1, NSPEC_AB
        do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          if ( is_CPML(ispec) ) then 
            dat_glob(iglob) = dat_bak(i,j,k,ispec)
          end if
        enddo;enddo;enddo
      enddo

      !dat = reshape(dat_glob(reshape(ibool,[5*5*5*NSPEC_AB])), [5,5,5,NSPEC_AB])
      do ispec = 1, NSPEC_AB
        do k = 1,NGLLZ;do j = 1,NGLLY;do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dat(i,j,k,ispec) = dat_glob(iglob)
        enddo;enddo;enddo
      enddo
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call synchronize_all()
    enddo
  else 
    call smooth_pde_cuda(is_sph_int,NSPEC_AB,NGLOB_AB,NPROC,ntstep, &
                          size(phase_ispec_inner_elastic,1),phase_ispec_inner_elastic, &
                          nspec_outer_elastic,nspec_inner_elastic,xixstore, &
                          xiystore,xizstore,etaxstore,etaystore,etazstore, &
                          gammaxstore,gammaystore,gammazstore, &
                          jacobian_all,rotate_r,wgllwgll_xy,wgllwgll_xz, &
                          wgllwgll_yz,hprime_xxT,hprimewgll_xx,ibool, &
                          is_CPML,cv,ch,num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                          my_neighbors_ext_mesh,nibool_interfaces_ext_mesh, &
                          ibool_interfaces_ext_mesh,dat_bak,dat,rvol,dat_glob,ddat_glob, &
                          irregular_element_number,xix_regular,NSPEC_IRREGULAR)
  endif
  call synchronize_all()
  ! call cpu_time(t2)
  ! if (myrank == 0) & 
    ! print *, 'Computation time with PDE-based smoothing on CPU:', t2-t1
  !! output
  ! file output
  ! smoothed kernel file name
  ! statistics
  ! min_val = minval(dat_glob)
  ! max_val = maxval(dat_glob)
  ! call min_all_cr(min_val, min_val_glob)
  ! call max_all_cr(max_val, max_val_glob)
  ! if (myrank == 0) then
  !   print *, 'After smoothing:'
  !   ! print *, '  '//trim(kernel_name)
  !   print *, '    minval:', min_val_glob
  !   print *, '    maxval:', max_val_glob
  ! endif
  ! write(ks_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_name)//'_smooth.bin'
  ! open(IOUT,file=trim(ks_file),status='unknown',form='unformatted',iostat=ier)
  ! if (ier /= 0) stop 'Error opening smoothed kernel file'
  ! write(IOUT) dat
  ! close(IOUT)
  ! if (myrank == 0) print *,'written: ',trim(ks_file)

  ! deallocate(ibool,irregular_element_number)
  ! deallocate(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian)
  ! deallocate(rotate_r)
  ! deallocate(dat, dat_glob, ddat_glob)
  ! deallocate(buffer_send_vector_ext_mesh_smooth, &
  !            buffer_recv_vector_ext_mesh_smooth)
  ! deallocate(rvol, rvol_local)

  ! call finalize_mpi()

end subroutine smooth_sem_pde

subroutine assemble_MPI_send_smooth(NPROC,NGLOB_AB, &
          array_val,buffer_send_vector_ext_mesh_smooth, &
          buffer_recv_vector_ext_mesh_smooth, &
          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
          my_neighbors_ext_mesh, &
          request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

    ! sends data

  use constants, only: CUSTOM_REAL, itag

  implicit none

  integer :: NPROC
  integer :: NGLOB_AB

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), &
    dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_vector_ext_mesh_smooth,buffer_recv_vector_ext_mesh_smooth

  integer, dimension(num_interfaces_ext_mesh) :: &
    nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh):: &
    ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: &
    request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer ipoin,iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        buffer_send_vector_ext_mesh_smooth(ipoin,iinterface) = &
             array_val(ibool_interfaces_ext_mesh(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      call isend_cr(buffer_send_vector_ext_mesh_smooth(1,iinterface), &
                    nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_send_vector_ext_mesh(iinterface))
      call irecv_cr(buffer_recv_vector_ext_mesh_smooth(1,iinterface), &
                    nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_recv_vector_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_send_smooth


  subroutine assemble_MPI_w_ord_smooth(NPROC,NGLOB_AB, &
          array_val,buffer_recv_vector_ext_mesh_smooth,num_interfaces_ext_mesh, &
          max_nibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
          request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
          my_neighbors_ext_mesh,myrank)

! waits for data to receive and assembles

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: NPROC
  integer :: NGLOB_AB
! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,myrank

  real(kind=CUSTOM_REAL), &
    dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh_smooth

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh)::&
    ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: &
    request_send_vector_ext_mesh,request_recv_vector_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbors_ext_mesh

  real(kind=CUSTOM_REAL), &
    dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
    mybuffer
  integer :: ipoin,iinterface,iglob
  logical :: need_add_my_contrib

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if (NPROC == 1) return

! move interface values of array_val to local buffers
  do iinterface = 1, num_interfaces_ext_mesh
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
      mybuffer(ipoin,iinterface) = array_val(iglob)
     ! set them to zero right away to avoid counting it more than once during
     ! assembly:
     ! buffers of higher rank get zeros on nodes shared with current buffer
      array_val(iglob) = 0._CUSTOM_REAL
    enddo
  enddo

! wait for communications completion (recv)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_recv_vector_ext_mesh(iinterface))
  enddo

! adding all contributions in order of processor rank
  need_add_my_contrib = .true.
  do iinterface = 1, num_interfaces_ext_mesh
    if (need_add_my_contrib .and. myrank < my_neighbors_ext_mesh(iinterface)) &
      call add_my_contrib()
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
      array_val(iglob) = array_val(iglob) + &
        buffer_recv_vector_ext_mesh_smooth(ipoin,iinterface)
    enddo
  enddo
  if (need_add_my_contrib) call add_my_contrib()

! wait for communications completion (send)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_send_vector_ext_mesh(iinterface))
  enddo

  contains

    subroutine add_my_contrib()

    integer :: my_iinterface,my_ipoin

    do my_iinterface = 1, num_interfaces_ext_mesh
      do my_ipoin = 1, nibool_interfaces_ext_mesh(my_iinterface)
        iglob = ibool_interfaces_ext_mesh(my_ipoin,my_iinterface)
        array_val(iglob) = array_val(iglob) + &
          mybuffer(my_ipoin,my_iinterface)
      enddo
    enddo
    need_add_my_contrib = .false.

    end subroutine add_my_contrib

  end subroutine assemble_MPI_w_ord_smooth


  subroutine get_gradient_element(s, dx_elem, dy_elem, dz_elem, &
     xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
  use constants, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ
  use specfem_par, only: hprime_xxT,hprime_yyT,hprime_zzT
  implicit none
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ), intent(in) :: s
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ), intent(out) :: &
          dx_elem, dy_elem, dz_elem
  integer :: i,j,k,l
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
  real(kind=CUSTOM_REAL) :: temp1l, temp2l, temp3l
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xix,xiy,xiz,etax,&
          etay,etaz,gammax,gammay,gammaz
  dx_elem(:,:,:) = 0.0
  dy_elem(:,:,:) = 0.0
  dz_elem(:,:,:) = 0.0
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        temp1l = 0.0
        temp2l = 0.0
        temp3l = 0.0
        do l = 1, NGLLX
          hp1 = hprime_xxT(l,i)
          temp1l = temp1l + s(l,j,k) * hp1
        enddo
        do l = 1, NGLLY
          hp2 = hprime_yyT(l,j)
          temp2l = temp2l + s(i,l,k) * hp2
        enddo
        do l = 1, NGLLZ
          hp3 = hprime_zzT(l,k)
          temp3l = temp3l + s(i,j,l) * hp3
        enddo
        dx_elem(i,j,k)=temp1l*xix(i,j,k)+&
                temp2l*etax(i,j,k)+temp3l*gammax(i,j,k)
        dy_elem(i,j,k)=temp1l*xiy(i,j,k)+&
                temp2l*etay(i,j,k)+temp3l*gammay(i,j,k)
        dz_elem(i,j,k)=temp1l*xiz(i,j,k)+&
                temp2l*etaz(i,j,k)+temp3l*gammaz(i,j,k)
      enddo
    enddo
  enddo
  end subroutine get_gradient_element

  subroutine get_gradient_element_regular(s, xix, dx_elem, dy_elem, dz_elem)
    use constants, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ
    use specfem_par, only: hprime_xxT,hprime_yyT,hprime_zzT
    implicit none
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ), intent(in) :: s
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ), intent(out) :: &
            dx_elem, dy_elem, dz_elem
    real(kind=CUSTOM_REAL), intent(in) :: xix
    integer :: i,j,k,l
    real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
    real(kind=CUSTOM_REAL) :: temp1l, temp2l, temp3l
    dx_elem(:,:,:) = 0.0
    dy_elem(:,:,:) = 0.0
    dz_elem(:,:,:) = 0.0
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          temp1l = 0.0
          temp2l = 0.0
          temp3l = 0.0
          do l = 1, NGLLX
            hp1 = hprime_xxT(l,i)
            temp1l = temp1l + s(l,j,k) * hp1
          enddo
          do l = 1, NGLLY
            hp2 = hprime_yyT(l,j)
            temp2l = temp2l + s(i,l,k) * hp2
          enddo
          do l = 1, NGLLZ
            hp3 = hprime_zzT(l,k)
            temp3l = temp3l + s(i,j,l) * hp3
          enddo
          dx_elem(i,j,k)=temp1l*xix+&
                  temp2l*xix+temp3l*xix
          dy_elem(i,j,k)=temp1l*xix+&
                  temp2l*xix+temp3l*xix
          dz_elem(i,j,k)=temp1l*xix+&
                  temp2l*xix+temp3l*xix
        enddo
      enddo
    enddo
  end subroutine get_gradient_element_regular
end module
