!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
!
! United States and French Government Sponsorship Acknowledged.


!==============================================================================
! \file save_adjoint_kernels
!
! TODO
! * Better doxygen documentation.
!==============================================================================


!> Save kernels.

  subroutine save_adjoint_kernels_fwat()

  use constants, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ
  use specfem_par, only: LOCAL_PATH, myrank, sigma_kl, NSPEC_AB, ADIOS_FOR_KERNELS, NOISE_TOMOGRAPHY, NSPEC_ADJOINT, &
                         APPROXIMATE_HESS_KL, ANISOTROPIC_KL, SAVE_TRANSVERSE_KL, SAVE_AZIMUTH_KL

  use specfem_par_acoustic, only: ACOUSTIC_SIMULATION
  use specfem_par_elastic, only: ELASTIC_SIMULATION
  use specfem_par_poroelastic, only: POROELASTIC_SIMULATION

  implicit none

  interface
    subroutine save_kernels_elastic_fwat(adios_handle, alphav_kl, alphah_kl, &
                                    betav_kl, betah_kl, eta_kl, &
                                    rhop_kl, alpha_kl, beta_kl, &
                                    gcp_kl, gsp_kl)

      use constants, only: CUSTOM_REAL

      integer(kind=8) :: adios_handle
      ! FIXME
      ! Break the CUSTOM_REAL stuff.
      ! put all this file in a module so interface is implicit
      ! OR
      ! redo what was done before SVN revision 22718
      !
      ! see other FIXME below (same than see one)
!! DK DK: sorry, we cannot afford to break the code; too many people use it; I thus put CUSTOM_REAL back
      real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
          alphav_kl,alphah_kl,betav_kl,betah_kl, &
          eta_kl, rhop_kl, alpha_kl, beta_kl, &
          gcp_kl, gsp_kl    !@ WK @!
    end subroutine save_kernels_elastic_fwat
  end interface

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: alphav_kl, &
                                                            alphah_kl, &
                                                            betav_kl, &
                                                            betah_kl, &
                                                            eta_kl

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: rhop_kl, &
                                                            alpha_kl, &
                                                            beta_kl
  !@ WK @!
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: gcp_kl, &
                                                            gsp_kl
  !@ WK @!

  integer(kind=8) :: adios_handle
  integer :: ier

  ! flag to save GLL weights
  logical,parameter :: SAVE_WEIGHTS = .false.
  
  if (ADIOS_FOR_KERNELS) then
    call define_kernel_adios_variables(adios_handle, SAVE_WEIGHTS)
  endif

  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    call save_kernels_acoustic(adios_handle)
  endif

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    ! allocates temporary transversely isotropic kernels
    if (ANISOTROPIC_KL) then
      if (SAVE_TRANSVERSE_KL) then
        allocate(alphav_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 alphah_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 betav_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 betah_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 eta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 stat=ier)
        if (ier /= 0) stop 'error allocating arrays alphav_kl,...'

        ! derived kernels
        ! vp kernel
        allocate(alpha_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
        if (ier /= 0) stop 'error allocating array alpha_kl'
        ! vs kernel
        allocate(beta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
        if (ier /= 0) stop 'error allocating array beta_kl'
      else if (SAVE_AZIMUTH_KL) then !@ WK @!
        allocate(betav_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 betah_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 gcp_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 gsp_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
                 stat=ier)
         if (ier /= 0) stop 'error allocating arrays betav_kl,...'

        ! derived kernels
        ! vs kernel
        allocate(beta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
        if (ier /= 0) stop 'error allocating array beta_kl'
        !@ WK @! 
      endif
    else
      ! derived kernels
      ! vp kernel
      allocate(alpha_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if (ier /= 0) stop 'error allocating array alpha_kl'
      ! vs kernel
      allocate(beta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if (ier /= 0) stop 'error allocating array beta_kl'
      ! density prime kernel
      allocate(rhop_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if (ier /= 0) stop 'error allocating array rhop_kl'
    endif

    call save_kernels_elastic_fwat(adios_handle, alphav_kl, alphah_kl, &
                              betav_kl, betah_kl, eta_kl, &
                              rhop_kl, alpha_kl, beta_kl, &
                              gcp_kl, gsp_kl)
  endif

  if (POROELASTIC_SIMULATION) then
    call save_kernels_poroelastic(adios_handle)
  endif

  ! save weights for volume integration,
  ! in order to benchmark the kernels with analytical expressions
  if (SAVE_WEIGHTS) then
    call save_weights_kernel()
  endif

  ! for noise simulations --- noise strength kernel
  if (NOISE_TOMOGRAPHY == 3) then
    call save_kernels_strength_noise(myrank,LOCAL_PATH,sigma_kl,NSPEC_AB)
  endif

  ! for preconditioner
  if (APPROXIMATE_HESS_KL) then
    call save_kernels_Hessian_fwat(adios_handle)
  endif

  if (ADIOS_FOR_KERNELS) then
    call perform_write_adios_kernels(adios_handle)
  endif

  if (ELASTIC_SIMULATION) then
    ! frees temporary arrays
    if (ANISOTROPIC_KL) then
      if (SAVE_TRANSVERSE_KL) then
        deallocate(alphav_kl,alphah_kl,betav_kl,betah_kl,eta_kl)
        deallocate(alpha_kl,beta_kl)
      else if (SAVE_AZIMUTH_KL) then !@ WK @!
        deallocate(betav_kl,betah_kl,beta_kl,gcp_kl,gsp_kl)
      endif
    else
      deallocate(rhop_kl,alpha_kl,beta_kl)
    endif
  endif

  end subroutine save_adjoint_kernels_fwat

!
!-------------------------------------------------------------------------------------------------
!

!> Save elastic related kernels

  subroutine save_kernels_elastic_fwat(adios_handle, alphav_kl, alphah_kl, &
                                betav_kl, betah_kl, eta_kl, &
                                rhop_kl, alpha_kl, beta_kl, &
                                gcp_kl, gsp_kl) !@ WK @!

  use specfem_par, only: CUSTOM_REAL,NSPEC_AB,ibool,mustore,kappastore,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL,FOUR_THIRDS, &
                         ADIOS_FOR_KERNELS,IOUT,prname,SAVE_MOHO_MESH,SAVE_AZIMUTH_KL  !@ WK @!
  use specfem_par_elastic

  implicit none

  interface
    subroutine save_kernels_elastic_adios(adios_handle, alphav_kl, alphah_kl, &
                                          betav_kl, betah_kl, eta_kl, &
                                          rhop_kl, alpha_kl, beta_kl)

      use constants, only: CUSTOM_REAL

      integer(kind=8), intent(in) :: adios_handle
      ! FIXME
      ! see other FIXME above.
!! DK DK: sorry, we cannot afford to break the code; too many people use it; I thus put CUSTOM_REAL back
      real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
          alphav_kl,alphah_kl,betav_kl,betah_kl, &
          eta_kl, rhop_kl, alpha_kl, beta_kl
    end subroutine save_kernels_elastic_adios
  end interface

  integer(kind=8) :: adios_handle
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    alphav_kl,alphah_kl,betav_kl,betah_kl, &
    eta_kl, rhop_kl, alpha_kl, beta_kl, &
    gcp_kl, gsp_kl

  ! local parameters
  integer:: ispec,i,j,k,iglob,ier
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal
  !!! WK added to avoid overwrite rho_kl mu_kl kappa_kl
  real(kind=CUSTOM_REAL) :: rho_kl_local,mu_kl_local,kappa_kl_local

  ! Transverse isotropic paramters
  real(kind=CUSTOM_REAL) :: A,N,C,L,F,eta
  real(kind=CUSTOM_REAL), dimension(21) :: cijkl_kl_local
  real(kind=CUSTOM_REAL), dimension(7) :: an_kl


  ! finalizes calculation of rhop, beta, alpha kernels
  do ispec = 1, NSPEC_AB

    ! elastic simulations
    if (ispec_is_elastic(ispec)) then

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)

            ! Store local material values
            rhol = rho_vs(i,j,k,ispec)*rho_vs(i,j,k,ispec) / mustore(i,j,k,ispec)
            mul = mustore(i,j,k,ispec)
            kappal = kappastore(i,j,k,ispec)

            if (ANISOTROPIC_KL) then
              if (SAVE_TRANSVERSE_KL) then
                cijkl_kl_local(:) = - cijkl_kl(:,i,j,k,ispec)

                ! Computes parameters for an isotropic model
                A = kappal + FOUR_THIRDS * mul
                C = A
                L = mul
                N = mul
                F = kappal - 2._CUSTOM_REAL/3._CUSTOM_REAL * mul
                eta = 1._CUSTOM_REAL

                ! note: cijkl_kl_local() is fully anisotropic C_ij kernel components (non-dimensionalized)
                !          for GLL point at (i,j,k,ispec)

                ! Purpose : compute the kernels for the An coeffs (an_kl)
                ! from the kernels for Cij (cijkl_kl_local)

                ! Definition of the input array cij_kl :
                ! cij_kl(1) = C11 ; cij_kl(2) = C12 ; cij_kl(3) = C13
                ! cij_kl(4) = C14 ; cij_kl(5) = C15 ; cij_kl(6) = C16
                ! cij_kl(7) = C22 ; cij_kl(8) = C23 ; cij_kl(9) = C24
                ! cij_kl(10) = C25 ; cij_kl(11) = C26 ; cij_kl(12) = C33
                ! cij_kl(13) = C34 ; cij_kl(14) = C35 ; cij_kl(15) = C36
                ! cij_kl(16) = C44 ; cij_kl(17) = C45 ; cij_kl(18) = C46
                ! cij_kl(19) = C55 ; cij_kl(20) = C56 ; cij_kl(21) = C66
                ! where the Cij (Voigt's notation) are defined as function of
                ! the components of the elastic tensor in spherical coordinates
                ! by eq. (A.1) of Chen & Tromp, GJI 168 (2007)

                ! From the relations giving Cij in function of An
                ! Checked with Min Chen's results (routine build_cij)

                an_kl(1) = cijkl_kl_local(1)+cijkl_kl_local(2)+cijkl_kl_local(7)  !A
                an_kl(2) = cijkl_kl_local(12)                                     !C
                an_kl(3) = -2*cijkl_kl_local(2)+cijkl_kl_local(21)                !N
                an_kl(4) = cijkl_kl_local(16)+cijkl_kl_local(19)                  !L
                an_kl(5) = cijkl_kl_local(3)+cijkl_kl_local(8)                    !F

                ! for parameterization: ( alpha_v, alpha_h, beta_v, beta_h, eta, rho )
                ! K_alpha_v
                alphav_kl(i,j,k,ispec) = 2.0 * C * an_kl(2)
                ! K_alpha_h
                alphah_kl(i,j,k,ispec) = 2.0 * A * an_kl(1) + 2.0 * A * eta * an_kl(5)
                ! K_beta_v
                betav_kl(i,j,k,ispec) = 2.0 * L * an_kl(4) - 4.0 * L * eta * an_kl(5)
                ! K_beta_h
                betah_kl(i,j,k,ispec) = 2.0 * N * an_kl(3)
                ! K_eta
                eta_kl(i,j,k,ispec) = F * an_kl(5)

                ! to check: isotropic kernels from transverse isotropic ones
                alpha_kl(i,j,k,ispec) = alphav_kl(i,j,k,ispec) &
                                                  + alphah_kl(i,j,k,ispec)
                beta_kl(i,j,k,ispec) = betav_kl(i,j,k,ispec) &
                                                  + betah_kl(i,j,k,ispec)
              else if (SAVE_AZIMUTH_KL) then
                cijkl_kl_local(:) = - cijkl_kl(:,i,j,k,ispec)

                ! Computes parameters for an isotropic model
                L = mul
                N = mul
                eta = 1._CUSTOM_REAL

                
                an_kl(1) = cijkl_kl_local(1)+cijkl_kl_local(2)+cijkl_kl_local(7)  !A
                an_kl(2) = cijkl_kl_local(12)                                     !C
                an_kl(3) = -2*cijkl_kl_local(2)+cijkl_kl_local(21)                !N
                an_kl(4) = cijkl_kl_local(16)+cijkl_kl_local(19)                  !L
                an_kl(5) = cijkl_kl_local(3)+cijkl_kl_local(8)                    !F
                an_kl(6) = -cijkl_kl_local(16)+cijkl_kl_local(19)                 !Gc
                an_kl(7) = -cijkl_kl_local(17)                                    !Gs

                ! for parameterization: ( alpha_v, alpha_h, beta_v, beta_h, eta, rho )
                ! K_beta_v
                betav_kl(i,j,k,ispec) = 2.0 * L * an_kl(4) - 4.0 * L * eta * an_kl(5)
                ! K_beta_h
                betah_kl(i,j,k,ispec) = 2.0 * N * an_kl(3)
                ! K_Gcp
                gcp_kl(i,j,k,ispec) = mul*an_kl(6)
                ! K_Gsp
                gsp_kl(i,j,k,ispec) = mul*an_kl(7)

                ! to check: isotropic kernels from transverse isotropic ones
                beta_kl(i,j,k,ispec) = betav_kl(i,j,k,ispec) &
                                                  + betah_kl(i,j,k,ispec)

              endif ! SAVE_TRANSVERSE_KL

            else

              ! isotropic kernels

              ! isotropic adjoint kernels (see e.g. Tromp et al. 2005)
              ! for a parameterization: (rho,mu,kappa) "primary" kernels
              ! density kernel
              ! multiplies with rho
              rho_kl_local = - rhol * rho_kl(i,j,k,ispec)

              ! shear modulus kernel
              mu_kl_local = - 2._CUSTOM_REAL * mul * mu_kl(i,j,k,ispec)

              ! bulk modulus kernel
              kappa_kl_local = - kappal * kappa_kl(i,j,k,ispec)

              ! for a parameterization: (rho,alpha,beta)
              ! density prime kernel
              rhop_kl(i,j,k,ispec) = rho_kl_local + kappa_kl_local + mu_kl_local
              ! if (isnan(rhop_kl(i,j,k,ispec))) rhop_kl(i,j,k,ispec) = 0._CUSTOM_REAL

              ! vs kernel
              beta_kl(i,j,k,ispec) = 2._CUSTOM_REAL * (mu_kl_local &
                    - 4._CUSTOM_REAL * mul / (3._CUSTOM_REAL * kappal) * kappa_kl_local)
              ! if (isnan(beta_kl(i,j,k,ispec))) beta_kl(i,j,k,ispec) = 0._CUSTOM_REAL

              ! vp kernel
              alpha_kl(i,j,k,ispec) = 2._CUSTOM_REAL * (1._CUSTOM_REAL &
                    + 4._CUSTOM_REAL * mul / (3._CUSTOM_REAL * kappal) ) * kappa_kl_local
              ! if (isnan(alpha_kl(i,j,k,ispec))) alpha_kl(i,j,k,ispec) = 0._CUSTOM_REAL

            endif

          enddo
        enddo
      enddo

    endif ! elastic

  enddo

  if (ADIOS_FOR_KERNELS) then
    call save_kernels_elastic_adios(adios_handle, alphav_kl, alphah_kl, &
                                      betav_kl, betah_kl, eta_kl, &
                                      rhop_kl, alpha_kl, beta_kl)
  else
    if (ANISOTROPIC_KL) then

      ! outputs transverse isotropic kernels only
      if (SAVE_TRANSVERSE_KL) then
        ! transverse isotropic kernels
        ! (alpha_v, alpha_h, beta_v, beta_h, eta, rho ) parameterization
        open(unit=IOUT,file=trim(prname)//'alphav_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) alphav_kl
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'alphah_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) alphah_kl
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'betav_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) betav_kl
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'betah_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) betah_kl
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'eta_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) eta_kl
        close(IOUT)

        ! transverse isotropic test kernels
        open(unit=IOUT,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT)  alpha_kl
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT)  beta_kl
        close(IOUT)

      else if (SAVE_AZIMUTH_KL) then
        open(unit=IOUT,file=trim(prname)//'betav_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) betav_kl
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'betah_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) betah_kl
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'gcp_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) gcp_kl
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'gsp_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) gsp_kl
        close(IOUT)

        ! transverse isotropic test kernels
        open(unit=IOUT,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT)  beta_kl
        close(IOUT)
       

      else
        ! fully anisotropic kernels
        ! note: the C_ij and density kernels are not for relative perturbations (delta ln( m_i) = delta m_i / m_i),
        !          but absolute perturbations (delta m_i = m_i - m_0).
        ! Kappa and mu are for absolute perturbations, can be used to check with purely isotropic versions.
        open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT)  - rho_kl
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'cijkl_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) - cijkl_kl
        close(IOUT)

      endif

    else

      ! save kernels to binary files
      open(unit=IOUT,file=prname(1:len_trim(prname))//'rho_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file rho_kernel.bin'
      write(IOUT) rho_kl
      close(IOUT)

      open(unit=IOUT,file=prname(1:len_trim(prname))//'mu_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file mu_kernel.bin'
      write(IOUT) mu_kl
      close(IOUT)

      open(unit=IOUT,file=prname(1:len_trim(prname))//'kappa_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file kappa_kernel.bin'
      write(IOUT) kappa_kl
      close(IOUT)

      open(unit=IOUT,file=prname(1:len_trim(prname))//'rhop_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file rhop_kernel.bin'
      write(IOUT) rhop_kl
      close(IOUT)

      open(unit=IOUT,file=prname(1:len_trim(prname))//'beta_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file beta_kernel.bin'
      write(IOUT) beta_kl
      close(IOUT)

      open(unit=IOUT,file=prname(1:len_trim(prname))//'alpha_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file alpha_kernel.bin'
      write(IOUT) alpha_kl
      close(IOUT)
    endif

    if (SAVE_MOHO_MESH) then
      open(unit=IOUT,file=prname(1:len_trim(prname))//'moho_kernel.bin',status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) stop 'error opening file moho_kernel.bin'
      write(IOUT) moho_kl
      close(IOUT)
    endif
  endif

  end subroutine save_kernels_elastic_fwat

!
!-------------------------------------------------------------------------------------------------
!

!> Save Hessians

  subroutine save_kernels_Hessian_fwat(adios_handle)

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none

  integer(kind=8) :: adios_handle
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: hess_kl_local, &
                                                            hess_ac_kl_local
 

  integer :: ier

  allocate(hess_kl_local(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
  if (ier /= 0) stop 'error allocating array hess_kl_local'
  allocate(hess_ac_kl_local(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
  if (ier /= 0) stop 'error allocating array hess_ac_kl_local'


  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    ! scales approximate Hessian
    hess_ac_kl_local(:,:,:,:) = 2._CUSTOM_REAL * hess_ac_kl(:,:,:,:)
  endif

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    ! scales approximate Hessian
    hess_kl_local(:,:,:,:) = 2._CUSTOM_REAL * hess_kl(:,:,:,:)
  endif

  if (ADIOS_FOR_KERNELS) then
    call save_kernels_Hessian_adios(adios_handle)
  else
    ! acoustic domains
    if (ACOUSTIC_SIMULATION) then
      ! stores into file
      open(unit=IOUT,file=trim(prname)//'hess_acoustic_kernel.bin', &
           status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0) stop 'error opening file hess_acoustic_kernel.bin'
      write(IOUT) hess_ac_kl_local
      close(IOUT)
    endif

    ! elastic domains
    if (ELASTIC_SIMULATION) then
      ! stores into file
      open(unit=IOUT,file=trim(prname)//'hess_kernel.bin', &
           status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0) stop 'error opening file hess_kernel.bin'
      write(IOUT) hess_kl_local
      close(IOUT)
    endif
  endif
  deallocate(hess_kl_local,hess_ac_kl_local)

  end subroutine save_kernels_Hessian_fwat


