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


subroutine read_kernels_iso()

! reads in smoothed kernels: bulk, beta, rho

  use tomography_kernels_iso
  ! use fullwave_adjoint_tomo_par, only: USE_RHO_SCALING_FWAT
  use fwat_input, only: tomo_par


  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: min_vp,min_vs,max_vp,max_vs,min_rho,max_rho
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'reading kernels...'

  ! allocate arrays for storing kernels and perturbations
  ! isotropic arrays
  allocate(kernel_bulk(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_beta(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_rho(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'error allocating kernel arrays'

  ! initializes arrays
  kernel_bulk = 0.0_CUSTOM_REAL
  kernel_beta = 0.0_CUSTOM_REAL
  kernel_rho = 0.0_CUSTOM_REAL

  ! reads in smoothed (& summed) event kernel
  if (USE_ALPHA_BETA_RHO) then
    ! reads in alpha kernel
    fname = 'alpha_kernel_smooth'
  else
    ! reads in bulk_c kernel
    fname = 'bulk_c_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_bulk(:,:,:,1:nspec)
  close(IIN)

  ! beta kernel
  if (USE_ALPHA_BETA_RHO) then
    ! reads in beta kernel
    fname = 'beta_kernel_smooth'
  else
    ! reads in bulk_beta kernel
    fname = 'bulk_beta_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_beta(:,:,:,1:nspec)
  close(IIN)

  ! rho kernel
  if (tomo_par%USE_RHO_SCALING_FWAT) then
    ! uses scaling relation with shear perturbations
    kernel_rho(:,:,:,:) = RHO_SCALING * kernel_beta(:,:,:,:)
    if (myrank == 0) print *,'  rho kernel uses scaling with shear kernel: scaling value = ',RHO_SCALING
  else
    ! uses rho kernel
    fname = 'rhop_kernel_smooth'
    write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) kernel_rho(:,:,:,1:nspec)
    close(IIN)
  endif

  ! statistics
  call min_all_cr(minval(kernel_bulk),min_vp)
  call max_all_cr(maxval(kernel_bulk),max_vp)

  call min_all_cr(minval(kernel_beta),min_vs)
  call max_all_cr(maxval(kernel_beta),max_vs)

  call min_all_cr(minval(kernel_rho),min_rho)
  call max_all_cr(maxval(kernel_rho),max_rho)

  if (myrank == 0) then
    print *
    print *,'initial kernels:'
    if (USE_ALPHA_BETA_RHO) then
      print *,'  alpha min/max    : ',min_vp,max_vp
      print *,'  beta min/max     : ',min_vs,max_vs
    else
      print *,'  bulk_c min/max   : ',min_vp,max_vp
      print *,'  bulk_beta min/max: ',min_vs,max_vs
    endif
    print *,'  rho min/max      : ',min_rho,max_rho
    print *
  endif
  call synchronize_all()

  ! statistics output
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_kernels_minmax',status='unknown')
    write(IOUT,*) '#min_vs #max_vs #min_vp #max_vp #min_rho #max_rho'
    write(IOUT,'(4e24.12)') min_vs, max_vs, min_vp, max_vp, min_rho, max_rho
    close(IOUT)
  endif

end subroutine read_kernels_iso

!
!-------------------------------------------------------------------------------------------------
!


subroutine read_kernels_tiso()

! reads in smoothed kernels: bulk, betav, betah, eta

  use tomography_kernels_tiso
  !use tomography_par, only: USE_ALPHA_BETA_RHO_TISO

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: min_vsv,min_vsh,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk
  real(kind=CUSTOM_REAL) :: min_vpv,min_vph,max_vpv,max_vph
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'reading kernels...'

  ! allocate arrays for storing kernels and perturbations
  ! transversely isotropic arrays
  if (USE_ALPHA_BETA_RHO_TISO) then 
     allocate(kernel_alphav(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_alphah(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_betav(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_betah(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_eta(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  else

      allocate(kernel_bulk(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_betav(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_betah(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_eta(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  endif
  if (ier /= 0) stop 'error allocating kernel arrays'

  ! initializes arrays
  if (USE_ALPHA_BETA_RHO_TISO) then
  kernel_alphav = 0.0_CUSTOM_REAL
  kernel_alphah = 0.0_CUSTOM_REAL
  else
  kernel_bulk = 0.0_CUSTOM_REAL
  endif
  kernel_betav = 0.0_CUSTOM_REAL
  kernel_betah = 0.0_CUSTOM_REAL
  kernel_eta = 0.0_CUSTOM_REAL

  if (USE_ALPHA_BETA_RHO_TISO) then
  ! alphav kernel
  fname = 'alphav_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_alphav(:,:,:,1:nspec)
  close(IIN)
  ! alphah kernel
  fname = 'alphah_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_alphah(:,:,:,1:nspec)
  close(IIN)
  ! betav kernel
  fname = 'betav_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betav(:,:,:,1:nspec)
  close(IIN)
  ! betah kernel
  fname = 'betah_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betah(:,:,:,1:nspec)
  close(IIN)

  else
  ! bulk kernel
  fname = 'bulk_c_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_bulk(:,:,:,1:nspec)
  close(IIN)
 
  ! bulk betav kernel
  fname = 'bulk_betav_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betav(:,:,:,1:nspec)
  close(IIN)

  ! bulk betah kernel
  fname = 'bulk_betah_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betah(:,:,:,1:nspec)
  close(IIN)
  
  endif

  ! eta kernel
  fname = 'eta_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(INPUT_KERNELS_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(INPUT_KERNELS_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_eta(:,:,:,1:nspec)
  close(IIN)


  ! statistics
  if (USE_ALPHA_BETA_RHO_TISO) then
  call min_all_cr(minval(kernel_alphah),min_vph)
  call max_all_cr(maxval(kernel_alphah),max_vph)

  call min_all_cr(minval(kernel_alphav),min_vpv)
  call max_all_cr(maxval(kernel_alphav),max_vpv)

  else 
  call min_all_cr(minval(kernel_bulk),min_bulk)
  call max_all_cr(maxval(kernel_bulk),max_bulk)
  endif
  
  call min_all_cr(minval(kernel_betah),min_vsh)
  call max_all_cr(maxval(kernel_betah),max_vsh)

  call min_all_cr(minval(kernel_betav),min_vsv)
  call max_all_cr(maxval(kernel_betav),max_vsv)

  call min_all_cr(minval(kernel_eta),min_eta)
  call max_all_cr(maxval(kernel_eta),max_eta)

  if (myrank == 0) then
    print *
    print *,'initial kernels:'
    if (USE_ALPHA_BETA_RHO_TISO) then
    print *,'  alphav min/max: ',min_vpv,max_vpv
    print *,'  alphah min/max: ',min_vph,max_vph
    else
    print *,'  bulk min/max : ',min_bulk,max_bulk
    endif
    print *,'  betav min/max: ',min_vsv,max_vsv
    print *,'  betah min/max: ',min_vsh,max_vsh
    print *,'  eta min/max  : ',min_eta,max_eta
    print *
  endif
  call synchronize_all()

  ! statistics output
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_kernels_minmax',status='unknown')
    if (USE_ALPHA_BETA_RHO_TISO) then
    write(IOUT,*) '#min_vpv #max_vpv $min_vph #max_vph #min_vsv #max_vsv #min_vsh #max_vsh #min_eta #max_eta'
    write(IOUT,'(4e24.12)') min_vpv, max_vpv, min_vph, max_vph, min_vsv, max_vsv, min_vsh, max_vsh, min_eta, max_eta
    else 
    write(IOUT,*) '#min_bulk #max_bulk #min_vsv #max_vsv #min_vsh #max_vsh #min_eta #max_eta'
    write(IOUT,'(4e24.12)') min_bulk, max_bulk, min_vsv, max_vsv, min_vsh, max_vsh, min_eta, max_eta
    endif
    close(IOUT)
  endif

end subroutine read_kernels_tiso

