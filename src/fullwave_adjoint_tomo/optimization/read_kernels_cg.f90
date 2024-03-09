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
subroutine read_kernels_cg_iso_old()

! reads in smoothed kernels from former iteration in OUTPUT_SUM.old/ : bulk, betav, betah, eta

  use tomography_kernels_iso_cg
  ! use fullwave_adjoint_tomo_par, only: USE_RHO_SCALING_FWAT
  use fwat_input, only: tomo_par


  implicit none
  real(kind=CUSTOM_REAL) :: min_vp,max_vp,min_vs,max_vs,min_rho,max_rho
  logical:: exist,exist_all,use_old_gradient_all
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'reading old kernels...'

  ! allocate arrays for storing kernels and perturbations
  allocate(kernel_bulk_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_beta_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_rho_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)

  if (ier /= 0) stop 'error allocating kernel arrays'

  ! initializes arrays
  kernel_bulk_old = 0.0_CUSTOM_REAL
  kernel_beta_old = 0.0_CUSTOM_REAL
  kernel_rho_old = 0.0_CUSTOM_REAL

  ! checks if files are available:
  if (USE_ALPHA_BETA_RHO) then
    ! alpha kernel
    fname = 'alpha_kernel_smooth'
  else
    fname = 'bulk_c_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if (.not. exist) then
    print *,'Error file does not exist: ',trim(m_file)
    call exit_mpi(myrank,'file not exist')
  endif
  ! makes sure all processes have same flag
  call any_all_l(exist,exist_all)
  if (.not. exist_all) then
    print *,'old kernels do not exist: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif


  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_bulk_old(:,:,:,1:nspec)
  close(IIN)

  ! beta kernel
  if (USE_ALPHA_BETA_RHO) then
    fname = 'beta_kernel_smooth'
  else
    fname = 'bulk_beta_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_beta_old(:,:,:,1:nspec)
  close(IIN)

  ! rho kernel
  if (tomo_par%USE_RHO_SCALING_FWAT) then
    ! uses scaling relation with shear perturbations
    kernel_rho_old(:,:,:,:) = RHO_SCALING * kernel_beta_old(:,:,:,:)
    if (myrank == 0) print *,'  rho kernel uses scaling with shear kernel: scaling value = ',RHO_SCALING
  else
    ! uses rho kernel
    fname = 'rhop_kernel_smooth'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) kernel_rho_old(:,:,:,1:nspec)
    close(IIN)
  endif

  ! statistics
  call min_all_cr(minval(kernel_bulk_old),min_vp)
  call max_all_cr(maxval(kernel_bulk_old),max_vp)

  call min_all_cr(minval(kernel_beta_old),min_vs)
  call max_all_cr(maxval(kernel_beta_old),max_vs)

  call min_all_cr(minval(kernel_rho_old),min_rho)
  call max_all_cr(maxval(kernel_rho_old),max_rho)

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

  ! statistics output
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_kernels_minmax',status='unknown')
    write(IOUT,*) '#min_vp #max_vp #min_vs #max_vs  #min_rho #max_rho'
    write(IOUT,'(4e24.12)') min_vp, max_vp, min_vs, max_vs, min_rho, max_rho
    close(IOUT)
  endif

  ! reads in old gradient directions (phi_(n-1))
  USE_OLD_GRADIENT = .true.

  ! checks if files are available:
  fname = 'dbulk'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if (.not. exist) then
 
    print *,'old kernel updates do not exist: ',trim(m_file)
    USE_OLD_GRADIENT = .false.
  endif

  ! makes sure all processes have same flag
  call any_all_l(exist,exist_all)
  if (.not. exist_all) then
    if (myrank == 0) print *,'old kernel updates do not exist for all: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! makes sure all processes have same flag
  use_old_gradient_all = .false.
  call synchronize_all()

  call any_all_l(USE_OLD_GRADIENT,use_old_gradient_all)
  if (.not. use_old_gradient_all) then
    print *,'old kernel updates exists, not consistent for all: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! reads in old gradient
  if (USE_OLD_GRADIENT) then

    ! user output
    if (myrank == 0) print *,'reading old gradient...'

    ! allocate arrays for storing old gradient
    allocate(model_dbulk_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
             model_dbeta_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
             model_drho_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
    if (ier /= 0) stop 'error allocating gradient arrays'

    ! initializes arrays
    model_dbulk_old = 0.0_CUSTOM_REAL
    model_dbeta_old = 0.0_CUSTOM_REAL
    model_drho_old = 0.0_CUSTOM_REAL

    ! dbulk kernel
    fname = 'dbulk'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbulk_old(:,:,:,1:nspec)
    close(IIN)

    ! beta kernel
    fname = 'dbeta'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbeta_old(:,:,:,1:nspec)
    close(IIN)

    ! rho kernel
    fname = 'drho'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_drho_old(:,:,:,1:nspec)
    close(IIN)

    ! statistics
    call min_all_cr(minval(model_dbulk_old),min_vp)
    call max_all_cr(maxval(model_dbulk_old),max_vp)

    call min_all_cr(minval(model_dbeta_old),min_vs)
    call max_all_cr(maxval(model_dbeta_old),max_vs)

    call min_all_cr(minval(model_drho_old),min_rho)
    call max_all_cr(maxval(model_drho_old),max_rho)

    if (myrank == 0) then
      print *
      print *,'old kernel updates:'
      print *,'  bulk min/max : ',min_vp,max_vp
      print *,'  beta min/max: ',min_vs,max_vs
      print *,'  rho min/max  : ',min_rho,max_rho
      print *
    endif

    ! statistics output
    if (PRINT_STATISTICS_FILES .and. myrank == 0) then
      open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_kernel_updates_minmax',status='unknown')
      write(IOUT,*) '#min_vp #max_vp #min_vs #max_vs #min_rho #max_rho'
      write(IOUT,'(4e24.12)') min_vp, max_vp, min_vs, max_vs, min_rho, max_rho
      close(IOUT)
    endif

  endif ! USE_OLD_GRADIENT
  call synchronize_all()

end subroutine read_kernels_cg_iso_old


subroutine read_kernels_cg_tiso_old()

! reads in smoothed kernels from former iteration in OUTPUT_SUM.old/ : bulk, betav, betah, eta

  use tomography_kernels_tiso_cg

  implicit none
  real(kind=CUSTOM_REAL) :: min_vsv,min_vsh,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk
  real(kind=CUSTOM_REAL) :: min_vpv,min_vph,max_vpv,max_vph
  logical:: exist,exist_all,use_old_gradient_all
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'reading old kernels...'

  ! allocate arrays for storing kernels and perturbations
  ! transversely isotropic arrays
  if (USE_ALPHA_BETA_RHO_TISO) then 
     allocate(kernel_alphav_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_alphah_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_betav_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_betah_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_eta_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  else

      allocate(kernel_bulk_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_betav_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_betah_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
           kernel_eta_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  endif

  if (ier /= 0) stop 'error allocating kernel arrays'

  ! initializes arrays
  if (USE_ALPHA_BETA_RHO_TISO) then
  kernel_alphav_old = 0.0_CUSTOM_REAL
  kernel_alphah_old = 0.0_CUSTOM_REAL
  else
  kernel_bulk_old = 0.0_CUSTOM_REAL
  endif
  kernel_betav_old = 0.0_CUSTOM_REAL
  kernel_betah_old = 0.0_CUSTOM_REAL
  kernel_eta_old = 0.0_CUSTOM_REAL

  ! checks if files are available:
  if (USE_ALPHA_BETA_RHO_TISO) then
  ! alphav kernel
  fname = 'alphav_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_alphav_old(:,:,:,1:nspec)
  close(IIN)

  ! alphah kernel
  fname = 'alphah_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_alphah_old(:,:,:,1:nspec)
  close(IIN)


  ! betav kernel
  fname = 'betav_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betav_old(:,:,:,1:nspec)
  close(IIN)

  ! betah kernel
  fname = 'betah_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betah_old(:,:,:,1:nspec)
  close(IIN)

  else

  fname = 'bulk_c_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if (.not. exist) then
    print *,'Error file does not exist: ',trim(m_file)
    call exit_mpi(myrank,'file not exist')
  endif
  ! makes sure all processes have same flag
  call any_all_l(exist,exist_all)
  if (.not. exist_all) then
    print *,'old kernels do not exist: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! bulk kernel
  fname = 'bulk_c_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_bulk_old(:,:,:,1:nspec)
  close(IIN)

  ! betav kernel
  fname = 'bulk_betav_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betav_old(:,:,:,1:nspec)
  close(IIN)

  ! betah kernel
  fname = 'bulk_betah_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betah_old(:,:,:,1:nspec)
  close(IIN)
  endif

  ! eta kernel
  fname = 'eta_kernel_smooth'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_eta_old(:,:,:,1:nspec)
  close(IIN)

  ! statistics
  if (USE_ALPHA_BETA_RHO_TISO) then
  call min_all_cr(minval(kernel_alphah_old),min_vph)
  call max_all_cr(maxval(kernel_alphah_old),max_vph)

  call min_all_cr(minval(kernel_alphav_old),min_vpv)
  call max_all_cr(maxval(kernel_alphav_old),max_vpv)
  else 

  call min_all_cr(minval(kernel_bulk_old),min_bulk)
  call max_all_cr(maxval(kernel_bulk_old),max_bulk)
  endif

  call min_all_cr(minval(kernel_betah_old),min_vsh)
  call max_all_cr(maxval(kernel_betah_old),max_vsh)

  call min_all_cr(minval(kernel_betav_old),min_vsv)
  call max_all_cr(maxval(kernel_betav_old),max_vsv)

  call min_all_cr(minval(kernel_eta_old),min_eta)
  call max_all_cr(maxval(kernel_eta_old),max_eta)

  if (myrank == 0) then
    print *
    print *,'old kernels:'
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

  ! reads in old gradient directions (phi_(n-1))
  USE_OLD_GRADIENT = .true.

  ! checks if files are available:

  if (USE_ALPHA_BETA_RHO_TISO) then
  fname = 'dalphav'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if (.not. exist) then
    print *,'old kernel updates do not exist: ',trim(m_file)
    USE_OLD_GRADIENT = .false.
  endif
  else 
  fname = 'dbulk_c'
  write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if (.not. exist) then
    print *,'old kernel updates do not exist: ',trim(m_file)
    USE_OLD_GRADIENT = .false.
  endif
  endif
  ! makes sure all processes have same flag
  call any_all_l(exist,exist_all)
  if (.not. exist_all) then
    if (myrank == 0) print *,'old kernel updates do not exist for all: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! makes sure all processes have same flag
  use_old_gradient_all = .false.
  call synchronize_all()

  call any_all_l(USE_OLD_GRADIENT,use_old_gradient_all)
  if (.not. use_old_gradient_all) then
    print *,'old kernel updates exists, not consistent for all: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! reads in old gradient
  if (USE_OLD_GRADIENT) then

    ! user output
    if (myrank == 0) print *,'reading old gradient...'

    ! allocate arrays for storing old gradient
    ! transversely isotropic arrays
    if (USE_ALPHA_BETA_RHO_TISO) then
    allocate(model_dalphav_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
             model_dalphah_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
             model_dbetav_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
             model_dbetah_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
             model_deta_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
    else 
    allocate(model_dbulk_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
             model_dbetav_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
             model_dbetah_old(NGLLX,NGLLY,NGLLZ,NSPEC), &
             model_deta_old(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
    endif
    if (ier /= 0) stop 'error allocating gradient arrays'

    ! initializes arrays
    if (USE_ALPHA_BETA_RHO_TISO) then
    model_dalphav_old = 0.0_CUSTOM_REAL
    model_dalphah_old = 0.0_CUSTOM_REAL
    else 
    model_dbulk_old = 0.0_CUSTOM_REAL
    endif
    model_dbetav_old = 0.0_CUSTOM_REAL
    model_dbetah_old = 0.0_CUSTOM_REAL
    model_deta_old = 0.0_CUSTOM_REAL

    if (USE_ALPHA_BETA_RHO_TISO) then
    ! alphav kernel
    fname = 'dalphav'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dalphav_old(:,:,:,1:nspec)
    close(IIN)

    ! alphah kernel
    fname = 'dalphah'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dalphah_old(:,:,:,1:nspec)
    close(IIN)
    else

    ! bulk kernel
    fname = 'dbulk_c'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbulk_old(:,:,:,1:nspec)
    close(IIN)
    endif

    ! betav kernel
    fname = 'dbetav'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbetav_old(:,:,:,1:nspec)
    close(IIN)

    ! betah kernel
    fname = 'dbetah'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbetah_old(:,:,:,1:nspec)
    close(IIN)

    ! eta kernel
    fname = 'deta'
    write(m_file,'(a,i6.6,a)') trim(KERNEL_OLD_DIR)//'/proc',myrank,trim(REG)//trim(fname)//'.bin'
    if (myrank == 0) print *,'  '//trim(KERNEL_OLD_DIR)//'/proc**'//trim(REG)//trim(fname)//'.bin'

    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_deta_old(:,:,:,1:nspec)
    close(IIN)

    ! statistics
    if (USE_ALPHA_BETA_RHO_TISO) then
    call min_all_cr(minval(model_dalphah_old),min_vph)
    call max_all_cr(maxval(model_dalphah_old),max_vph)

    call min_all_cr(minval(model_dalphav_old),min_vpv)
    call max_all_cr(maxval(model_dalphav_old),max_vpv)
    else

    call min_all_cr(minval(model_dbulk_old),min_bulk)
    call max_all_cr(maxval(model_dbulk_old),max_bulk)
    endif

    call min_all_cr(minval(model_dbetah_old),min_vsh)
    call max_all_cr(maxval(model_dbetah_old),max_vsh)

    call min_all_cr(minval(model_dbetav_old),min_vsv)
    call max_all_cr(maxval(model_dbetav_old),max_vsv)

    call min_all_cr(minval(model_deta_old),min_eta)
    call max_all_cr(maxval(model_deta_old),max_eta)

    if (myrank == 0) then
      print *
      print *,'old kernel updates:'
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

    ! statistics output
    if (PRINT_STATISTICS_FILES .and. myrank == 0) then
      open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_kernel_updates_minmax',status='unknown')
      if (USE_ALPHA_BETA_RHO_TISO) then
      write(IOUT,*) '#min_vpv #max_vpv #min_vph #max_vph #min_vsv #max_vsv #min_vsh #max_vsh #min_eta #max_eta'
      write(IOUT,'(4e24.12)') min_vpv, max_vpv, min_vph, max_vph, min_vsv, max_vsv, min_vsh, max_vsh, min_eta, max_eta
      else 
      write(IOUT,*) '#min_bulk #max_bulk #min_vsv #max_vsv #min_vsh #max_vsh #min_eta #max_eta'
      write(IOUT,'(4e24.12)') min_bulk, max_bulk, min_vsv, max_vsv, min_vsh, max_vsh, min_eta, max_eta
      endif
      close(IOUT)
    endif

  endif ! USE_OLD_GRADIENT
  call synchronize_all()

end subroutine read_kernels_cg_tiso_old

