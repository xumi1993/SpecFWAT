module post_processing
  use config
  use fwat_mpi
  use fwat_constants
  use logger, only: log
  use input_params, fpar => fwat_par_global
  use utils, only: zeros
  use kernel_io
  use taper3d
  use common_lib, only: get_dat_type, get_kernel_names, mkdir
  use smooth_mod
  use zprecond
  use shared_parameters
  use model_grid_data, only: read_grid_kernel_smooth, create_grid,&
                              write_grid_kernel_smooth, write_grid, gll2grid

  implicit none
  character(len=MAX_STRING_LEN), private :: msg

  type PostFlow
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: ker_names
    real(kind=cr), dimension(:,:,:,:,:), allocatable :: ker_data
    real(kind=cr), dimension(:,:,:,:), pointer :: ker_data_smooth
    real(kind=cr), dimension(:,:,:), allocatable :: hess_smooth
    integer :: nker, ks_win
    character(len=MAX_STRING_LEN) :: kernel_path
    contains
    procedure :: sum_kernel, init=>init_post_flow, sum_precond, init_for_type,apply_precond, &
                 write_gradient_grid, remove_ekernel, pde_smooth, taper_kernel_grid, finalize, &
                 read_sum_kernel
  end type PostFlow
contains

  subroutine init_post_flow(this, is_read_database)
    class(PostFlow), intent(inout) :: this
    logical, intent(in) :: is_read_database
    integer :: itype

    call get_kernel_names()

    call create_grid()
    local_path_backup = LOCAL_PATH

    if (fpar%update%MODEL_TYPE == 1) then
      ANISOTROPIC_KL = .false.
      ANISOTROPY = .false.
    else
      ANISOTROPIC_KL = .true.
      ANISOTROPY = .true.
    endif

    call log%init('output_post_processing_'//trim(model_name)//'.log')
    call log%write('*******************************************', .false.)

  end subroutine init_post_flow

  subroutine init_for_type(this, itype)
    class(PostFlow), intent(inout) :: this
    integer, intent(in) :: itype

    ! set simu type
    simu_type = INV_TYPE_NAMES(itype)
    LOCAL_PATH = local_path_backup
    call fpar%select_simu_type()
    LOCAL_PATH = local_path_fwat

    ! read src_rec for this data type
    call get_dat_type()

    ! read src_rec for this data type
    call fpar%acqui%read()
    
    ! setup mesh
    call read_mesh_databases_for_init()
    
    call log%write('Simulation type: '//trim(simu_type), .false.)
    write(msg, '(a,F10.1,F10.1)') 'PDE smoothing; SIGMA_H, SIGMA_V: ', fpar%sim%SIGMA_H, &
                                fpar%sim%SIGMA_V
    call log%write(msg, .false.)
    if (.not. fpar%postproc%IS_PRECOND .and. is_joint) then
      write(msg, '(a)') 'Preconditioned L-BFGS'
    else
      write(msg, '(a, L5)') 'Preconditioning: ', fpar%postproc%IS_PRECOND
    endif
    call log%write(msg, .false.)
    if (.not. fpar%postproc%IS_PRECOND .and. is_joint) then
      write(msg, '(a,I3)') 'Precondition type: ', tele_par%PRECOND_TYPE
    else
      write(msg, '(a,I3)') 'Precondition type: ', fpar%sim%PRECOND_TYPE
    endif
    call log%write(msg, .false.)
    write(msg, '(a, L5)') 'USE_RHO_SCALING: ', fpar%sim%USE_RHO_SCALING
    call log%write(msg, .false.)
    call log%write('-------------------------------------------', .false.)

    if (is_joint) then
      this%kernel_path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)//'_'//trim(simu_type)
    else
      this%kernel_path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)
    endif
    call mkdir(this%kernel_path)
    call synchronize_all()

    this%ker_data = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_FWAT, nkernel)
    call prepare_shm_array_cr_4d(this%ker_data_smooth, ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel, this%ks_win)
  end subroutine init_for_type
  
  subroutine sum_kernel(this)
    class(PostFlow), intent(inout) :: this
    real(kind=cr), dimension(:,:,:,:), allocatable :: ker
    real(kind=cr), dimension(:,:,:,:,:), allocatable :: ker_aniso
    integer :: iker, ievt
    real(kind=cr) :: max_loc, min_loc, max_glob, min_glob

    call log%write('This is taking sum of kernels...', .true.)
    if (ANISOTROPIC_KL) then
      ! take sum of anisotropic kernels
      if (fpar%update%MODEL_TYPE == 2) then
        do ievt = 1, fpar%acqui%nevents
          call kernel_cijkl2hti(ievt, ker_aniso)
          this%ker_data = this%ker_data + ker_aniso
        enddo
      endif
    else
      do iker = 1, nkernel
        if(fpar%sim%USE_RHO_SCALING .and. (kernel_names(iker) == 'rhop')) cycle
        do ievt = 1, fpar%acqui%nevents
          call read_event_kernel(ievt, trim(kernel_names(iker))//'_kernel', ker)
          this%ker_data(:,:,:,:,iker) = this%ker_data(:,:,:,:,iker) + ker
        enddo
      enddo
    endif
    if (is_output_sum_kernel .or. fpar%update%DO_LS) then
      call log%write('This is writing sum of kernels...', .true.)
      do iker = 1, nkernel
        call write_kernel(this%kernel_path, trim(kernel_names(iker))//'_kernel', this%ker_data(:,:,:,:,iker))
      enddo
    endif
    max_loc = maxval(this%ker_data)
    min_loc = minval(this%ker_data)
    call max_all_all_cr(max_loc, max_glob)
    call min_all_all_cr(min_loc, min_glob)
    write(msg, '(a,G0.9)') 'Max kernel value: ', max_glob
    call log%write(msg, .false.)
    write(msg, '(a,G0.9)') 'Min kernel value: ', min_glob
    call log%write(msg, .false.)
    call synchronize_all()
  
  end subroutine sum_kernel

  subroutine read_sum_kernel(this)
    class(PostFlow), intent(inout) :: this
    integer :: iker
    character(len=MAX_STRING_LEN) :: fname
    real(kind=cr), dimension(:,:,:,:), allocatable :: ker

    call log%write('This is reading sum of kernels...', .true.)
    do iker = 1, nkernel
      call read_kernel(trim(this%kernel_path), trim(kernel_names(iker))//'_kernel', ker)
      this%ker_data(:,:,:,:,iker) = ker
    enddo
    call synchronize_all()
  end subroutine

  subroutine apply_precond(this)
    class(PostFlow), intent(inout) :: this
    integer :: iker,  ievt, i, j, k, iglob, ispec
    real(kind=cr), dimension(:,:,:,:), allocatable :: total_hess, ker

    call log%write('This is applying preconditioned kernels...', .true.)

    if (fpar%sim%PRECOND_TYPE <= DEFAULT_PRECOND) then
      if (run_mode == 1) then
        total_hess = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_FWAT)
        do ievt = 1, fpar%acqui%nevents
          call read_event_kernel(ievt, 'hess_kernel', ker)
          if (fpar%sim%PRECOND_TYPE == DEFAULT_PRECOND) then
            total_hess = total_hess + abs(ker)
          else
            total_hess = total_hess + ker
          endif
        enddo
        if (fpar%sim%PRECOND_TYPE == 0) total_hess = abs(total_hess)
        if (is_output_sum_kernel) call write_kernel(this%kernel_path, 'hess_kernel', total_hess)
        call invert_hess(total_hess)
        if (is_output_sum_kernel) call write_kernel(this%kernel_path, 'inv_hess_kernel', total_hess)
      elseif (run_mode == 2) then
        call read_kernel(this%kernel_path, 'inv_hess_kernel', total_hess)
      endif
    else
      call zprecond_gll(total_hess) 
    endif

    do iker = 1, nkernel
      this%ker_data(:,:,:,:,iker) = this%ker_data(:,:,:,:,iker) * total_hess
    enddo
    call synchronize_all()
  end subroutine apply_precond

  subroutine sum_precond(this)
    class(PostFlow), intent(inout) :: this
    character(len=MAX_STRING_LEN), parameter :: hess_name = 'hess'
    real(kind=cr), dimension(:,:,:,:), allocatable :: total_hess, total_hess_smooth, ker
    real(kind=cr), dimension(:), allocatable :: z_precond
    real(kind=cr), dimension(:,:), allocatable :: gk
    real(kind=cr), dimension(:,:,:), allocatable :: gm
    real(kind=cr) :: precond
    character(len=MAX_STRING_LEN) :: fname
    integer :: iker, ievt, i, j, k, iglob

    call log%write('This is saving preconditioned kernels...', .true.)
    if (fpar%sim%PRECOND_TYPE <= 1) then
      total_hess = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_FWAT)
      do ievt = 1, fpar%acqui%nevents
        call read_event_kernel(ievt, trim(hess_name)//'_kernel', ker)
        if (fpar%sim%PRECOND_TYPE == DEFAULT_PRECOND) then
          total_hess = total_hess + abs(ker)
        else
          total_hess = total_hess + ker
        endif
      enddo
      if (fpar%sim%PRECOND_TYPE == 0) total_hess = abs(total_hess)
      if (is_output_sum_kernel) call write_kernel(this%kernel_path, trim(hess_name)//'_kernel', total_hess)
      call invert_hess(total_hess)
      if (is_output_sum_kernel) call write_kernel(this%kernel_path, trim(hess_name)//'_inv_kernel', total_hess)
      call smooth_sem_pde(total_hess, fpar%sim%sigma_h, fpar%sim%sigma_v, total_hess_smooth, .false.)
      call gll2grid(total_hess_smooth, gm)
      this%hess_smooth = gm/maxval(abs(gm))
    else
      call zprecond_grid(this%hess_smooth)
    endif
    call synchronize_all()
    if (worldrank == 0) then
      fname = trim(OPT_DIR)//'/'//trim(HESS_PREFIX)//'_'//trim(model_name)//'.h5'
      call write_grid(fname, HESS_PREFIX, this%hess_smooth)
    endif
    call synchronize_all()
    
  end subroutine sum_precond

  function set_z_precond(z) result(precond)
    use utils, only: arange
    real(kind=cr) :: precond
    real(kind=cr) :: z_max, min_z_pre=1.0_cr
    real(kind=cr), intent(in) :: z
    integer :: ndep

    if (z >= 0) then
      precond = min_z_pre
    else
      precond = z-min_z_pre
    endif
    if (tele_par%PRECOND_TYPE == Z_PRECOND) then
      precond = sqrt(precond**2)
      z_max = abs(z_min_glob)
    elseif (tele_par%PRECOND_TYPE == Z_SQRT_PRECOND) then
      precond = sqrt(sqrt(precond**2))
      z_max = sqrt(abs(z_min_glob))
    endif
    precond = precond / z_max 
  end function set_z_precond

  subroutine pde_smooth(this)
    class(PostFlow), intent(inout) :: this
    integer :: iker
    real(kind=cr), dimension(:,:,:,:), allocatable :: gk
    real(kind=cr), dimension(:,:,:), allocatable :: gm

    do iker = 1, nkernel
      if((.not. ANISOTROPIC_KL) .and. fpar%sim%USE_RHO_SCALING .and. (kernel_names(iker) == 'rhop')) then
        call log%write('This is scaling for rhop kernels...', .true.)
        if (worldrank == 0) this%ker_data_smooth(:,:,:,iker) = this%ker_data_smooth(:,:,:,2) * RHO_SCALING_FAC
      else
        call log%write('This is smoothing of '//trim(kernel_names(iker))//' kernels...', .true.)
        call smooth_sem_pde(this%ker_data(:,:,:,:,iker), fpar%sim%sigma_h, fpar%sim%sigma_v, gk, .false.)
        call gll2grid(gk, gm)
        if (worldrank == 0) this%ker_data_smooth(:,:,:,iker) = gm
      endif
    end do
    call synchronize_all()
    call sync_from_main_rank_cr_4d(this%ker_data_smooth, ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel)
  end subroutine pde_smooth

  subroutine write_gradient_grid(this)
    class(PostFlow), intent(inout) :: this
    integer :: iker, ievt
    character(len=MAX_STRING_LEN) :: fname, suffix
    logical :: is_simu_type
    
    if (is_joint) then
      suffix = trim(model_name)//'_'//trim(simu_type)
    else
      suffix = trim(model_name)
    endif
    fname = trim(OPT_DIR)//'/gradient_'//trim(suffix)//'.h5'

    call write_grid_kernel_smooth(this%ker_data_smooth, fname)

  end subroutine write_gradient_grid

  subroutine taper_kernel_grid(this)
    class(PostFlow), intent(inout) :: this
    type(taper_cls) :: tap
    integer :: iker, i, j, k
    real(kind=cr) :: dh, val
    real(kind=cr), dimension(:,:,:), allocatable :: taper
    call log%write('This is tapering kernels...', .true.)

    ! create regluar densey taper grid
    if (worldrank == 0) then
      dh = (distance_max_glob+distance_min_glob)/2.0_cr
      call tap%create(x_min_glob, x_max_glob,&
                    y_min_glob, y_max_glob,&
                    z_min_glob, z_max_glob, &
                    dh, dh, dh, &
                    fpar%postproc%TAPER_H_SUPPRESS, &
                    fpar%postproc%TAPER_H_BUFFER, &
                    fpar%postproc%TAPER_V_SUPPRESS, &
                    fpar%postproc%TAPER_V_BUFFER)

      ! create mask grid
      taper = zeros(ext_grid%nx, ext_grid%ny, ext_grid%nz)
      do i = 1, ext_grid%nx
        do j = 1, ext_grid%ny
          do k = 1, ext_grid%nz
            taper(i,j,k) = tap%Interp(ext_grid%x(i), ext_grid%y(j), ext_grid%z(k))
          enddo
        enddo
      enddo

      ! apply tape for each kernel
      do iker = 1, nkernel
        this%ker_data_smooth(:,:,:,iker) = this%ker_data_smooth(:,:,:,iker) * taper
      enddo
    endif
    call sync_from_main_rank_cr_4d(this%ker_data_smooth, ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel)
    call synchronize_all()

  end subroutine taper_kernel_grid

  subroutine invert_hess( hess_matrix )
  ! inverts the Hessian matrix
  ! the approximate Hessian is only defined for diagonal elements: like
  ! H_nn = \frac{ \partial^2 \chi }{ \partial \rho_n \partial \rho_n }
  ! on all GLL points, which are indexed (i,j,k,ispec)

    real(kind=cr), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_FWAT) :: hess_matrix

    ! local parameters
    real(kind=cr) :: maxh,maxh_all

    ! maximum value of Hessian
    maxh = maxval( abs(hess_matrix) )

    ! determines maximum from all slices on main
    call max_all_all_cr(maxh, maxh_all)

    ! normalizes Hessian
    if (maxh_all < 1.e-18) then
      ! Hessian is zero, re-initializes
      hess_matrix = 1.0_cr
      !stop 'Error Hessian too small'
    else
      ! since Hessian has absolute values, this scales between [0,1]
      hess_matrix = hess_matrix / maxh_all
    endif


    ! inverts Hessian values
    where( abs(hess_matrix(:,:,:,:)) > THRESHOLD_HESS )
      hess_matrix = 1.0_cr / hess_matrix
    elsewhere
      hess_matrix = 1.0_cr / THRESHOLD_HESS
    endwhere

    maxh = maxval( abs(hess_matrix) )
    call max_all_all_cr(maxh, maxh_all)
    ! normalizes Hessian
    hess_matrix = hess_matrix / maxh_all

    ! rescales Hessian
    !hess_matrix = hess_matrix * maxh_all

  end subroutine invert_hess

  subroutine remove_ekernel(this)
    class(PostFlow), intent(inout) :: this
    integer :: ievt

    if (is_output_event_kernel) return

    do ievt = 1, fpar%acqui%nevents
      if (ANISOTROPIC_KL) then
        if (SAVE_TRANSVERSE_KL) then
          call remove_event_kernel(ievt, 'alphav_kernel')
          call remove_event_kernel(ievt, 'alphah_kernel')
          call remove_event_kernel(ievt, 'betav_kernel')
          call remove_event_kernel(ievt, 'betah_kernel')
          call remove_event_kernel(ievt, 'eta_kernel')
        else
          call remove_event_kernel(ievt, 'rho_kernel')
          call remove_event_kernel(ievt, 'c11_kernel')
          call remove_event_kernel(ievt, 'c12_kernel')
          call remove_event_kernel(ievt, 'c13_kernel')
          call remove_event_kernel(ievt, 'c14_kernel')
          call remove_event_kernel(ievt, 'c15_kernel')
          call remove_event_kernel(ievt, 'c16_kernel')
          call remove_event_kernel(ievt, 'c22_kernel')
          call remove_event_kernel(ievt, 'c23_kernel')
          call remove_event_kernel(ievt, 'c24_kernel')
          call remove_event_kernel(ievt, 'c25_kernel')
          call remove_event_kernel(ievt, 'c26_kernel')
          call remove_event_kernel(ievt, 'c33_kernel')
          call remove_event_kernel(ievt, 'c34_kernel')
          call remove_event_kernel(ievt, 'c35_kernel')
          call remove_event_kernel(ievt, 'c36_kernel')
          call remove_event_kernel(ievt, 'c44_kernel')
          call remove_event_kernel(ievt, 'c45_kernel')
          call remove_event_kernel(ievt, 'c46_kernel')
          call remove_event_kernel(ievt, 'c55_kernel')
          call remove_event_kernel(ievt, 'c56_kernel')
          call remove_event_kernel(ievt, 'c66_kernel')
        endif
      else
        call remove_event_kernel(ievt, 'rho_kernel')
        call remove_event_kernel(ievt, 'mu_kernel')
        call remove_event_kernel(ievt, 'kappa_kernel')
        call remove_event_kernel(ievt, 'rhop_kernel')
        call remove_event_kernel(ievt, 'beta_kernel')
        call remove_event_kernel(ievt, 'alpha_kernel')
      endif
      call remove_event_kernel(ievt, 'hess_kernel')
    enddo
    call synchronize_all()

  end subroutine remove_ekernel

  subroutine sum_joint_kernel_grid()
    integer :: itype, iker
    character(len=MAX_STRING_LEN) :: type_name, fname
    real(kind=cr), dimension(:,:,:,:), allocatable :: kernel, total_kernel
    real(kind=cr) :: norm_val(NUM_INV_TYPE), std_val, linf_val, mean_val

    call log%write('This is taking sum of kernels for joint inversion...', .true.)
    norm_val = 0._cr
    ! calculate normalization value
    ! set normalization value 
    do itype = 1, NUM_INV_TYPE
      if (.not. fpar%postproc%INV_TYPE(itype)) cycle
      if (fpar%postproc%NORM_TYPE == 1 .or. fpar%postproc%NORM_TYPE == 3) then
        if (worldrank == 0) then
          call calc_kernel0_weight_grid(itype, linf_val, std_val)
          if (fpar%postproc%NORM_TYPE == 1) then
            norm_val(itype) = linf_val
          else
            norm_val(itype) = std_val
          endif
        endif
      elseif (fpar%postproc%NORM_TYPE == 2) then
        call calc_misfit0_weight(itype, norm_val(itype))
      else
        call log%write('ERROR: Unknown normalization type for joint inversion', .false.)
        call exit_MPI(0, 'Unknown normalization type for joint inversion')
      endif
    enddo
    call bcast_all_cr(norm_val, NUM_INV_TYPE)

    if (worldrank == 0) then
      total_kernel = zeros(ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel)
      do itype = 1, NUM_INV_TYPE
        if (.not. fpar%postproc%INV_TYPE(itype)) cycle
        type_name = INV_TYPE_NAMES(itype)

        ! read gradient for itype
        fname = trim(OPT_DIR)//'/gradient_'//trim(model_name)//'_'//trim(type_name)//'.h5'
        call read_grid_kernel_smooth(fname, kernel)

        ! take sum of kernels
        kernel = fpar%postproc%JOINT_WEIGHT(itype)*kernel/norm_val(itype)
        total_kernel = total_kernel + kernel

        ! calculate infinity norm and std for logger
        mean_val = sum(kernel) / (ext_grid%nx * ext_grid%ny * ext_grid%nz * nkernel)
        std_val = sqrt(sum((kernel - mean_val)**2) / (ext_grid%nx * ext_grid%ny * ext_grid%nz * nkernel))
        linf_val = maxval(abs(kernel))
        write(msg, '(a,G18.9,"  ",G18.9)') 'infinity-norm, std for '//trim(type_name)//': ', linf_val, std_val
        print *, linf_val, std_val
        call log%write(msg, .false.)
      enddo
      fname = trim(OPT_DIR)//'/gradient_'//trim(model_name)//'.h5'
      call write_grid_kernel_smooth(total_kernel, fname)
    endif

  end subroutine sum_joint_kernel_grid

  subroutine calc_kernel0_weight_grid(itype, linf_val, std_val)
    integer, intent(in) :: itype
    real(kind=cr), dimension(:,:,:,:), allocatable :: kernel_data
    real(kind=cr) :: mean_val
    real(kind=cr), intent(out) :: linf_val, std_val
    character(len=MAX_STRING_LEN) :: fname, modname

    ! read gradient of iter_start for itype 
    if (fpar%update%OPT_METHOD == 2) then
      write(modname, '("M",I2.2)') fpar%update%ITER_START
    else
      modname = 'M00'
    endif
    fname = trim(OPT_DIR)//'/gradient_'//trim(modname)//'_'//trim(INV_TYPE_NAMES(itype))//'.h5'
    call read_grid_kernel_smooth(fname, kernel_data)

    ! calculate infinity norm and std for normalization
    linf_val = maxval(abs(kernel_data))
    mean_val = sum(kernel_data) / (ext_grid%nx * ext_grid%ny * ext_grid%nz * nkernel)
    std_val = sqrt(sum((kernel_data - mean_val)**2) / (ext_grid%nx * ext_grid%ny * ext_grid%nz * nkernel))

  end subroutine calc_kernel0_weight_grid

  subroutine calc_misfit0_weight(itype, weight)
    use window_chi, only: read_model_misfit
    use common_lib, only: get_dat_type
    integer, intent(in) :: itype
    real(kind=cr), intent(out) :: weight
    real(kind=dp) :: misfit_loc, misfit_all
    character(len=MAX_STRING_LEN) :: modname
    integer :: ievt

    ! read total misfit of iter_start for itype
    if (fpar%update%OPT_METHOD == 2) then
      write(modname, '("M",I2.2)') fpar%update%ITER_START
    else
      modname = 'M00'
    endif

    ! set simu type
    simu_type = INV_TYPE_NAMES(itype)
    LOCAL_PATH = local_path_backup
    call fpar%select_simu_type()
    LOCAL_PATH = local_path_fwat

    call get_dat_type()

    call fpar%acqui%read()

    ! calculate total misfit
    misfit_all = 0._dp
    do ievt = 1, fpar%acqui%nevents
      call read_model_misfit(modname, ievt, misfit_loc)
      misfit_all = misfit_all + misfit_loc
    enddo
    weight = real(misfit_all)

    call fpar%acqui%finalize()

  end subroutine calc_misfit0_weight

  subroutine finalize(this)
    class(PostFlow), intent(inout) :: this
    call log%write('This is finalizing post-processing of '//trim(simu_type)//'...', .true.)
    call log%write('*******************************************', .false.)
    call free_shm_array(this%ks_win)
  end subroutine finalize

end module post_processing