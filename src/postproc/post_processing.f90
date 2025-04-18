module post_processing
  use config
  use fwat_mpi
  use fwat_constants
  use logger, only: log
  use input_params, fpar => fwat_par_global
  use utils, only: zeros
  use kernel_io
  use taper3d
  use common_lib, only: get_dat_type, get_kernel_names
  use multigrid
  use shared_parameters
  use model_grid_data, only: read_grid_kernel_smooth, create_grid,&
                              write_grid_kernel_smooth, write_grid

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
    procedure :: sum_kernel, init=>init_post_flow, sum_precond, init_for_type,&
                 write_gradient_grid, remove_ekernel, multigrid_smooth, taper_kernel_grid, finalize
  end type PostFlow
contains

  subroutine init_post_flow(this, is_read_database)
    class(PostFlow), intent(inout) :: this
    logical, intent(in) :: is_read_database
    integer :: itype

    ! call read_mesh_databases_minimum(is_read_database)

    call get_kernel_names()

    call create_grid()

    ! call system('mkdir -p '//trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name))

    call log%init('output_post_processing_'//trim(model_name)//'.log')
    call log%write('*******************************************', .false.)

  end subroutine init_post_flow

  subroutine init_for_type(this, itype)
    class(PostFlow), intent(inout) :: this
    integer, intent(in) :: itype

    ! set simu type
    simu_type = INV_TYPE_NAMES(itype)
    call fpar%select_simu_type()

    ! read src_rec for this data type
    call get_dat_type()

    ! read src_rec for this data type
    call fpar%acqui%read()
    
    ! setup mesh
    call get_mesh_coord()
    
    call log%write('Simulation type: '//trim(simu_type), .false.)
    write(msg, '(a,I4,I4,I4)') 'Multi-grid smoothing: ', fpar%postproc%ninv(1), &
                                fpar%postproc%ninv(2), fpar%postproc%ninv(3)
    call log%write(msg, .false.)
    if (is_joint) then
      write(msg, '(a)') 'Preconditioned L-BFGS'
    else
      write(msg, '(a, L5)') 'Preconditioning: ', fpar%postproc%IS_PRECOND
    endif
    call log%write(msg, .false.)
    if (is_joint) then
      write(msg, '(a,I3)') 'Precondition type: ', tele_par%PRECOND_TYPE
    else
      write(msg, '(a,I3)') 'Precondition type: ', fpar%sim%PRECOND_TYPE
    endif
    call log%write(msg, .false.)
    write(msg, '(a, L5)') 'USE_RHO_SCALING: ', fpar%sim%USE_RHO_SCALING
    call log%write(msg, .false.)
    call log%write('-------------------------------------------', .false.)

    if (worldrank == 0) then
      if (is_joint) then
        this%kernel_path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)//'_'//trim(simu_type)
      else
        this%kernel_path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)
      endif
      call system('mkdir -p '//trim(this%kernel_path))
    endif
    call bcast_all_ch_array(this%kernel_path, 1, MAX_STRING_LEN)
    call synchronize_all()

    this%ker_data = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_FWAT, nkernel)
    call prepare_shm_array_cr_4d(this%ker_data_smooth, MEXT_V%nx, MEXT_V%ny, MEXT_V%nz, nkernel, this%ks_win)
  end subroutine init_for_type
  
  subroutine sum_kernel(this)
    class(PostFlow), intent(inout) :: this
    real(kind=cr), dimension(:,:,:,:), allocatable :: ker
    integer :: iker, ievt
    real(kind=cr) :: max_loc, min_loc, max_glob, min_glob

    call log%write('This is taking sum of kernels...', .true.)
    do iker = 1, nkernel
      if((.not. ANISOTROPIC_KL) .and. fpar%sim%USE_RHO_SCALING .and. (kernel_names(iker) == 'rhop')) cycle
      do ievt = 1, fpar%acqui%nevents
        call read_event_kernel(ievt, trim(kernel_names(iker))//'_kernel', ker)
        this%ker_data(:,:,:,:,iker) = this%ker_data(:,:,:,:,iker) + ker
      enddo
      if (is_output_event_kernel) then
        call write_kernel(this%kernel_path, trim(kernel_names(iker))//'_kernel', this%ker_data(:,:,:,:,iker))
      endif
    enddo
    max_loc = maxval(this%ker_data)
    min_loc = minval(this%ker_data)
    call max_all_all_cr(max_loc, max_glob)
    call min_all_all_cr(min_loc, min_glob)
    write(msg, '(a,E0.9)') 'Max kernel value: ', max_glob
    call log%write(msg, .false.)
    write(msg, '(a,E0.9)') 'Min kernel value: ', min_glob
    call log%write(msg, .false.)
    call synchronize_all()
  
  end subroutine sum_kernel

  subroutine sum_precond(this)
    class(PostFlow), intent(inout) :: this
    character(len=MAX_STRING_LEN), parameter :: hess_name = 'hess'
    real(kind=cr), dimension(:,:,:,:), allocatable :: total_hess, ker
    real(kind=cr), dimension(:), allocatable :: z_precond
    real(kind=cr), dimension(:,:), allocatable :: gk
    real(kind=cr), dimension(:,:,:), allocatable :: gm
    real(kind=cr) :: precond
    character(len=MAX_STRING_LEN) :: fname
    integer :: iker, ievt, i, j, k, iglob
    type(InvGrid) :: inv

    call log%write('This is applying preconditioned kernels...', .true.)
    if (fpar%sim%PRECOND_TYPE == DEFAULT_PRECOND) then
      total_hess = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_FWAT)
      do ievt = 1, fpar%acqui%nevents
        call read_event_kernel(ievt, trim(hess_name)//'_kernel', ker)
        total_hess = total_hess + abs(ker)
      enddo
      call invert_hess(total_hess)
      ! if (IS_OUTPUT_HESS_INV) call write_kernel(trim(hess_name)//'_inv_kernel', total_hess, .false.)

      ! do iker = 1, nkernel
      !   if((.not. ANISOTROPIC_KL) .and. fpar%sim%USE_RHO_SCALING .and. (kernel_names(iker) == 'rhop')) cycle
      !   this%ker_data(:,:,:,:,iker) = this%ker_data(:,:,:,:,iker) * total_hess
      ! enddo
      call inv%init()
      call inv%sem2inv(total_hess, gk)
      call inv%inv2grid(gk, gm)
      this%hess_smooth = gm
    else
      if (worldrank == 0) then
        this%hess_smooth = zeros(MEXT_V%nx, MEXT_V%ny, MEXT_V%nz)
        do i = 1, MEXT_V%nx
          do j = 1, MEXT_V%ny
            do k = 1, MEXT_V%nz
              this%hess_smooth(i,j,k) = set_z_precond(MEXT_V%z(k))
            enddo
          enddo
        enddo
      endif
    endif
    call synchronize_all()
    if (worldrank == 0) then
      fname = trim(OPT_DIR)//'/'//trim(HESS_PREFIX)//'_'//trim(model_name)//'.h5'
      call write_grid(fname, HESS_PREFIX, this%hess_smooth)
    endif
    
  end subroutine sum_precond

  function set_z_precond(z) result(precond)
    use utils, only: arange
    real(kind=cr) :: precond
    real(kind=cr) :: z_max
    real(kind=cr), intent(in) :: z
    integer :: ndep

    if (z > -THRESHOLD_HESS) then
      precond = THRESHOLD_HESS
    else
      precond = z
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

  subroutine multigrid_smooth(this)
    class(PostFlow), intent(inout) :: this
    type(InvGrid) :: inv
    integer :: iker
    real(kind=cr), dimension(:,:), allocatable :: gk
    real(kind=cr), dimension(:,:,:), allocatable :: gm

    call inv%init()

    do iker = 1, nkernel
      if((.not. ANISOTROPIC_KL) .and. fpar%sim%USE_RHO_SCALING .and. (kernel_names(iker) == 'rhop')) then
        call log%write('This is scaling for rhop kernels...', .true.)
        this%ker_data_smooth(:,:,:,iker) = this%ker_data_smooth(:,:,:,2) * RHO_SCALING_FAC
      else
        call log%write('This is multi-grid smoothing of '//trim(kernel_names(iker))//' kernels...', .true.)
        call inv%sem2inv(this%ker_data(:,:,:,:,iker), gk)
        call inv%inv2grid(gk, gm)
        this%ker_data_smooth(:,:,:,iker) = gm
      endif
    end do
    call synchronize_all()
  end subroutine

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
    call log%write('This is tapering kernels...', .true.)

    ! dh = (distance_max_glob+distance_min_glob)/2.0_cr
    call tap%create(MEXT_V%x(1), MEXT_V%x(MEXT_V%nx),&
                    MEXT_V%y(1), MEXT_V%y(MEXT_V%ny),&
                    MEXT_V%z(1), MEXT_V%z(MEXT_V%nz), &
                    fpar%grid%regular_grid_interval(1),&
                    fpar%grid%regular_grid_interval(2),&
                    fpar%grid%regular_grid_interval(3),&
                    fpar%postproc%TAPER_H_SUPPRESS, &
                    fpar%postproc%TAPER_H_BUFFER, &
                    fpar%postproc%TAPER_V_SUPPRESS, &
                    fpar%postproc%TAPER_V_BUFFER)
    
    if (worldrank == 0) then
      do iker = 1, nkernel
      !   do i = 1, MEXT_V%nx
      !     do j = 1, MEXT_V%ny
      !       do k = 1, MEXT_V%nz
      !         val = tap%interp(MEXT_V%x(i), MEXT_V%y(j), MEXT_V%z(k))
      !         this%ker_data_smooth(i,j,k,iker) = this%ker_data_smooth(i,j,k,iker) * val
      !       enddo
      !     enddo
      !   enddo
        this%ker_data_smooth(:,:,:,iker) = this%ker_data_smooth(:,:,:,iker) * tap%taper
      enddo
    endif
    call sync_from_main_rank_cr_4d(this%ker_data_smooth, MEXT_V%nx, MEXT_V%ny, MEXT_V%nz, nkernel)
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
          call remove_event_kernel(ievt, 'cijkl_kernel')
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
    real(kind=cr) :: norm_val

    call log%write('This is taking sum of kernels for joint inversion...', .true.)

    if (worldrank == 0) then
      total_kernel = zeros(MEXT_V%nx, MEXT_V%ny, MEXT_V%nz, nkernel)
      do itype = 1, NUM_INV_TYPE
        if (.not. fpar%postproc%INV_TYPE(itype)) cycle
        type_name = INV_TYPE_NAMES(itype)
        call calc_kernel0_weight_grid(itype, norm_val)
        fname = trim(OPT_DIR)//'/gradient_'//trim(model_name)//'_'//trim(type_name)//'.h5'
        call read_grid_kernel_smooth(fname, kernel)
        total_kernel = total_kernel + fpar%postproc%JOINT_WEIGHT(itype)*kernel/norm_val
        write(msg, '(a,F18.6)') 'Max gradient of '//trim(type_name), &
              fpar%postproc%JOINT_WEIGHT(itype)*maxval(abs(kernel))/norm_val
        call log%write(msg, .false.)
      enddo
      fname = trim(OPT_DIR)//'/gradient_'//trim(model_name)//'.h5'
      call write_grid_kernel_smooth(total_kernel, fname)
    endif

  end subroutine sum_joint_kernel_grid

  subroutine calc_kernel0_weight_grid(itype, max_global)
    integer, intent(in) :: itype
    real(kind=cr), dimension(:,:,:,:), allocatable :: kernel_data
    real(kind=cr), intent(out) :: max_global
    character(len=MAX_STRING_LEN) :: fname

    fname = trim(OPT_DIR)//'/gradient_M00_'//trim(INV_TYPE_NAMES(itype))//'.h5'
    call read_grid_kernel_smooth(fname, kernel_data)
    max_global = maxval( abs(kernel_data))

  end subroutine calc_kernel0_weight_grid

  subroutine finalize(this)
    class(PostFlow), intent(inout) :: this
    call log%write('This is finalizing post-processing of '//trim(simu_type)//'...', .true.)
    call log%write('*******************************************', .false.)
    call free_shm_array(this%ks_win)
  end subroutine finalize

end module post_processing