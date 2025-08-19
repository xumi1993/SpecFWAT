module optimize_grid
  use config
  use fwat_mpi
  use input_params, fpar => fwat_par_global
  use opt_io
  use logger, only: log
  use line_search
  use utils, only: zeros, ones, inner_product
  use model_grid_data
  use shared_input_parameters, only: TOMOGRAPHY_PATH
  use hdf5_interface
  use common_lib, only: get_kernel_names
  implicit none

  character(len=MAX_STRING_LEN), private :: msg
  integer, private :: ier

  type :: OptGridFlow
    real(kind=cr), dimension(:,:,:,:), allocatable :: model, model_iso, model_tmp, gradient, direction, hess
    integer :: iter_current, iter_prev, iter_next
    real(kind=cr) :: angle
    character(len=MAX_STRING_LEN) :: output_model_path, model_fname
    contains 
      procedure :: init => init_optimize, model_update, get_SD_direction, get_lbfgs_direction, get_CG_direction, run_linesearch
      procedure, private :: get_model_idx, model_update_tmp, interp_initial_model, read_hess_inv
  end type OptGridFlow

contains
  subroutine init_optimize(this)
    class(OptGridFlow), intent(inout) :: this
  
    this%model_fname = trim(TOMOGRAPHY_PATH)//'/tomography_model.h5'

    if ( fpar%update%MODEL_TYPE > 1) then
      ANISOTROPY = .true.
      ANISOTROPIC_KL = .true. 
    else
      ANISOTROPY = .false.
      ANISOTROPIC_KL = .false.
    end if

    call create_grid()

    call log%init('output_optimize_'//trim(model_name)//'.log')

    call get_kernel_names()
    step_len = fpar%update%MAX_SLEN

    call this%get_model_idx()

    if (.not. fpar%postproc%IS_PRECOND) then
      call this%read_hess_inv()
    else
      this%hess = ones(ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel)
    endif

    this%output_model_path = trim(OPT_DIR)//'/model_'//trim(model_next)//'.h5'

    if (this%iter_current == 0) then
      call this%interp_initial_model()
      call write_grid_model(trim(OPT_DIR)//'/model_M00.h5', this%model)
    else
      call read_model_grid(this%iter_current, this%model)
    endif
    if (this%iter_current < fpar%update%ITER_START) then
      call log%write('ERROR: Iteration '//trim(model_current)//&
                     ' is less than ITER_START', .true.)
      call exit_MPI(0, 'ERROR: Iteration '//trim(model_current)//&
                    ' is less than ITER_START')
    endif
    call read_gradient_grid(this%iter_current, this%gradient)
    call synchronize_all()
  end subroutine init_optimize

  subroutine interp_initial_model(this)
    class(OptGridFlow), intent(inout) :: this
    real(kind=cr), dimension(:), allocatable :: x, y, z
    real(kind=cr), dimension(:,:,:), allocatable :: rm, gm
    type(hdf5_file) :: h5file
    integer :: ipar
    character(len=256), dimension(:), allocatable :: keys

    if(worldrank == 0) then
      this%model = zeros(ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel)      
      call h5file%open(fpar%update%INIT_MODEL_PATH, status='old', action='read')

      call h5file%get('/x', x)
      call h5file%get('/y', y)
      call h5file%get('/z', z)

      do ipar = 1, nkernel
        call h5file%get('/'//trim(parameter_names(ipar)), gm)
        gm = transpose_3(gm)
        call model_interpolation(x, y, z, gm, rm)
        this%model(:,:,:,ipar) = rm
      enddo
      call h5file%close(finalize=.true.)
    endif

  end subroutine interp_initial_model

  subroutine read_hess_inv(this)
    class(OptGridFlow), intent(inout) :: this
    real(kind=cr), dimension(:,:,:), allocatable :: hess_loc
    integer :: i
    
    if(worldrank == 0) then
      call read_grid(trim(OPT_DIR)//'/'//trim(HESS_PREFIX)//'_'//trim(model_current)//'.h5', &
                      HESS_PREFIX, hess_loc)
      this%hess = zeros_dp(ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel)
      do i = 1, nkernel
        this%hess(:,:,:,i) = hess_loc
      enddo
    endif
    call synchronize_all()
  end subroutine read_hess_inv

  subroutine get_model_idx(this)
    class(OptGridFlow), intent(inout) :: this

    ! get model index
    read(model_name(2:3),'(I2.2)', iostat=ier) this%iter_current
    if (ier /= 0) call exit_MPI(0, 'Error reading model name of '//trim(model_name))
    if (this%iter_current < fpar%update%ITER_START) &
      call exit_MPI(0, 'ERROR: Iteration '//trim(model_name)//' is less than ITER_START')
    model_current = trim(model_name)
    model_name = trim(model_name)//'_ls'

    ! get model prev and next
    this%iter_next = this%iter_current+1
    this%iter_prev = this%iter_current-1

    write(model_next,'(A1,I2.2)') 'M', this%iter_next
    write(model_start,'(A1,I2.2)') 'M', fpar%update%ITER_START
    if (this%iter_prev < fpar%update%ITER_START) then
      model_prev='none'
    else
      write(model_prev,'(A1,I2.2)') 'M',this%iter_prev
    endif
  end subroutine get_model_idx

  subroutine model_update(this)
    class(OptGridFlow), intent(inout) :: this
    real(kind=cr) :: max_dir_loc, max_dir
    integer :: ipar

    if (worldrank == 0) then
      if (is_output_direction) then
        call write_grid_model(trim(OPT_DIR)//'/direction_'//trim(model_current)//'.h5', this%direction)
      endif
      write(msg, '(a,F10.8)') 'Update model parameter with step length: ', step_len
      call log%write(msg, .true.)
      if (fpar%update%model_type == 1) then
      ! update model
        this%model = this%model * exp(step_len*this%direction)
        call alpha_scaling(this%model)
      elseif (fpar%update%model_type == 2) then
        this%model(:,:,:,1:3) = this%model(:,:,:,1:3) * exp(step_len*this%direction(:,:,:,1:3))
        call alpha_scaling(this%model)
        this%model(:,:,:,4:5) = this%model(:,:,:,4:5) + step_len*this%direction(:,:,:,4:5)
      else
        call exit_MPI(0, 'Unknown model type')
      endif
    endif
    call synchronize_all()
  end subroutine model_update

  subroutine model_update_tmp(this)
    class(OptGridFlow), intent(inout) :: this
    integer :: ipar

    if (worldrank == 0) then
      this%model_tmp = zeros(ext_grid%nx,ext_grid%ny,ext_grid%nz, nkernel)
      if (fpar%update%model_type == 1) then
        this%model_tmp = this%model * exp(step_len*this%direction)
        call alpha_scaling(this%model_tmp)
      elseif (fpar%update%model_type == 2) then
        this%model_tmp(:,:,:,1:3) = this%model(:,:,:,1:3) * exp(step_len*this%direction(:,:,:,1:3))
        call alpha_scaling(this%model_tmp)
        this%model_tmp(:,:,:,4:5) = this%model(:,:,:,4:5) + step_len*this%direction(:,:,:,4:5)
      else
        call exit_MPI(0, 'Unknown model type')
      endif
    endif
    call synchronize_all()

  end subroutine model_update_tmp

  subroutine get_SD_direction(this)
    class(OptGridFlow), intent(inout) :: this
    real(kind=cr) :: max_dir

    call log%write('Starting steepest descent direction', .true.)
    if (worldrank == 0) then
      this%direction = -this%gradient * this%hess
      max_dir = maxval(abs(this%direction))
      this%direction = this%direction / max_dir
    endif

    call synchronize_all()
  end subroutine get_SD_direction

  subroutine get_lbfgs_direction(this)
    class(OptGridFlow), intent(inout) :: this

    integer :: iter_store, istore
    real(kind=cr) :: max_dir_loc, max_dir
    real(kind=cr), dimension(:,:,:,:), allocatable :: model0, model1, gradient0, gradient1, grad_bak
    real(kind=cr), dimension(:,:,:,:), allocatable :: q_vector, r_vector, gradient_diff, model_diff
    real(kind=cr) :: p_sum, a_sum, b_sum, p(1000), a(1000), b, p_k_up_sum, p_k_down_sum, p_k, angle

    ! get direction
    if (this%iter_current == fpar%update%ITER_START) then
      call this%get_SD_direction()
      return
    endif

    call log%write('Starting L-BFGS direction', .true.)

    if (worldrank == 0) then
      iter_store = this%iter_current-fpar%update%LBFGS_M_STORE
      if ( iter_store <= fpar%update%ITER_START ) then
        iter_store = fpar%update%ITER_START
      endif

      call read_gradient_grid(this%iter_current, q_vector)
      grad_bak = q_vector

      do istore=this%iter_current-1,iter_store,-1
        call read_gradient_grid(istore+1,gradient1)
        call read_gradient_grid(istore,gradient0)
        call read_model_grid(istore+1,model1,.true.)
        call read_model_grid(istore,model0,.true.)

        gradient_diff = gradient1-gradient0
        model_diff = model1-model0

        ! call Parallel_ComputeInnerProduct(gradient_diff, model_diff, p_sum)
        ! p_sum = sum(gradient_diff*model_diff)
        p_sum = inner_product(gradient_diff, model_diff, ext_grid%dx, ext_grid%dy, ext_grid%dz)
        p(istore) = 1.0_cr / p_sum
        ! a_sum = sum(model_diff*q_vector)
        a_sum = inner_product(model_diff, q_vector, ext_grid%dx, ext_grid%dy, ext_grid%dz)
        a(istore) = p(istore)*a_sum
        q_vector = q_vector - a(istore)*gradient_diff
      enddo

      istore = this%iter_current - 1
      call read_gradient_grid(istore+1, gradient1)
      call read_gradient_grid(istore, gradient0)
      call read_model_grid(istore+1, model1,.true.)
      call read_model_grid(istore, model0,.true.)
      gradient_diff = gradient1 - gradient0
      model_diff = model1 - model0

      ! p_k_up_sum = sum(gradient_diff*model_diff)
      p_k_up_sum = inner_product(gradient_diff, model_diff, ext_grid%dx, ext_grid%dy, ext_grid%dz)
      ! p_k_down_sum = sum(gradient_diff*gradient_diff)
      p_k_down_sum = inner_product(gradient_diff, gradient_diff, ext_grid%dx, ext_grid%dy, ext_grid%dz)
      p_k = p_k_up_sum / p_k_down_sum
      r_vector=p_k * this%hess * q_vector

      ! forward store
      do istore = iter_store, this%iter_current-1, 1
        call read_gradient_grid(istore+1, gradient1)
        call read_gradient_grid(istore, gradient0)
        call read_model_grid(istore+1, model1,.true.)
        call read_model_grid(istore, model0,.true.)

        gradient_diff = gradient1-gradient0
        model_diff = model1-model0

        ! call Parallel_ComputeInnerProduct(gradient_diff, r_vector, b_sum)
        ! b_sum = sum(gradient_diff*r_vector)
        b_sum = inner_product(gradient_diff, r_vector, ext_grid%dx, ext_grid%dy, ext_grid%dz)
        b = p(istore)*b_sum
        r_vector = r_vector + (a(istore)-b)*model_diff
      enddo
      this%direction = -r_vector

      max_dir = maxval(abs(this%direction))
      this%direction = this%direction / max_dir
      
      ! check the angle between search direction and negative grad
      block
        real(kind=cr) :: grad_sum, grad_norm, direc_norm
        grad_bak = -grad_bak
        grad_sum = sum(grad_bak*this%direction)
        grad_norm = sqrt(sum(grad_bak*grad_bak))
        direc_norm = sqrt(sum(this%direction*this%direction))
        this%angle = acos(grad_sum/(grad_norm*direc_norm))
        this%angle = this%angle * rad2deg
      end block
      write(msg, '(a,F8.4,a)') 'Angle between search direction and negative gradient: ', this%angle, ' degrees'
      call log%write(msg, .true.)
    endif ! worldrank == 0
    call synchronize_all()
    call bcast_all_singlecr(this%angle)
  end subroutine get_lbfgs_direction

  subroutine get_CG_direction(this)
    class(OptGridFlow), intent(inout) :: this
    real(kind=cr) :: max_dir_loc, max_dir, alpha
    real(kind=cr), dimension(:,:,:,:), allocatable :: gradient0, gradient1, r_vector
    real(kind=cr), dimension(:), allocatable :: norm0, norm1, ratio
    integer :: i

    ! get direction
    if (this%iter_current == fpar%update%ITER_START) then
      call this%get_SD_direction()
      return
    endif

    call log%write('Starting conjugate gradient direction', .true.)
    if (worldrank == 0) then
      ! read gradient
      gradient1 = this%gradient
      call read_gradient_grid(this%iter_prev, gradient0)

      ! compute powell ratio
      norm0 = zeros(nkernel)
      norm1 = zeros(nkernel)
      ratio = zeros(nkernel)
      do i = 1, nkernel
        norm0(i) = sum(gradient0(:,:,:,i)*gradient0(:,:,:,i))
        norm1(i) = sum(gradient0(:,:,:,i)*gradient1(:,:,:,i))
        ratio(i) = norm1(i) / norm0(i)
      enddo

      ! check the ratio
      ! if (all(ratio > 0.2)) then
        ! call log%write('Powell ratio: ( > 0.2 then restart with steepest descent)', .true.)
        ! alpha = 0.0_cr
      ! else
        ! difference kernel/gradient
        ! length ( ( gamma_n - gamma_(n-1))^T * lambda_n )
      r_vector = zeros(ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel)
        ! PR-CG
      do i = 1, nkernel
        norm1(i) = sum((gradient1(:,:,:,i)-gradient0(:,:,:,i))*gradient1(:,:,:,i))
      enddo
      alpha = sum(norm1) / sum(norm0)
        ! if (alpha < 0.0_cr) alpha = 0.0_cr
      ! endif
      ! compute direction
      write(msg, '(a,f0.6)') 'Alpha: ', alpha
      call log%write(msg, .true.)
      r_vector = - gradient1 - alpha * gradient0
      this%direction = this%hess * r_vector
      max_dir = maxval(abs(this%direction))
      this%direction = this%direction / max_dir
    endif
    call synchronize_all()
  end subroutine get_CG_direction

  subroutine alpha_scaling(model_inout)
    real(kind=cr), dimension(:,:,:,:), intent(inout) :: model_inout
    real(kind=cr), dimension(:,:,:), allocatable :: ratio

    ratio = model_inout(:,:,:,1) / model_inout(:,:,:,2)

    where (ratio < fpar%update%VPVS_RATIO_RANGE(1))
      model_inout(:,:,:,1) = fpar%update%VPVS_RATIO_RANGE(1)*model_inout(:,:,:,2)
    end where
    where (ratio > fpar%update%VPVS_RATIO_RANGE(2))
      model_inout(:,:,:,1) = fpar%update%VPVS_RATIO_RANGE(2)*model_inout(:,:,:,2)
    end where
  end subroutine alpha_scaling

  subroutine run_linesearch(this)
    use meshfem3D_subs
    use generate_databases_subs
    class(OptGridFlow), intent(inout) :: this
    real(kind=dp), dimension(NUM_INV_TYPE) :: total_misfit, misfit_prev
    integer :: itype, isub
    logical :: break_flag

    run_mode = FORWARD_MEASADJ
    local_path_backup = trim(LOCAL_PATH)
    do isub = 1, fpar%update%MAX_SUB_ITER
      total_misfit = 0.0_dp
      misfit_prev = 0.0_dp
      write(msg, '("Starting ",I0,"th sub-iteration with step_length: ",F10.8)') isub, step_len
      call log%write(msg, .true.)
      call this%model_update_tmp()
      call write_grid_model(this%model_fname, this%model_tmp)
      call synchronize_all()
      do itype = 1, NUM_INV_TYPE
        if (.not. fpar%postproc%INV_TYPE(itype)) cycle
        ! setup simulation type
        simu_type = INV_TYPE_NAMES(itype)
        call fpar%select_simu_type()

        ! generate mesh and database for this simu_type
        call meshfem3D_fwat(fpar%sim%Mesh_Par_file)
        call generate_databases_fwat(.true.)
        
        ! Forward simulation and measure misfits
        call forward_for_simu_type(total_misfit(itype), misfit_prev(itype))
        write(msg, '(A,F20.6)') 'Total misfit for '//trim(simu_type)//': ', total_misfit(itype)
        call log%write(msg, .true.)
      enddo
      call synchronize_all()

      call backtracking(this%direction, total_misfit, misfit_prev, break_flag)

      if (break_flag) exit      
      
    enddo
    if (isub > fpar%update%MAX_SUB_ITER) then
      call log%write('ERROR: Reached maximum number of line search iterations', .true.)
      call log%write('Please try to use more sub-iterations or reset the L-BFGS', .true.)
      call exit_MPI(0, 'ERROR: Reached maximum number of line search iterations')
    endif

  end subroutine run_linesearch
end module optimize_grid