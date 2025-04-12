module optimize_grid
  use config
  use fwat_mpi
  use input_params, fpar => fwat_par_global
  use opt_io
  use kernel_io, only: read_mesh_databases_minimum
  use logger, only: log
  use line_search
  use utils, only: zeros
  use model_grid_data
  use shared_input_parameters, only: TOMOGRAPHY_PATH
  use hdf5_interface
  use common_lib, only: get_kernel_names


  implicit none


  character(len=MAX_STRING_LEN) :: msg
  integer, private :: ier

  type :: OptGridFlow
    real(kind=cr), dimension(:,:,:,:), allocatable :: model, model_tmp, gradient, direction
    integer :: iter_current, iter_prev, iter_next
    character(len=MAX_STRING_LEN) :: output_model_path, model_fname
    contains 
      procedure :: init => init_optimize, model_update, get_SD_direction, get_lbfgs_direction, run_linesearch
      procedure, private :: get_model_idx, model_update_tmp, interp_initial_model
  end type OptGridFlow

contains
  subroutine init_optimize(this, is_read_database)
    class(OptGridFlow), intent(inout) :: this
    logical, intent(in) :: is_read_database
  
    this%model_fname = trim(TOMOGRAPHY_PATH)//'/tomography_model.h5'

    call read_mesh_databases_minimum(is_read_database)

    call create_grid()

    call log%init('output_fwat_optimize_'//trim(model_name)//'.log')

    call get_kernel_names()
    step_len = fpar%update%MAX_SLEN

    call this%get_model_idx()
    this%output_model_path = trim(OPT_DIR)//'/model_'//trim(model_next)//'.h5'

    if (this%iter_current == 0) then
      call this%interp_initial_model()
      call write_grid_model(trim(OPT_DIR)//'/model_M00.h5', this%model)
    else
      call read_model_grid(this%iter_current, this%model)
    endif
    call read_gradient_grid(this%iter_current, this%gradient)
    call synchronize_all()
  end subroutine init_optimize

  subroutine interp_initial_model(this)
    class(OptGridFlow), intent(inout) :: this
    real(kind=cr), dimension(:), allocatable :: x, y, z
    real(kind=cr), dimension(:,:,:), allocatable :: rm, gm
    integer :: ipar, i, j, k, nx, ny, nz
    type(hdf5_file) :: h5file

    if(worldrank == 0) then
      this%model = zeros(MEXT_V%nx, MEXT_V%ny, MEXT_V%nz, nkernel)
      
      call h5file%open(this%model_fname, status='old', action='read')
      call h5file%get('/x', x)
      call h5file%get('/y', y)
      call h5file%get('/z', z)
      nx = size(x)
      ny = size(y)
      nz = size(z)
      do ipar = 1, nkernel
        call h5file%get('/'//trim(parameter_names(ipar)), gm)
        gm = transpose_3(gm)
        do i = 1, MEXT_V%nx
          do j = 1, MEXT_V%ny
            do k = 1, MEXT_V%nz
              if (MEXT_V%x(i) < x(1) .or. MEXT_V%x(i) > x(nx) .or. &
                  MEXT_V%y(j) < y(1) .or. MEXT_V%y(j) > y(ny) .or. &
                  MEXT_V%z(k) < z(1) .or. MEXT_V%z(k) > z(nz)) then
                this%model(i,j,k,ipar) = interp3_nearest_simple(x, y, z, gm, MEXT_V%x(i), MEXT_V%y(j), MEXT_V%z(k))
              else
                this%model(i,j,k,ipar) = interp3(x, y, z, gm, MEXT_V%x(i), MEXT_V%y(j), MEXT_V%z(k))
              endif
            enddo
          enddo
        enddo
      enddo
      call h5file%close(finalize=.true.)
    endif

  end subroutine interp_initial_model

  subroutine get_model_idx(this)
    class(OptGridFlow), intent(inout) :: this

    ! get model index
    read(model_name(2:3),'(I2.2)', iostat=ier) this%iter_current
    if (ier /= 0) call exit_MPI(0, 'Error reading model name of '//trim(model_name))
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
      if (fpar%update%model_type == 1) then
      ! update model
        this%model = this%model * exp(step_len*this%direction)
        call alpha_scaling(this%model)
      elseif (fpar%update%model_type == 2) then
        this%model = this%model + step_len*this%direction
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
      this%model_tmp = zeros(MEXT_V%nx,MEXT_V%ny,MEXT_V%nz, nkernel)
      if (fpar%update%model_type == 1) then
        this%model_tmp = this%model * exp(step_len*this%direction)
        call alpha_scaling(this%model_tmp)
      elseif (fpar%update%model_type == 2) then
        this%model_tmp(:,:,:,ipar) = this%model(:,:,:,ipar) + step_len*this%direction(:,:,:,ipar)
      else
        call exit_MPI(0, 'Unknown model type')
      endif
    endif
    call synchronize_all()
    if (fpar%update%model_type == 1) then
    endif
    call synchronize_all()

  end subroutine model_update_tmp

  subroutine get_SD_direction(this)
    class(OptGridFlow), intent(inout) :: this
    real(kind=cr) :: max_dir

    call log%write('Starting steepest descent direction', .true.)
    if (worldrank == 0) then
      this%direction = -this%gradient
      max_dir = maxval(abs(this%direction))
      this%direction = this%direction / max_dir
    endif

    call synchronize_all()
  end subroutine get_SD_direction

  subroutine get_lbfgs_direction(this)
    class(OptGridFlow), intent(inout) :: this

    integer :: iter_store, istore
    real(kind=cr) :: max_dir_loc, max_dir
    real(kind=cr), dimension(:,:,:,:), allocatable :: model0, model1, gradient0, gradient1
    real(kind=cr), dimension(:,:,:,:), allocatable :: q_vector, r_vector, gradient_diff, model_diff
    real(kind=cr) :: p_sum, a_sum, b_sum, p(1000), a(1000), b, p_k_up_sum, p_k_down_sum, p_k

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

      do istore=this%iter_current-1,iter_store,-1
        call read_gradient_grid(istore+1,gradient1)
        call read_gradient_grid(istore,gradient0)
        call read_model_grid(istore+1,model1)
        call read_model_grid(istore,model0)

        gradient_diff = gradient1-gradient0
        model_diff = model1-model0

        ! call Parallel_ComputeInnerProduct(gradient_diff, model_diff, p_sum)
        p_sum = sum(gradient_diff*model_diff)
        p(istore) = 1.0_cr / p_sum
        ! call Parallel_ComputeInnerProduct(model_diff, q_vector, a_sum)
        a_sum = sum(model_diff*q_vector)
        a(istore) = p(istore)*a_sum
        q_vector = q_vector - a(istore)*gradient_diff
      enddo

      istore = this%iter_current - 1
      call read_gradient_grid(istore+1, gradient1)
      call read_gradient_grid(istore, gradient0)
      call read_model_grid(istore+1, model1)
      call read_model_grid(istore, model0)
      gradient_diff = gradient1 - gradient0
      model_diff = model1 - model0

      p_k_up_sum = sum(gradient_diff*model_diff)
      p_k_down_sum = sum(gradient_diff*gradient_diff)
      p_k = p_k_up_sum / p_k_down_sum
      r_vector=p_k*q_vector

      ! forward store
      do istore = iter_store, this%iter_current-1, 1
        call read_gradient_grid(istore+1, gradient1)
        call read_gradient_grid(istore, gradient0)
        call read_model_grid(istore+1, model1)
        call read_model_grid(istore, model0)

        gradient_diff = gradient1-gradient0
        model_diff = model1-model0

        ! call Parallel_ComputeInnerProduct(gradient_diff, r_vector, b_sum)
        b_sum = sum(gradient_diff*r_vector)
        b = p(istore)*b_sum
        r_vector = r_vector + (a(istore)-b)*model_diff
      enddo
      this%direction = -r_vector

      max_dir = maxval(abs(this%direction))
      this%direction = this%direction / max_dir
    endif ! worldrank == 0
    call synchronize_all()

  end subroutine get_lbfgs_direction

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
    use generate_databases_subs
    class(OptGridFlow), intent(inout) :: this
    real(kind=dp), dimension(NUM_INV_TYPE) :: total_misfit, misfit_start, misfit_prev
    integer :: itype, isub

    run_mode = FORWARD_MEASADJ
    do isub = 1, fpar%update%MAX_SUB_ITER
      total_misfit = 0.0_dp
      misfit_start = 0.0_dp
      misfit_prev = 0.0_dp
      write(msg, '("Starting ",I0,"th sub-iteration with step_length: ",F10.8)') isub, step_len
      call log%write(msg, .true.)
      call this%model_update_tmp()
      call write_grid_model(this%model_fname, this%model_tmp)
      call generate_databases_fwat()
      call synchronize_all()
      do itype = 1, NUM_INV_TYPE
        if (.not. fpar%postproc%INV_TYPE(itype)) cycle

        simu_type = INV_TYPE_NAMES(itype)
        
        call forward_for_simu_type(total_misfit(itype), misfit_start(itype), misfit_prev(itype))

        total_misfit(itype) = fpar%postproc%JOINT_WEIGHT(itype)*total_misfit(itype)/misfit_start(itype)
        misfit_prev(itype) = fpar%postproc%JOINT_WEIGHT(itype)*misfit_prev(itype)/misfit_start(itype)
      enddo
      call synchronize_all()
      
      if (sum(total_misfit) < sum(misfit_prev)) then
        if (worldrank == 0) then
          write(msg, '("Sum of misfit reduced from ",F22.8," to ",F22.8)') sum(misfit_prev), sum(total_misfit)
          call log%write(msg, .true.)
        endif
        exit
      else
        if (myrank == 0) then
          write(msg, '("Sum of misfit increased from ",F22.8," to ",F22.8)') sum(misfit_prev), sum(total_misfit)
          call log%write(msg, .true.)
        endif
        step_len = step_len * fpar%update%MAX_SHRINK
      endif
    enddo

  end subroutine run_linesearch
end module optimize_grid