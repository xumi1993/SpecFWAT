module optimize
  use config
  use fwat_mpi
  use input_params, fpar => fwat_par_global
  use opt_io
  use kernel_io, only: read_mesh_databases_minimum

  implicit none

  type :: OptFlow
    real(kind=cr), dimension(:,:,:,:,:), allocatable :: model, gradient, direction
    integer :: iter_current, iter_prev, iter_next
    character(len=MAX_STRING_LEN) :: model_next, model_prev, output_model_path
    contains
      procedure :: init, model_update, get_SD_direction, get_lbfgs_direction
      procedure, private :: get_model_idx, alpha_scaling
  end type OptFlow

contains
  subroutine init(this, is_read_database)
    class(OptFlow), intent(inout) :: this
    logical, intent(in) :: is_read_database
    
    call read_mesh_databases_minimum(is_read_database)

    if (fpar%update%model_type == 1) then
      nkernel = size(KERNEL_ISO)
      kernel_names = KERNEL_ISO
      parameter_names = MODEL_ISO
    elseif (fpar%update%model_type == 2) then
      nkernel = size(KERNEL_AZI_ANI)
      kernel_names = KERNEL_AZI_ANI
      parameter_names = MODEL_AZI_ANI
    else
      call exit_MPI(0, 'Unknown model type')
    endif
    step_len = fpar%update%MAX_SLEN

    call this%get_model_idx()
    this%output_model_path = trim(OPT_DIR)//'/MODEL_'//trim(this%model_next)
    if (worldrank == 0) then
      if (this%iter_current == 0) then
        call system('mkdir -p '//trim(OPT_DIR)//'/MODEL_M00 && cp model_initial/* '&
                    //trim(OPT_DIR)//'/MODEL_M00/')
      endif
      call system('mkdir -p '//trim(this%output_model_path))
    endif
    call synchronize_all()
    call read_model(this%iter_current, this%model)
    call read_gradient(this%iter_current, this%gradient)

  end subroutine init

  subroutine get_model_idx(this)
    class(OptFlow), intent(inout) :: this

    ! get model index
    read(model_name(2:),'(I2.2)') this%iter_current

    ! get model prev and next
    this%iter_next = this%iter_current+1
    this%iter_prev = this%iter_current-1

    write(this%model_next,'(A1,I2.2)') 'M', this%iter_next
    if (this%iter_prev < fpar%update%ITER_START) then
      this%model_prev='none'
    else
      write(this%model_prev,'(A1,I2.2)') 'M',this%iter_prev
    endif
  end subroutine get_model_idx

  subroutine model_update(this)
    class(OptFlow), intent(inout) :: this
    real(kind=cr) :: max_dir_loc, max_dir
    integer :: ipar

    ! normalize direction
    max_dir_loc = maxval(abs(this%direction))
    call max_all_all_cr(max_dir_loc, max_dir)
    this%direction = this%direction / max_dir

    if (fpar%update%model_type == 1) then
    ! update model
      do ipar = 1, nkernel
        this%model(:,:,:,:,ipar) = this%model(:,:,:,:,ipar) * exp(step_len*this%direction(:,:,:,:,ipar))
      enddo
      call this%alpha_scaling()
    elseif (fpar%update%model_type == 2) then
      do ipar = 1, nkernel
        this%model(:,:,:,:,ipar) = this%model(:,:,:,:,ipar) + step_len*this%direction(:,:,:,:,ipar)
      enddo
    else
      call exit_MPI(0, 'Unknown model type')
    endif
    call synchronize_all()
  end subroutine model_update

  subroutine get_SD_direction(this)
    class(OptFlow), intent(inout) :: this

    this%direction = this%gradient
  end subroutine get_SD_direction

  subroutine get_lbfgs_direction(this)
    class(OptFlow), intent(inout) :: this

    integer :: iter_store, istore
    real(kind=cr), dimension(:,:,:,:,:), allocatable :: model0, model1, gradient0, gradient1
    real(kind=cr), dimension(:,:,:,:,:), allocatable :: q_vector, r_vector, gradient_diff, model_diff
    real(kind=cr) :: p_sum, a_sum, b_sum, p(1000), a(1000), b, p_k_up_sum, p_k_down_sum, p_k

    ! get direction
    if (this%iter_current == fpar%update%ITER_START) then
      this%direction = this%gradient
      return
    endif

    iter_store = this%iter_current-fpar%update%LBFGS_M_STORE
    if ( iter_store <= fpar%update%ITER_START ) then
          iter_store = fpar%update%ITER_START
    endif

    call read_gradient(this%iter_current, q_vector)

    do istore=this%iter_current-1,iter_store,-1
      call read_gradient(istore+1,gradient1)
      call read_gradient(istore,gradient0)
      call read_model(istore+1,model1)
      call read_model(istore,model0)

      gradient_diff = gradient1-gradient0
      model_diff = model1-model0

      call Parallel_ComputeInnerProduct(gradient_diff, model_diff, p_sum)
      p(istore) = 1.0_cr / p_sum
      call Parallel_ComputeInnerProduct(model_diff, q_vector, a_sum)
      a(istore) = p(istore)*a_sum
      q_vector = q_vector - a(istore)*gradient_diff
    enddo

    istore = this%iter_current - 1
    call read_gradient(istore+1, gradient1)
    call read_gradient(istore, gradient0)
    call read_model(istore+1, model1)
    call read_model(istore, model0)
    gradient_diff = gradient1 - gradient0
    model_diff = model1 - model0

    call Parallel_ComputeInnerProduct(gradient_diff, model_diff, p_k_up_sum)
    call Parallel_ComputeInnerProduct(gradient_diff, gradient_diff, p_k_down_sum)
    p_k = p_k_up_sum / p_k_down_sum
    r_vector=p_k*q_vector

    ! forward store
    do istore = iter_store, this%iter_current-1, 1
      call read_gradient(istore+1, gradient1)
      call read_gradient(istore, gradient0)
      call read_model(istore+1, model1)
      call read_model(istore, model0)

      gradient_diff = gradient1-gradient0
      model_diff = model1-model0

      call Parallel_ComputeInnerProduct(gradient_diff, r_vector, b_sum)
      b = p(istore)*b_sum
      r_vector = r_vector + (a(istore)-b)*model_diff
    enddo
    this%direction = -r_vector
    call synchronize_all()

  end subroutine get_lbfgs_direction

  subroutine alpha_scaling(this)
    class(OptFlow), intent(inout) :: this

    where (this%model(:,:,:,:,1) < fpar%update%VPVS_RATIO_RANGE(1))
      this%model(:,:,:,:,1) = fpar%update%VPVS_RATIO_RANGE(1)
    end where
    where (this%model(:,:,:,:,1) > fpar%update%VPVS_RATIO_RANGE(2))
      this%model(:,:,:,:,1) = fpar%update%VPVS_RATIO_RANGE(2)
    end where
  end subroutine alpha_scaling
end module optimize