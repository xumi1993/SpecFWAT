module optimize_gll
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use config
  use fwat_mpi
  use input_params, only: fpar => fwat_par_global
  use common_lib, only: mkdir
  use logger, only: runlog => log
  implicit none
  private
  public :: OptGLLFlow, read_gll_vector, write_gll_vector, gll_model_path

  type :: OptGLLFlow
    real(kind=cr), allocatable :: model(:,:,:,:,:), gradient(:,:,:,:,:), direction(:,:,:,:,:)
    real(kind=cr), allocatable :: hess(:,:,:,:), trial(:,:,:,:,:)
    real(kind=dp), allocatable :: weights(:,:,:,:)
    integer :: iter
    logical :: converged = .false.
  contains
    procedure :: init, archive_model, get_direction, update_trial, run_linesearch, dot, finalize
    procedure, private :: history_pair, normalize_direction
  end type OptGLLFlow
contains
  function gll_model_path(iter) result(path)
    integer, intent(in) :: iter
    character(len=MAX_STRING_LEN) :: path, name
    write(name, '(A,I2.2)') 'M', iter
    path = trim(OPT_DIR)//'/model_'//trim(name)
  end function gll_model_path

  subroutine read_gll_vector(path, names, vector, nspec)
    character(len=*), intent(in) :: path, names(:)
    real(kind=cr), allocatable, intent(out) :: vector(:,:,:,:,:)
    integer, optional, intent(in) :: nspec
    character(len=MAX_STRING_LEN) :: prefix, filename
    integer :: ipar, unit, ios, nlocal
    integer(kind=8) :: nbytes, record_bytes, element_bytes
    logical :: exists

    write(prefix, '(A,I6.6,A)') trim(path)//'/proc', worldrank, '_'
    element_bytes = int(NGLLX,8)*NGLLY*NGLLZ*(storage_size(0.0_cr)/8)
    do ipar = 1, size(names)
      filename = trim(prefix)//trim(names(ipar))//'.bin'
      inquire(file=filename, exist=exists, size=nbytes, iostat=ios)
      if (ios /= 0 .or. .not. exists) call exit_MPI(worldrank, 'Missing GLL file: '//trim(filename))
      ! SPECFEM binary files contain one sequential record with 4-byte markers.
      record_bytes = nbytes - 8_8
      if (record_bytes <= 0 .or. mod(record_bytes, element_bytes) /= 0) &
        call exit_MPI(worldrank, 'Invalid GLL record size: '//trim(filename))
      nlocal = int(record_bytes / element_bytes)
      if (present(nspec)) then
        if (nlocal /= nspec) call exit_MPI(worldrank, 'GLL mesh size mismatch: '//trim(filename))
      endif
      if (.not. allocated(vector)) allocate(vector(NGLLX,NGLLY,NGLLZ,nlocal,size(names)))
      if (nlocal /= size(vector,4)) call exit_MPI(worldrank, 'Inconsistent GLL parameter sizes: '//trim(filename))
      open(newunit=unit, file=filename, status='old', action='read', form='unformatted', iostat=ios)
      if (ios /= 0) call exit_MPI(worldrank, 'Cannot open GLL file: '//trim(filename))
      read(unit, iostat=ios) vector(:,:,:,:,ipar)
      close(unit)
      if (ios /= 0) call exit_MPI(worldrank, 'Cannot read GLL file: '//trim(filename))
      if (.not. all(ieee_is_finite(vector(:,:,:,:,ipar)))) &
        call exit_MPI(worldrank, 'Non-finite GLL values: '//trim(filename))
    enddo
  end subroutine read_gll_vector

  subroutine write_gll_vector(path, names, vector)
    character(len=*), intent(in) :: path, names(:)
    real(kind=cr), intent(in) :: vector(:,:,:,:,:)
    character(len=MAX_STRING_LEN) :: prefix, filename
    integer :: ipar, unit, ios

    call mkdir(path)
    call synchronize_all()
    write(prefix, '(A,I6.6,A)') trim(path)//'/proc', worldrank, '_'
    do ipar = 1, size(names)
      filename = trim(prefix)//trim(names(ipar))//'.bin'
      open(newunit=unit, file=filename, status='replace', action='write', form='unformatted', iostat=ios)
      if (ios /= 0) call exit_MPI(worldrank, 'Cannot write GLL file: '//trim(filename))
      write(unit, iostat=ios) vector(:,:,:,:,ipar)
      close(unit)
      if (ios /= 0) call exit_MPI(worldrank, 'Failed writing GLL file: '//trim(filename))
    enddo
    call synchronize_all()
  end subroutine write_gll_vector

  subroutine archive_model(this, iter)
    use specfem_par, only: NSPEC_AB, kappastore, mustore, rhostore, IMODEL, ATTENUATION
    use external_model, only: EXTERNAL_QMU, EXTERNAL_QKAPPA
    class(OptGLLFlow), intent(inout) :: this
    integer, intent(in) :: iter
    real(kind=cr), allocatable :: h(:,:,:,:,:)

    if (IMODEL == IMODEL_USER_EXTERNAL) then
      ! Postprocessing has reloaded the original material arrays from the solver
      ! database, before attenuation rescaling. No mesh-file export is needed.
      if (any(rhostore <= 0.0_cr) .or. any(mustore <= 0.0_cr) .or. any(kappastore <= 0.0_cr)) &
        call exit_MPI(worldrank, 'Initial elastic density and moduli must be positive')
      allocate(this%model(NGLLX,NGLLY,NGLLZ,NSPEC_AB,3))
      this%model(:,:,:,:,1) = sqrt((kappastore + (4.0_cr/3.0_cr)*mustore)/rhostore)
      this%model(:,:,:,:,2) = sqrt(mustore/rhostore)
      this%model(:,:,:,:,3) = rhostore
    else
      ! Preserve the exact accepted GLL values used by this database build.
      call read_gll_vector(local_path_fwat, parameter_names, this%model, NSPEC_AB)
    endif
    if (.not. all(ieee_is_finite(this%model)) .or. any(this%model <= 0.0_cr)) &
      call exit_MPI(worldrank, 'Current GLL vp, vs and rho must be finite and positive')
    if (any(this%model(:,:,:,:,1)/this%model(:,:,:,:,2) < fpar%update%VPVS_RATIO_RANGE(1)) .or. &
        any(this%model(:,:,:,:,1)/this%model(:,:,:,:,2) > fpar%update%VPVS_RATIO_RANGE(2))) &
      call exit_MPI(worldrank, 'Current GLL Vp/Vs is outside VPVS_RATIO_RANGE')
    call write_gll_vector(gll_model_path(iter), parameter_names, this%model)
    if (IMODEL == IMODEL_USER_EXTERNAL) &
      call write_gll_vector(local_path_fwat, parameter_names, this%model)

    ! Export only the fixed quality factors needed by later GLL database builds.
    if (ATTENUATION .and. iter == 0) then
      if (IMODEL == IMODEL_USER_EXTERNAL) then
        allocate(h(NGLLX,NGLLY,NGLLZ,NSPEC_AB,2))
        h(:,:,:,:,1) = EXTERNAL_QMU
        h(:,:,:,:,2) = EXTERNAL_QKAPPA
        call write_gll_vector(local_path_fwat, ['qmu   ', 'qkappa'], h)
      else
        call read_gll_vector(local_path_fwat, ['qmu   ', 'qkappa'], h, NSPEC_AB)
      endif
      call write_gll_vector(gll_model_path(0), ['qmu   ', 'qkappa'], h)
      deallocate(h)
    endif
  end subroutine archive_model

  subroutine init(this, iter)
    use specfem_par, only: NSPEC_AB, irregular_element_number, jacobian_regular, jacobianstore, wxgll, wygll, wzgll
    class(OptGLLFlow), intent(inout) :: this
    integer, intent(in) :: iter
    real(kind=cr), allocatable :: h(:,:,:,:,:)
    integer :: i, j, k, ispec, irreg
    real(kind=dp) :: jac
    character(len=MAX_STRING_LEN) :: path

    this%iter = iter
    this%converged = .false.
    call this%archive_model(iter)
    path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)
    call read_gll_vector(path, trim_names('_kernel_smooth'), this%gradient, NSPEC_AB)
    if (allocated(this%hess)) deallocate(this%hess, this%weights)
    allocate(this%hess(NGLLX,NGLLY,NGLLZ,NSPEC_AB), this%weights(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    this%hess = 1.0_cr
    if (.not. fpar%postproc%IS_HESS_PRECOND) then
      call read_gll_vector(path, [HESS_PREFIX], h, NSPEC_AB)
      this%hess = h(:,:,:,:,1)
    endif
    call setup_GLL_points()
    do ispec = 1, NSPEC_AB
      irreg = irregular_element_number(ispec)
      do k = 1, NGLLZ; do j = 1, NGLLY; do i = 1, NGLLX
        jac = jacobian_regular
        if (irreg /= 0) jac = jacobianstore(i,j,k,irreg)
        this%weights(i,j,k,ispec) = jac*wxgll(i)*wygll(j)*wzgll(k)
      enddo; enddo; enddo
    enddo
  end subroutine init

  function trim_names(suffix) result(names)
    character(len=*), intent(in) :: suffix
    character(len=MAX_STRING_LEN) :: names(nkernel)
    integer :: i
    do i = 1, nkernel
      names(i) = trim(kernel_names(i))//suffix
    enddo
  end function trim_names

  real(kind=dp) function dot(this, a, b) result(value)
    class(OptGLLFlow), intent(in) :: this
    real(kind=cr), intent(in) :: a(:,:,:,:,:), b(:,:,:,:,:)
    real(kind=dp) :: local_sum
    integer :: ipar
    local_sum = 0.0_dp
    do ipar = 1, size(a,5)
      local_sum = local_sum + sum(dble(a(:,:,:,:,ipar))*dble(b(:,:,:,:,ipar))*this%weights)
    enddo
    call sum_all_all_dp(local_sum, value)
  end function dot

  subroutine history_pair(this, iter, s, y)
    class(OptGLLFlow), intent(in) :: this
    integer, intent(in) :: iter
    real(kind=cr), allocatable, intent(out) :: s(:,:,:,:,:), y(:,:,:,:,:)
    real(kind=cr), allocatable :: a(:,:,:,:,:), b(:,:,:,:,:)
    character(len=MAX_STRING_LEN) :: name
    integer :: idx

    do idx = iter, iter+1
      call read_gll_vector(gll_model_path(idx), parameter_names, a, size(this%model,4))
      if (any(a <= 0.0_cr)) call exit_MPI(worldrank, 'Non-positive model in GLL L-BFGS history')
      write(name, '(A,I2.2)') 'M', idx
      call read_gll_vector(trim(OPT_DIR)//'/SUM_KERNELS_'//trim(name), &
                           trim_names('_kernel_smooth'), b, size(this%model,4))
      if (idx == iter) then
        s = -log(a)
        y = -b
      else
        s = s + log(a)
        y = y + b
      endif
    enddo
  end subroutine history_pair

  subroutine get_direction(this)
    class(OptGLLFlow), intent(inout) :: this
    real(kind=cr), allocatable :: q(:,:,:,:,:), s(:,:,:,:,:), y(:,:,:,:,:)
    real(kind=dp), allocatable :: rho(:), alpha(:)
    real(kind=dp) :: ys, yy, ss, gamma, beta
    integer :: first, idx, ipar
    character(len=MAX_STRING_LEN) :: msg

    q = this%gradient
    gamma = 1.0_dp
    first = max(fpar%update%ITER_START, this%iter-fpar%update%LBFGS_M_STORE)
    allocate(rho(first:this%iter-1), alpha(first:this%iter-1))
    rho = 0.0_dp
    alpha = 0.0_dp
    if (fpar%update%OPT_METHOD == 2) then
      do idx = this%iter-1, first, -1
        call this%history_pair(idx, s, y)
        ys = this%dot(y,s)
        yy = this%dot(y,y)
        ss = this%dot(s,s)
        if (ys <= 1.e-8_dp*sqrt(ss)*sqrt(yy) .or. yy <= tiny(yy)) then
          write(msg, '(A,I0)') 'Skipping L-BFGS pair with insufficient positive curvature: ', idx
          call runlog%write(msg)
          cycle
        endif
        if (all(rho == 0.0_dp)) gamma = ys/yy
        rho(idx) = 1.0_dp/ys
        alpha(idx) = rho(idx)*this%dot(s,q)
        q = q - real(alpha(idx),cr)*y
      enddo
    endif
    do ipar = 1, nkernel
      q(:,:,:,:,ipar) = real(gamma,cr)*this%hess*q(:,:,:,:,ipar)
    enddo
    if (fpar%update%OPT_METHOD == 2) then
      do idx = first, this%iter-1
        if (rho(idx) == 0.0_dp) cycle
        call this%history_pair(idx, s, y)
        beta = rho(idx)*this%dot(y,q)
        q = q + real(alpha(idx)-beta,cr)*s
      enddo
    endif
    this%direction = -q
    if (this%dot(this%gradient,this%direction) >= 0.0_dp) then
      call runlog%write('Restarting with preconditioned steepest descent')
      do ipar = 1, nkernel
        this%direction(:,:,:,:,ipar) = -this%hess*this%gradient(:,:,:,:,ipar)
      enddo
    endif
    call this%normalize_direction()
  end subroutine get_direction

  subroutine normalize_direction(this)
    class(OptGLLFlow), intent(inout) :: this
    real(kind=cr) :: local_max, global_max
    if (.not. all(ieee_is_finite(this%direction))) call exit_MPI(worldrank, 'Non-finite GLL search direction')
    local_max = maxval(abs(this%direction))
    call max_all_all_cr(local_max, global_max)
    this%converged = global_max <= tiny(global_max)
    if (.not. this%converged) this%direction = this%direction/global_max
  end subroutine normalize_direction

  subroutine update_trial(this)
    class(OptGLLFlow), intent(inout) :: this
    this%trial = this%model*exp(step_len*this%direction)
    this%trial(:,:,:,:,1) = max(fpar%update%VPVS_RATIO_RANGE(1)*this%trial(:,:,:,:,2), &
                          min(fpar%update%VPVS_RATIO_RANGE(2)*this%trial(:,:,:,:,2), this%trial(:,:,:,:,1)))
    if (.not. all(ieee_is_finite(this%trial)) .or. any(this%trial <= 0.0_cr)) &
      call exit_MPI(worldrank, 'Invalid updated GLL model')
  end subroutine update_trial

  subroutine run_linesearch(this, current_misfit)
    use generate_databases_subs, only: generate_databases_fwat
    use line_search, only: forward_for_simu_type
    use specfem_par, only: LOCAL_PATH, OUTPUT_FILES, SIMULATION_TYPE, SAVE_FORWARD
    class(OptGLLFlow), intent(inout) :: this
    real(kind=cr), allocatable :: raw(:,:,:,:,:), displacement(:,:,:,:,:)
    real(kind=dp), intent(in) :: current_misfit
    real(kind=dp) :: f1, f0, slope
    integer :: isub
    character(len=MAX_STRING_LEN) :: msg, base_name, base_output

    base_name = model_name
    base_output = OUTPUT_FILES
    model_current = base_name
    call read_gll_vector(trim(OPT_DIR)//'/SUM_KERNELS_'//trim(base_name), trim_names('_kernel'), raw, size(this%model,4))
    model_name = trim(base_name)//'_ls'
    run_mode = FORWARD_MEASADJ
    do isub = 1, fpar%update%MAX_SUB_ITER
      call this%update_trial()
      displacement = log(this%trial/this%model)
      slope = this%dot(raw,displacement)
      if (.not. ieee_is_finite(slope) .or. slope >= 0.0_dp) exit
      call write_gll_vector(local_path_fwat, parameter_names, this%trial)
      LOCAL_PATH = local_path_fwat
      OUTPUT_FILES = base_output
      SIMULATION_TYPE = 1
      SAVE_FORWARD = .true.
      call generate_databases_fwat()
      call forward_for_simu_type(f1, f0, current_misfit)
      write(msg, '(A,I0,A,ES12.4,A,ES18.8,A,ES18.8)') &
        'GLL line search ', isub, ': step=', step_len, ', misfit=', f1, ', previous=', f0
      call runlog%write(msg, .true.)
      if (ieee_is_finite(f1) .and. f1 <= f0 + fpar%update%C1*slope) then
        model_name = base_name
        LOCAL_PATH = local_path_fwat
        OUTPUT_FILES = base_output
        return
      endif
      step_len = step_len*fpar%update%MAX_SHRINK
    enddo
    ! Restore the accepted model and its solver database on search failure.
    call write_gll_vector(local_path_fwat, parameter_names, this%model)
    LOCAL_PATH = local_path_fwat
    OUTPUT_FILES = base_output
    SIMULATION_TYPE = 1
    SAVE_FORWARD = .true.
    call generate_databases_fwat()
    call exit_MPI(worldrank, 'GLL line search failed; current model restored. Reduce MAX_SLEN or reset ITER_START.')
  end subroutine run_linesearch
  subroutine finalize(this)
    class(OptGLLFlow), intent(inout) :: this
    if (allocated(this%model)) deallocate(this%model)
    if (allocated(this%gradient)) deallocate(this%gradient)
    if (allocated(this%direction)) deallocate(this%direction)
    if (allocated(this%hess)) deallocate(this%hess)
    if (allocated(this%weights)) deallocate(this%weights)
    if (allocated(this%trial)) deallocate(this%trial)
  end subroutine finalize
end module optimize_gll
