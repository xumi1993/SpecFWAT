module line_search
  use fwat_mpi
  use config
  use input_params, fpar => fwat_par_global
  use opt_io
  use preproc_fwd
  use misfit_mod, only : read_evt_misfit, rename_evt_misfits
  use common_lib, only : get_dat_type, run_shell, rename_path
  use logger, only: log
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  character(len=MAX_STRING_LEN), private :: msg
  ! Bracket of the weak Wolfe search: trial_lo at step_lo met the Armijo but
  ! not the curvature condition, and step_hi failed Armijo (< 0: none yet).
  real(kind=cr), private :: step_lo, step_hi
  integer, private :: trial_lo

contains
  subroutine init_step_len(iter)
    ! The line search starts from the step length accepted in the previous
    ! iteration (key M<iter-1> of optimize/step_len.h5), not larger than MAX_SLEN.
    integer, intent(in) :: iter
    real(kind=cr) :: slen_prev
    logical :: found

    step_len = fpar%update%MAX_SLEN
    if (.not. fpar%update%DO_LS .or. iter <= 0) return
    call read_step_len(iter - 1, slen_prev, found)
    if (found .and. slen_prev > 0.0_cr) then
      step_len = min(slen_prev, fpar%update%MAX_SLEN)
      write(msg, '(a,F10.8)') 'Line search starts from the step length accepted in the previous iteration: ', step_len
      call log%write(msg, .true.)
    endif
  end subroutine init_step_len

  subroutine forward_for_simu_type(total_misfit, misfit_prev, current_misfit)
    type(PrepareFWD) :: ffwd
    logical :: is_output_backup
    real(kind=dp), intent(out) :: total_misfit, misfit_prev
    real(kind=dp), optional, intent(in) :: current_misfit
    real(kind=dp) :: misfit_loc
    integer :: ievt

    is_output_backup = IS_OUTPUT_PREPROC
    IS_OUTPUT_PREPROC = .false.
    total_misfit = 0.0_dp
    misfit_prev = 0.0_dp
    if (present(current_misfit)) misfit_prev = current_misfit

    call get_dat_type()

    call fpar%acqui%read()

    ! initialize fwd
    call ffwd%init(.false.)

    do ievt = 1, fpar%acqui%nevents
      ffwd%ievt = ievt

      ! prepare simulation
      call ffwd%prepare_for_event()

      ! run forward simulation
      call ffwd%simulation()
      
      ! take sum of misfit
      total_misfit = total_misfit + ffwd%obj_func
      if (.not. present(current_misfit)) then
        misfit_loc = read_evt_misfit(model_current, ievt)
        misfit_prev = misfit_prev + misfit_loc
      endif

      call synchronize_all()
    enddo

    call ffwd%destroy()

    call fpar%acqui%finalize()

    IS_OUTPUT_PREPROC = is_output_backup
    call synchronize_all()
  end subroutine forward_for_simu_type

  function trial_name(isub) result(name)
    ! Each line-search trial keeps its solver outputs under its own model name.
    integer, intent(in) :: isub
    character(len=MAX_STRING_LEN) :: name

    write(name, '(A,"_ls",I2.2)') trim(model_current), isub
  end function trial_name

  subroutine wolfe_init()
    step_lo = 0.0_cr
    step_hi = -1.0_cr
    trial_lo = 0
  end subroutine wolfe_init

  subroutine wolfe_step(isub, f0, f1, g0, g1, accepted)
    ! Test the weak Wolfe conditions for trial isub at step_len, where g0 and
    ! g1 are the directional derivatives of the objective at steps 0 and
    ! step_len. Return the accepted trial, or choose the next step_len:
    !   Armijo fails:           shrink, or bisect once a lower bound exists;
    !   curvature fails:        expand up to MAX_SLEN, or bisect once an upper
    !                           bound exists;
    !   both hold:              accept.
    ! Kernels of trials that can no longer be accepted are removed.
    integer, intent(in) :: isub
    real(kind=dp), intent(in) :: f0, f1, g0, g1
    integer, intent(out) :: accepted
    logical :: armijo, curvature

    armijo = ieee_is_finite(f1) .and. f1 <= f0 + fpar%update%C1*step_len*g0
    curvature = g1 >= fpar%update%C2*g0
    call bcast_all_singlel(armijo)
    call bcast_all_singlel(curvature)
    write(msg, '(a,ES12.4,a,ES16.8,a,ES16.8,a,ES14.6,a,ES14.6)') 'Wolfe: step = ', step_len, &
      ', f1 = ', f1, ', f0+c1*step*g0 = ', f0 + fpar%update%C1*step_len*g0, &
      ', g1 = ', g1, ', c2*g0 = ', fpar%update%C2*g0
    call log%write(msg, .true.)

    accepted = 0
    if (.not. armijo) then
      call remove_trial(isub)
      step_hi = step_len
      if (trial_lo > 0) then
        step_len = 0.5_cr*(step_lo + step_hi)
      else
        step_len = step_len*fpar%update%MAX_SHRINK
      endif
      write(msg, '(a,F10.8)') 'Armijo condition not met; next step length: ', step_len
    elseif (curvature) then
      if (trial_lo > 0) call remove_trial(trial_lo)
      accepted = isub
      write(msg, '(a,F10.8)') 'Weak Wolfe conditions met; accept step length: ', step_len
    else
      if (trial_lo > 0) call remove_trial(trial_lo)
      step_lo = step_len
      trial_lo = isub
      if (step_hi > 0.0_cr) then
        step_len = 0.5_cr*(step_lo + step_hi)
      elseif (step_lo >= fpar%update%MAX_SLEN) then
        accepted = isub
      else
        step_len = min(step_len/fpar%update%MAX_SHRINK, fpar%update%MAX_SLEN)
      endif
      if (accepted > 0) then
        write(msg, '(a,F10.8)') 'Curvature condition not met at MAX_SLEN; accept step length: ', step_len
      else
        write(msg, '(a,F10.8)') 'Curvature condition not met; next step length: ', step_len
      endif
    endif
    call log%write(msg, .true.)
  end subroutine wolfe_step

  subroutine wolfe_fallback(accepted)
    ! When MAX_SUB_ITER is used up, accept the largest trial that satisfied
    ! the Armijo condition, if any.
    integer, intent(out) :: accepted

    accepted = trial_lo
    if (accepted > 0) then
      step_len = step_lo
      write(msg, '(a,F10.8)') 'MAX_SUB_ITER reached; accept the last step length satisfying Armijo: ', step_len
      call log%write(msg, .true.)
    endif
  end subroutine wolfe_fallback

  subroutine remove_trial(isub)
    ! Delete the event kernels and other solver outputs of a rejected trial.
    integer, intent(in) :: isub

    call run_shell('rm -rf '//trim(SOLVER_DIR)//'/'//trim(trial_name(isub))//'.*')
    call log%write('Removed solver outputs of '//trim(trial_name(isub)), .true.)
  end subroutine remove_trial

  subroutine adopt_trial(isub, next_name)
    ! The accepted trial already ran the forward and adjoint simulations of
    ! the next model, so rename its solver outputs and misfits to next_name.
    ! LS_KERNEL_FLAG lets xfwat_fwd_measure_adj -r 3 skip that model.
    integer, intent(in) :: isub
    character(len=*), intent(in) :: next_name
    character(len=MAX_STRING_LEN) :: trial, solver_next
    integer :: itype, ievt, unit

    trial = trial_name(isub)
    do itype = 1, NUM_INV_TYPE
      if (.not. fpar%postproc%INV_TYPE(itype)) cycle
      simu_type = INV_TYPE_NAMES(itype)
      LOCAL_PATH = local_path_backup
      call fpar%select_simu_type()
      LOCAL_PATH = local_path_fwat
      call get_dat_type()

      ! One solver directory per data type, renamed on the main rank.
      solver_next = trim(SOLVER_DIR)//'/'//trim(next_name)//'.'//trim(dat_type)
      call run_shell('rm -rf '//trim(solver_next))
      if (worldrank == 0) then
        call rename_path(trim(SOLVER_DIR)//'/'//trim(trial)//'.'//trim(dat_type), solver_next)
        open(newunit=unit, file=trim(solver_next)//'/'//trim(LS_KERNEL_FLAG), status='replace')
        close(unit)
      endif

      ! Misfit files of the events, distributed over the ranks.
      call fpar%acqui%read()
      do ievt = worldrank + 1, fpar%acqui%nevents, worldsize
        call rename_evt_misfits(trial, next_name, ievt)
      enddo
      call fpar%acqui%finalize()
      call synchronize_all()
    enddo
    call log%write('Solver outputs of '//trim(trial)//' are kept as '//trim(next_name), .true.)
  end subroutine adopt_trial

  subroutine joint_norm_weights(norm_val)
    ! Weight each data type as sum_joint_kernel_grid does, so the objective
    ! sum_i w_i*f_i/n_i matches the joint gradient sum_i w_i*g_i/n_i.
    use post_processing, only: calc_misfit0_weight, calc_kernel0_weight_grid
    real(kind=cr), intent(out) :: norm_val(NUM_INV_TYPE)
    real(kind=cr) :: std_val, linf_val
    integer :: itype

    if (.not. is_joint) then
      norm_val = 1.0_cr
      return
    endif
    norm_val = 0.0_cr
    do itype = 1, NUM_INV_TYPE
      if (.not. fpar%postproc%INV_TYPE(itype)) cycle
      if (fpar%postproc%NORM_TYPE == 2) then
        call calc_misfit0_weight(itype, norm_val(itype))
      elseif (worldrank == 0) then
        ! The reference gradient is read on the main rank only.
        call calc_kernel0_weight_grid(itype, linf_val, std_val)
        if (fpar%postproc%NORM_TYPE == 1) then
          norm_val(itype) = linf_val
        else
          norm_val(itype) = std_val
        endif
      endif
    enddo
    call bcast_all_cr(norm_val, NUM_INV_TYPE)
    where (fpar%postproc%INV_TYPE) norm_val = fpar%postproc%JOINT_WEIGHT / norm_val
    norm_val = norm_val / maxval(norm_val)
  end subroutine joint_norm_weights

  function grid_slope(direction, norm_val, from_trial) result(slope)
    ! Directional derivative sum_i w_i <K_i, d> of the weighted objective
    ! along the grid direction d, using the raw (unsmoothed) kernels K_i of
    ! the current model, or of the trial named model_name if from_trial.
    use kernel_io
    use post_processing, only: sum_event_kernels
    use model_grid_data, only: grid2gll

    real(kind=cr), dimension(:,:,:,:), intent(in) :: direction
    real(kind=cr), intent(in) :: norm_val(NUM_INV_TYPE)
    logical, intent(in) :: from_trial
    real(kind=dp) :: slope
    real(kind=cr) :: qp
    real(kind=cr), dimension(:,:,:,:,:), allocatable :: ker_sum
    real(kind=cr), dimension(:,:,:,:), allocatable :: ker, direc_gll
    real(kind=cr), dimension(:,:,:), allocatable :: direc_grid
    character(len=MAX_STRING_LEN) :: kernel_path
    integer :: itype, iker

    slope = 0.0_dp
    do itype = 1, NUM_INV_TYPE
      if (.not. fpar%postproc%INV_TYPE(itype)) cycle
      ! Each data type has its own mesh, on which its kernels were written.
      simu_type = INV_TYPE_NAMES(itype)
      LOCAL_PATH = local_path_backup
      call fpar%select_simu_type()
      LOCAL_PATH = local_path_fwat
      call read_mesh_databases_for_init()
      ! Quadrature weights of the inner product; no simulation may have set them yet.
      call setup_GLL_points()
      if (from_trial) then
        call get_dat_type()
        call fpar%acqui%read()
        allocate(ker_sum(NGLLX, NGLLY, NGLLZ, NSPEC_FWAT, nkernel))
        call sum_event_kernels(ker_sum)
        call fpar%acqui%finalize()
      elseif (is_joint) then
        kernel_path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_current)//'_'//trim(simu_type)
      else
        kernel_path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_current)
      endif

      do iker = 1, nkernel
        if (from_trial) then
          ker = ker_sum(:,:,:,:,iker)
        else
          call read_kernel(trim(kernel_path), trim(kernel_names(iker))//'_kernel', ker)
        endif
        direc_grid = zeros(ext_grid%nx, ext_grid%ny, ext_grid%nz)
        if (worldrank == 0) direc_grid = direction(:,:,:,iker)
        call bcast_all_cr(direc_grid, ext_grid%nx*ext_grid%ny*ext_grid%nz)
        call grid2gll(direc_grid, direc_gll)
        call Parallel_ComputeInnerProduct(direc_gll, ker, qp)
        slope = slope + qp * norm_val(itype)
        call synchronize_all()
      enddo
      if (allocated(ker_sum)) deallocate(ker_sum)
    enddo
  end function grid_slope
end module line_search