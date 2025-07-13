module line_search
  use fwat_mpi
  use config
  use input_params, fpar => fwat_par_global
  use opt_io
  use preproc_fwd
  use window_chi, only : read_model_misfit
  use common_lib, only : get_dat_type
  use logger, only: log

  implicit none

  character(len=MAX_STRING_LEN), private :: msg

contains
  subroutine forward_for_simu_type(total_misfit, misfit_prev)
    type(PrepareFWD) :: ffwd
    logical :: is_output_backup
    real(kind=dp), intent(out) :: total_misfit, misfit_prev
    real(kind=dp) :: misfit_loc
    integer :: ievt

    is_output_backup = IS_OUTPUT_PREPROC
    IS_OUTPUT_PREPROC = .false.
    total_misfit = 0.0_dp
    misfit_prev = 0.0_dp

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

      total_misfit = total_misfit + ffwd%obj_func

      if (model_prev /= 'none') then
        call read_model_misfit(model_current, ievt, misfit_loc)
        misfit_prev = misfit_prev + misfit_loc
      else
        misfit_prev = total_misfit
      endif

      call synchronize_all()
    enddo

    call ffwd%destroy()

    call fpar%acqui%finalize()

    IS_OUTPUT_PREPROC = is_output_backup
    call synchronize_all()
  end subroutine forward_for_simu_type

  subroutine backtracking(direction, misfit_ls, misfit_prev, break_flag)
    use kernel_io
    use post_processing, only: calc_misfit0_weight, calc_kernel0_weight_grid
    use model_grid_data, only: gll2grid

    real(kind=cr), dimension(:,:,:,:), intent(in) :: direction
    real(kind=dp), dimension(NUM_INV_TYPE), intent(in) :: misfit_prev, misfit_ls
    logical, intent(out) :: break_flag
    real(kind=dp) :: f1, f0, qt
    real(kind=cr) :: norm_val(NUM_INV_TYPE), std_val, linf_val, mean_val
    real(kind=cr), dimension(:,:,:,:), allocatable :: ker, q0
    real(kind=cr), dimension(:,:,:), allocatable :: gker
    character(len=MAX_STRING_LEN) :: kernel_path
    integer :: itype, iker

    ! calculate norm scale
    if (is_joint) then
      norm_val = 0.0_cr
      do itype = 1, NUM_INV_TYPE
        if (.not. fpar%postproc%INV_TYPE(itype)) cycle
        if (fpar%postproc%NORM_TYPE == 2) then
          call calc_misfit0_weight(itype, norm_val(itype))
        else
          call calc_kernel0_weight_grid(itype, linf_val, std_val)
          if (worldrank == 0) then
            if (fpar%postproc%NORM_TYPE == 1) then
              norm_val(itype) = linf_val
            else
              norm_val(itype) = std_val
            endif
          endif
        endif
        norm_val(itype) = fpar%postproc%JOINT_WEIGHT(itype) / norm_val(itype)
      enddo
      call bcast_all_cr(norm_val, NUM_INV_TYPE)
      norm_val = norm_val / maxval(norm_val)
    else
      norm_val = 1.0_cr
    endif

    ! calculate f0 and f1
    f0 = 0.0_dp
    f1 = 0.0_dp
    kernel_path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_current)
    if (worldrank == 0) q0 = zeros(ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel)
    do itype = 1, NUM_INV_TYPE
      if (.not. fpar%postproc%INV_TYPE(itype)) cycle
      f0 = f0 + misfit_prev(itype) * norm_val(itype)
      f1 = f1 + misfit_ls(itype) * norm_val(itype)
  
      ! read gradient and project to regular grid
      if (is_joint) then
        kernel_path = trim(kernel_path)//'_'//trim(INV_TYPE_NAMES(itype))
      endif
      do iker = 1, nkernel
        call read_kernel(trim(kernel_path), trim(kernel_names(iker))//'_kernel', ker)
        call gll2grid(ker, gker)
        if (worldrank == 0) q0(:,:,:,iker) = q0(:,:,:,iker) + gker * norm_val(itype)
        call synchronize_all()
      enddo
    enddo

    if (worldrank == 0) then
      qt = sum(q0*direction*ext_grid%dx*ext_grid%dy*ext_grid%dz)*step_len*fpar%update%C1
      write(msg, '(a,f20.6,a,f20.6,a,f20.6)') 'Backtracking: f1 = ', f1, ', f0 = ', f0, ', f0+qt = ', f0+qt
      call log%write(msg, .true.)
      if (f1 <= f0 + qt) then
        break_flag = .true.
        write(msg, '(a,f20.6)') 'Backtracking: Accept step length: ', step_len
        call log%write(msg, .true.)
      else
        step_len = step_len * fpar%update%MAX_SHRINK
        break_flag = .false.
        write(msg, '(a,f20.6)') 'Backtracking: Try to shrink step length: ', step_len
        call log%write(msg, .true.)
      endif
    endif
    call synchronize_all()
    call bcast_all_singlecr(step_len)
    call bcast_all_singlel(break_flag)

  end subroutine backtracking
end module line_search