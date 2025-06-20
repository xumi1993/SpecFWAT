module line_search
  use fwat_mpi
  use config
  use input_params, fpar => fwat_par_global
  use opt_io
  use preproc_fwd
  use window_chi, only : read_model_misfit
  use common_lib, only : get_dat_type

  implicit none

contains
  subroutine forward_for_simu_type(total_misfit, misfit_start, misfit_prev)
    type(PrepareFWD) :: ffwd
    logical :: is_output_backup
    real(kind=dp), intent(out) :: total_misfit, misfit_start, misfit_prev
    real(kind=dp) :: misfit_loc
    integer :: ievt

    is_output_backup = IS_OUTPUT_PREPROC
    IS_OUTPUT_PREPROC = .false.
    total_misfit = 0.0_dp
    misfit_start = 0.0_dp
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

      call read_model_misfit(model_start, ievt, misfit_loc)
      misfit_start = misfit_start + misfit_loc

      if (model_prev /= 'none') then
        call read_model_misfit(model_current, ievt, misfit_loc)
        misfit_prev = misfit_prev + misfit_loc
      else
        misfit_prev = misfit_start
      endif

      call synchronize_all()
    enddo

    call ffwd%destroy()

    call fpar%acqui%finalize()

    IS_OUTPUT_PREPROC = is_output_backup
    call synchronize_all()
  end subroutine forward_for_simu_type

end module line_search