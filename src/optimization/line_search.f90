module line_search
  use fwat_mpi
  use config
  use input_params, fpar => fwat_par_global
  use opt_io
  use preproc_fwd
  use window_chi, only : WindowChi

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
    if (simu_type == SIMU_TYPE_NOISE) then
      dat_type = 'noise'
    elseif (simu_type == SIMU_TYPE_TELE) then
      if (fpar%postproc%TELE_TYPE == 1) then
        dat_type = 'tele'
      elseif (fpar%postproc%TELE_TYPE == 2) then
        dat_type = 'rf'
      else
        call exit_MPI(0, 'LINE SEARCH: Unknown teleseismic type')
      endif
    else
      call exit_MPI(0, 'LINE SEARCH: Unknown simulation type')
    endif

    call fpar%select_simu_type()

    call fpar%acqui%read()

    ! initialize fwd
    call ffwd%init(.false.)

    do ievt = 1, fpar%acqui%nevents
      ffwd%ievt = ievt

      ! prepare simulation
      call ffwd%prepare_for_event()

      ! prepare fk wavefield
      call ffwd%calc_or_read_fk_wavefield()

      ! run forward simulation
      call ffwd%simulation()

      total_misfit = total_misfit + ffwd%obj_func

      call read_model_misfit(model_start, ievt, misfit_loc)
      misfit_start = misfit_start + misfit_loc

      if (model_prev /= 'none') then
        call read_model_misfit(model_prev, ievt, misfit_loc)
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

  subroutine read_model_misfit(model_name_in, ievt, misfit)
    use common_lib, only: get_band_name
    character(len=MAX_STRING_LEN), intent(in) :: model_name_in
    integer, intent(in) :: ievt
    real(kind=dp), intent(out) :: misfit
    character(len=MAX_STRING_LEN) :: band_name
    integer :: iflt
    type(WindowChi) :: wchi

    misfit = 0.0_dp
    do iflt = 1, fpar%sim%NUM_FILTER
      if (dat_type /= 'rf') then
        call get_band_name(fpar%sim%SHORT_P(iflt), fpar%sim%LONG_P(iflt), band_name)
      else
        write(band_name, '("F",F3.1)') fpar%sim%rf%f0(iflt)
      endif

      call wchi%read(model_name_in, ievt, band_name)

      misfit = misfit + wchi%sum_chi(29)

      call wchi%finalize()
    enddo

  end subroutine read_model_misfit

end module line_search