subroutine run_forward_only_or_inversion(evtid, simu_opt)
  use config
  use fwat_mpi
  use common_lib, only: find_string
  use input_params, fpar => fwat_par_global
  use obs_data, fdat => fwat_evt_data_global
  use logger
  implicit none

  character(len=MAX_STRING_LEN), intent(in) :: evtid ! event name
  integer, intent(in) :: simu_opt ! 1: forward only, 2: forward and measure_adj, 3: forward and adjoint
  integer :: ievt
  type(logger_type) :: log

  ievt = find_string(fpar%acqui%evtid_names, evtid)

  call log%init(trim(fpar%acqui%out_fwd_path(ievt))//'/fwd_measure_adj.log')
  call log%write('Running forward simulation for event '//trim(evtid))

  ! create directory for forward simulation
  if (worldrank == 0) then
    call system('mkdir -p '//trim(fpar%acqui%out_fwd_path(ievt)))
    call system('mkdir -p '//trim(fpar%acqui%out_fwd_path(ievt))//'/'//trim(OUTPUT_PATH))
    if (simu_opt > 1) call system('mkdir -p '//trim(fpar%acqui%out_fwd_path(ievt))//'/'//trim(ADJOINT_PATH))
    if (simu_opt > 2) call system('mkdir -p '//trim(fpar%acqui%out_fwd_path(ievt))//'/'//trim(EKERNEL_PATH))
  endif
  call synchronize_all()


end subroutine run_forward_only_or_inversion