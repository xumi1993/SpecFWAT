program fwat_fwd_measure_adj
use config
use fwat_mpi
use common_lib, only: get_simu_type
use imput_params, fpar => fwat_par_global
use obs_data, fdat => fwat_evt_data_global
use preproc_fwd
use specfem_par, only: DT, NSTEP

implicit none
integer :: nargs
integer, parameter :: max_num_args = 3
type(PrepareFWD) :: ffwd
character(len=MAX_STRING_LEN) :: usage, evt_index_str

call init_mpi()
call init_mpi_fwat()

usage = 'Usage: fwat_fwd_measure_adj model data_type [event_index]'
! read command line arguments
nargs = command_argument_count()
if (nargs > max_num_args) then
  if (worldrank == 0) then
    print *, usage
    call exit_MPI(0, 'ERROR: Too more arguments')
  endif
else if (nargs < max_num_args - 1 ) then
  if (worldrank == 0) then
    print *, usage
    call exit_MPI(0, 'ERROR: Too few arguments')
  endif
else
  call get_command_argument(1, model_name)
  call get_command_argument(2, dat_type)
endif
if (nargs == max_num_args) then
  call get_command_argument(3, evt_index_str)
  read(evt_index_str, *) ffwd%ievt
  single_run = .true.
endif
call get_simu_type()

! read input parameters
call fpar%read(FWAT_PAR_FILE)

! select simu_type
call fpar%select_simu_type()

! read src_rec for this data type
call fpar%acqui%read()

! read observed data
! call fdat%read_stations(fpar%acqui%evtid_names(1))

! initialize fwd
call ffwd%init(FORWARD_ADJOINT)

call ffwd%calc_fk_wavefield()

call ffwd%prepare_for_event()

call finalize_mpi()

end program fwat_fwd_measure_adj