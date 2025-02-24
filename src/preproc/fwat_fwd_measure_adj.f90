program fwat_fwd_measure_adj
use config
use fwat_mpi
use common_lib, only: get_simu_type
use imput_params, fpar => fwat_par_global
use obs_data, fdat => fwat_evt_data_global

implicit none
integer :: nargs
integer, parameter :: num_args = 3

call init_mpi()
call init_mpi_fwat()

! read command line arguments
nargs = command_argument_count()
if (nargs /= num_args) then
  if (worldrank == 0) call exit_MPI(0, 'Usage: fwat_fwd_measure_adj model set_name data_type')
else
  call get_command_argument(1, model_name)
  call get_command_argument(2, set_name)
  call get_command_argument(3, dat_type)
endif
call get_simu_type()

! read input parameters
call fpar%read(FWAT_PAR_FILE)

! select simu_type
call fpar%select_simu_type()

! read src_rec
call fpar%acqui%read()

! read observed data
call fdat%read_stations(fpar%acqui%evtid_names(1))

call finalize_mpi()

end program fwat_fwd_measure_adj