program fwat_fwd_measure_adj
use config
use fwat_mpi
use common_lib, only: get_dat_type
use input_params, only: fpar => fwat_par_global
use preproc_fwd
use argparse, only: parse_args_fwd_meas_adj
use param_check, only: check_model_name, check_preproc_params, check_nevents

implicit none
integer :: nsim, iter
integer, parameter :: max_num_args = 4
type(PrepareFWD) :: ffwd
logical :: BROADCAST_AFTER_READ = .true.

call init_mpi()
call init_mpi_fwat()

call parse_args_fwd_meas_adj(ffwd%ievt)

! read input parameters
call fpar%read(FWAT_PAR_FILE)
call read_parameter_file(BROADCAST_AFTER_READ)
local_path_backup = trim(LOCAL_PATH)

! select simu_type
call fpar%select_simu_type()

! Reject an inconsistent setup before running any simulation.
call check_model_name(iter)
call check_preproc_params(iter)

call get_dat_type()

! read src_rec for this data type
call fpar%acqui%read()
call check_nevents()

! initialize fwd
call ffwd%init()

if (single_run) then
  nsim = ffwd%ievt
else
  nsim = fpar%acqui%nevents
  ffwd%ievt = 1
endif

do while (ffwd%ievt <= nsim)
  ! prepare simulation
  call ffwd%prepare_for_event()
  
  ! run forward simulation
  call ffwd%simulation()

  ffwd%ievt = ffwd%ievt + 1
enddo

! Free shared memory
call fpar%acqui%finalize()

call log%write('*******************************************', .false.)
call log%write('*********** PRE-PROCESSING DONE ***********', .false.)
call log%write('*******************************************', .false.)
call log%finalize()

call synchronize_all()

call finalize_mpi()

end program fwat_fwd_measure_adj