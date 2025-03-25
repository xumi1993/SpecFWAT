program fwat_fwd_measure_adj
use config
use fwat_mpi
use common_lib, only: get_simu_type
use input_params, fpar => fwat_par_global
use obs_data, fdat => fwat_evt_data_global
use preproc_fwd
use specfem_par, only: DT, NSTEP
use argparse, only: parse_args_fwd_meas_adj

implicit none
integer :: nargs, nsim
integer, parameter :: max_num_args = 4
type(PrepareFWD) :: ffwd
character(len=MAX_STRING_LEN) :: usage, evt_index_str, run_mode_str

call init_mpi()
call init_mpi_fwat()

call parse_args_fwd_meas_adj(ffwd%ievt)

call get_simu_type()

! read input parameters
call fpar%read(FWAT_PAR_FILE)

! select simu_type
call fpar%select_simu_type()

! read src_rec for this data type
call fpar%acqui%read()

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

  ! prepare fk wavefield
  call ffwd%calc_or_read_fk_wavefield()
  
  ! run forward simulation
  call ffwd%simulation()

  ffwd%ievt = ffwd%ievt + 1
enddo

call log%write('*******************************************', .false.)
call log%write('*********** PRE-PROCESSING DONE ***********', .false.)
call log%write('*******************************************', .false.)

call finalize_mpi()

end program fwat_fwd_measure_adj