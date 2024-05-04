!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
!

 program fwat1_fwd_measure_adj
  use fullwave_adjoint_tomo_par
  use fwat_input

  implicit none

  
  character(len=MAX_STRING_LEN)                   :: model 
  character(len=MAX_STRING_LEN)                   :: evtset 
  character(len=MAX_STRING_LEN)                   :: simu_type 
  character(len=MAX_STRING_LEN)                   :: run_opt 
  integer :: myrank, run_opt_num
 

! MPI initialization
  call init_mpi()
  call world_rank(myrank)
  if (command_argument_count() /= 4) then
     if (myrank == 0) then
       print *,'USAGE:  mpirun -np NPROC bin/xfwat1_fwd_measure_adj model set simu_type run_opt'
       print *,'model     --- name of current model, such as M00, M01, ...'
       print *,'set       --- name of event set, such as set1, set2, ...'
       print *,'simu_type --- simulation type: noise, tele, rf, telecd'
       print *,'run_opt   --- run_opt= 1 (fwd), 2 (fwd_meas), 3 (fwd_meas_adj)'
       stop 'Please check command line arguments'
     endif
   endif

  call get_command_argument(1, model)     
  call get_command_argument(2, evtset)     
  call get_command_argument(3, simu_type)     
  call get_command_argument(4, run_opt)    
  read(run_opt,*) run_opt_num 

  call read_parameter_file(myrank, .true.)
  call read_fwat_par_file()

! run the main program
  call run_fwat1_fwd_measure_adj(model,evtset,simu_type,run_opt_num)

  ! MPI finish
  call finalize_mpi()

end program fwat1_fwd_measure_adj

