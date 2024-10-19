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

 program fwat2_postproc_opt

  use fullwave_adjoint_tomo_par
  use fwat_input

  character(len=MAX_STRING_LEN)                   :: model 
  ! character(len=MAX_STRING_LEN)                   :: evtsetb,evtsete 
  ! character(len=MAX_STRING_LEN)                   :: is_smooth 
  integer myrank
  real :: t1, t2
! MPI initialization
  call init_mpi()
  call world_rank(myrank)  
  if (myrank == 0) print *,"Running XFWAT2_POSTPROC_OPT"
  call synchronize_all()
  call cpu_time(t1)

  ! parse command line arguments
  if (command_argument_count() /= 1) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xfwat2_postproc_opt model'
      stop 'Please check command line arguments'
    endif
  endif
  call synchronize_all()
  call get_command_argument(1, model)     
  ! call get_command_argument(2, evtsetb)     
  ! call get_command_argument(3, evtsete)     
  ! call get_command_argument(4, is_smooth)     
  call read_parameter_file(myrank, .true.)
  call read_fwat_par_file()
 
  ! call run_fwat2_postproc_opt(model,evtsetb,evtsete,is_smooth)
  call run_fwat2_postproc_opt(model)
  call run_optim_man()
  call cpu_time(t2)

  if(myrank==0) print *,'Computation time with CPU:',t2-t1


! MPI finish
  call finalize_mpi()

end program fwat2_postproc_opt

