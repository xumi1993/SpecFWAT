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

  character(len=MAX_STRING_LEN)                   :: model, simu_type
  integer myrank, nargs
  real :: t1, t2
! MPI initialization
  call init_mpi()
  call world_rank(myrank)  
  if (myrank == 0) print *,"Running XFWAT2_POSTPROC_OPT"
  call synchronize_all()
  call cpu_time(t1)

  nargs = command_argument_count()
  ! parse command line arguments
  if (nargs < 1 .or. nargs > 2) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xfwat2_postproc_opt model [simu_type]'
      stop 'Please check command line arguments'
    endif
  endif
  call synchronize_all()
  call get_command_argument(1, model)
  if (nargs == 2) then
    call get_command_argument(2, simu_type)
  endif
  call read_parameter_file(myrank, .true.)
  call read_fwat_par_file()
 
  if (nargs == 1) then
    call run_fwat2_postproc_opt(model)
  else
    call run_fwat2_postproc_single(model, simu_type)
  endif
  call cpu_time(t2)

  if(myrank==0) print *,'Computation time with CPU:',t2-t1


! MPI finish
  call finalize_mpi()

end program fwat2_postproc_opt

