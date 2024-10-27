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

 program fwat3_linesearch
  use fullwave_adjoint_tomo_par
  use fwat_input
  use postproc_sub

  implicit none

  character(len=MAX_STRING_LEN)                   :: evtset 
  character(len=MAX_STRING_LEN)                   :: simu_type 
 
! MPI initialization
  call init_mpi()
  call world_rank(myrank)
  if (command_argument_count() /= 1) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xfwat3_linesearch model_step'
      print *,'model_step --- name of trial models, such as M00 M01, ...'
      stop 'Please check command line arguments'
    endif
  endif

  call get_command_argument(1, model)     

  call read_parameter_file(myrank, .true.)
  call read_fwat_par_file()
  call read_database()

  call run_optim()

  ! MPI finish
  call finalize_mpi()

end program fwat3_linesearch

