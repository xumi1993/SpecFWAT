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

  implicit none

  
  character(len=MAX_STRING_LEN)                   :: model 
  character(len=MAX_STRING_LEN)                   :: evtset 
  character(len=MAX_STRING_LEN)                   :: simu_type 
  integer :: myrank
 
! MPI initialization
  call init_mpi()
  call world_rank(myrank)
  if (command_argument_count() /= 3) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xfwat3_linesearch model_step ls simu_type'
      print *,'model_step --- name of trial models, such as M00_step0.010, M00_step0.020, ...'
      print *,'ls        --- name of event set: ls (read src_rec/sources_ls.dat)'
      print *,'simu_type  --- simulation type: noise, tele'
      stop 'Please check command line arguments'
    endif
  endif

  call get_command_argument(1, model)     
  call get_command_argument(2, evtset)     
  call get_command_argument(3, simu_type)     

  call read_parameter_file(myrank, .true.)
  call read_fwat_par_file()

! run the main program
  call run_linesearch(model,evtset,simu_type)

  ! MPI finish
  call finalize_mpi()

end program fwat3_linesearch

