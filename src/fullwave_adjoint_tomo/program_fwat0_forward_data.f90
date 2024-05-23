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

 program fwat0_forward_data
  use fullwave_adjoint_tomo_par

  implicit none

  
  character(len=MAX_STRING_LEN)                   :: model 
  character(len=MAX_STRING_LEN)                   :: evtset 
  character(len=MAX_STRING_LEN)                   :: simu_type 
  character(len=MAX_STRING_LEN)                   :: simu_type_list( 4 )
  integer :: myrank
 
! MPI initialization
  call init_mpi()
  call world_rank(myrank)
  if (command_argument_count() /= 3) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xfwat0_forward_data model_step evt_set simu_type'
      print *,'model_step --- name of trial models, such as M00, M01, ...'
      print *,'set        --- name of event set: such as set1, set2, ...'
      print *,'simu_type  --- simulation type: noise, tele, leq, rf'
      stop 'Please check command line arguments'
    endif
  endif

  call get_command_argument(1, model)     
  call get_command_argument(2, evtset)     
  call get_command_argument(3, simu_type)     
  
  simu_type_list(1) = 'noise'
  simu_type_list(2) = 'tele'
  simu_type_list(3) = 'rf'
  simu_type_list(4) = 'leq'
  
  if (.not. any(trim(simu_type) == simu_type_list)) then
    stop 'Please check simu_type. noise, tele and rf are available'
  endif

  call read_parameter_file(myrank, .true.)

! run the main program
  call run_forward_data(model,evtset,simu_type)

  ! MPI finish
  call finalize_mpi()

end program fwat0_forward_data

