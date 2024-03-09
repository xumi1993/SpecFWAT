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

 program dynmbat_sele_ctrlgrp

  use fullwave_adjoint_tomo_par

  character(len=MAX_STRING_LEN)                   :: model 
  character(len=MAX_STRING_LEN)                   :: evtsetb,evtsete,str 
  integer myrank,nctrl
  real :: t1, t2
! MPI initialization
  call init_mpi()
  call world_rank(myrank)  
  if (myrank == 0) print *,"Running XMINIBATCH_SELECT_CONTROL_GROUP"
  call synchronize_all()
  call cpu_time(t1)

  ! parse command line arguments
  if (command_argument_count() /= 4) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xdynbat_sele_ctrlgrp model setb sete nctrl'
      stop 'Please check command line arguments'
    endif
  endif
  call synchronize_all()
  call get_command_argument(1, model)     
  call get_command_argument(2, evtsetb)     
  call get_command_argument(3, evtsete)     
  call get_command_argument(4, str) 
  read(str,*) nctrl
 
  call run_sele_ctrlgrp(model,evtsetb,evtsete,nctrl)
  call cpu_time(t2)

  if(myrank==0) print *,'Computation time with CPU:',t2-t1


! MPI finish
  call finalize_mpi()

end program dynmbat_sele_ctrlgrp

