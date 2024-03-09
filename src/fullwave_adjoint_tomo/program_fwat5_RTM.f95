!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Mijian Xu
!                              Nanyang Technological University
!                              Singapore
!                           (c) Martch 2020
!
!=====================================================================
!

program fwat4_cmt3d_fwd
  use fullwave_adjoint_tomo_par

  character(len=MAX_STRING_LEN)                   :: evtsetb,evtsete 
  character(len=MAX_STRING_LEN)                   :: is_smooth 
  integer myrank
  real :: t1, t2
! MPI initialization
  call init_mpi()
  call world_rank(myrank)  
  if (myrank == 0) print *,"Running XFWAT5_RTM"
  call synchronize_all()
  call cpu_time(t1)

  ! parse command line arguments
  if (command_argument_count() /= 3) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xfwat5_RTM setb sete is_smooth'
      stop 'Please check command line arguments'
    endif
  endif
  call synchronize_all()
  call get_command_argument(1, evtsetb)     
  call get_command_argument(2, evtsete)     
  call get_command_argument(3, is_smooth)     
 
  call run_fwat5_RTM(evtsetb,evtsete,is_smooth)
  call cpu_time(t2)

  if(myrank==0) print *,'Computation time with CPU:',t2-t1


! MPI finish
  call finalize_mpi()
end program