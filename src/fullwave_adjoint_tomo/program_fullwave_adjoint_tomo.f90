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

 program fullwave_adjoint_tomo

! MPI initialization
  call init_mpi()

! run the main program
  call fullwave_adjoint_tomo_main()

  ! MPI finish
  call finalize_mpi()

end program fullwave_adjoint_tomo

