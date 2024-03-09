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
subroutine fullwave_adjoint_tomo_main()

  use fullwave_adjoint_tomo_par
  use fwat_input
  use postproc_sub, only: select_set_range, type_name, isetb, isete
  
  implicit none

  integer, parameter :: max_arg = 2
  character(len=MAX_STRING_LEN)                   :: args(max_arg)  
  character(len=MAX_STRING_LEN)                   :: model
  character(len=MAX_STRING_LEN)                   :: evtset
  integer :: myrank, maxit, iter_start, iter_end, iter, i, j

  ! initialize MPI
  call world_rank(myrank)
  !!!##############################################################################################################################
  !!! ---------------------------------------------- INITIALIZE RUNTIME ----------------------------------------------------------
  !!!##############################################################################################################################
  if (command_argument_count() /= max_arg) then
     if (myrank == 0) then
       print *,'USAGE:  mpirun -np NPROC bin/xfullwave_adjoint_tomo iter_start iter_end set simu_type'
       print *,'iter_start   --- starting iteration number'
       print *,'iter_end     --- ending iteration number'
       stop 'Please check command line arguments'
     endif
   endif

  ! read command line arguments
  do i = 1, max_arg
    call get_command_argument(i, args(i))
  end do
  read(args(1),*) iter_start
  read(args(2),*) iter_end
  call synchronize_all()
 
  !!!##############################################################################################################################
  !!! -------------------------------------------  initialize starting model -----------------------------------------------------
  !!!##############################################################################################################################
  ! read parameter file
  call read_parameter_file(myrank,.true.)
  call read_fwat_par_file()
  is_read_database = .false.
  ! read mesh parameter file
  call read_mesh_parameter_file()
  call synchronize_all()


  !!!##############################################################################################################################
  !!! -------------------------------  different running mode : forward or FWI ----------------------------------------------------
  !!!##############################################################################################################################
  do iter = iter_start, iter_end
    write(model,'("M",I2.2)') iter
   ! initialize starting model mesh
    call meshfem3d_fwat()
    ! generate database for forward simulation
    call generate_database_fwat()
    ! run forward and adjoint simulation
    do i = 1, 3
      if (tomo_par%INV_TYPE(i)) then
        type_name = tomo_par%INV_TYPE_NAME(i)    
        call select_set_range()
        do j = isetb, isete
          write(evtset,'("set",I0)') j
          call run_fwat1_fwd_measure_adj(model,evtset,type_name,3)
        enddo
      endif
    enddo
    ! run postprocessing
    call run_fwat2_postproc_opt(model)
  end do
end subroutine fullwave_adjoint_tomo_main

