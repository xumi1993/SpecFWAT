module axisem_coupling
  use specfem_par
  use specfem_par_coupling
  use input_params, only: fpar => fwat_par_global
  implicit none

contains
  subroutine setup_axisem_coupling()
    ! local parameters
    integer :: ier
    character(len=MAX_STRING_LEN) :: prname_trac

      ! AxiSEM coupling
    allocate(Veloc_axisem(3,NGLLSQUARE*num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2192')
    Veloc_axisem(:,:) = 0._CUSTOM_REAL

    allocate(Tract_axisem(3,NGLLSQUARE*num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2193')
    Tract_axisem(:,:) = 0._CUSTOM_REAL

    ! user output
    ! if (myrank == 0) then
    !   write(IMAIN,*) '  tractions: opening files ', trim(prname_trac) // 'sol_axisem'
    !   write(IMAIN,*)
    !   call flush_IMAIN()
    ! endif
    ! debug
    !write(*,*) 'OPENING ', trim(prname_trac) // 'sol_axisem'
    call create_name_database(prname_trac,myrank,TRACTION_PATH)

    open(unit=IIN_veloc_dsm,file=trim(prname_trac) // 'sol_axisem',status='old', &
          action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open file ',trim(prname_trac) // 'sol_axisem'
      print *,'Please check if traction file exists for coupling with AxiSEM...'
      stop 'Error opening tractions file proc****_sol_axisem'
    endif
  end subroutine setup_axisem_coupling
end module axisem_coupling