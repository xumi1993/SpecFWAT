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
  use fwat_utils, only: get_mesh_file_path
  
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
  call read_parameter_file(myrank, .true.)
  call prepare_inversion()
  call read_fwat_par_file()
  is_read_database = .false.
  call init_inversion()

  !!!##############################################################################################################################
  !!! -------------------------------  different running mode : forward or FWI ----------------------------------------------------
  !!!##############################################################################################################################
  do iter = iter_start, iter_end
    write(model,'("M",I2.2)') iter

    ! run forward and adjoint simulation
    do i = 1, NUM_INV_TYPE
      if (tomo_par%INV_TYPE(i)) then
        ! read mesh parameter file
        call read_mesh_parameter_file_fwat(get_mesh_file_path(i))
        ! initialize starting model mesh
        call meshfem3d_fwat()
        ! generate database for forward simulation
        call generate_database_fwat(USE_H5)
        
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

    ! line search
    ! call run_optim()
    
  end do
end subroutine fullwave_adjoint_tomo_main


subroutine prepare_inversion()
  use constants
  use shared_parameters

  call system('mkdir -p optimize')
  call system('mkdir -p '//trim(OUTPUT_FILES_BASE))
  call system('mkdir -p '//trim(LOCAL_PATH))
end subroutine prepare_inversion


subroutine init_inversion()

  use fullwave_adjoint_tomo_par
  use fwat_input
  use generate_databases_par, only: IMODEL
  use fwat_utils, only: get_mesh_file_path
  use specfem_par

  implicit none

  integer :: itype
  logical :: exist

  if (myrank == 0) then
    ! check IMODEL
    if (IMODEL /= 7 .and. IMODEL /= 6) then
      print *, 'ERROR: MODEL must be gll or external'
      stop
    endif

    ! check for joint
    if (count(tomo_par%INV_TYPE) > 1) then
      is_joint = .true.
      do itype = 1, NUM_INV_TYPE
        if (tomo_par%INV_TYPE(itype)) then
          inquire(file=get_mesh_file_path(itype), exist=exist)
          if (.not. exist) then
            print *, 'ERROR: No found mesh file of ', trim(get_mesh_file_path(itype))
            stop
          endif
        endif
      end do
    endif ! is_joint
  endif ! myrank == 0
  call bcast_all_singlel(is_joint)
  
  call read_mesh_parameter_file_fwat(get_mesh_file_path(2))
  call meshfem3d_fwat()
  call generate_database_fwat(USE_H5)
  block
    use meshfem3D_par, only: NEX_XI, NEX_ETA, NER
    use generate_databases_par, only: xstore, ystore, zstore
    real(kind=CUSTOM_REAL) :: x_min_loc, x_max_loc, y_min_loc, y_max_loc, z_min_loc, z_max_loc
    real(kind=CUSTOM_REAL) :: x_min, x_max, y_min, y_max, z_min, z_max, deltax, deltay, deltaz

    x_min_loc = minval(xstore)
    x_max_loc = maxval(xstore)
    y_min_loc = minval(ystore)
    y_max_loc = maxval(ystore)
    z_min_loc = minval(zstore)
    z_max_loc = maxval(zstore)
    call min_all_all_cr(x_min_loc, x_min)
    call max_all_all_cr(x_max_loc, x_max)
    call min_all_all_cr(y_min_loc, y_min)
    call max_all_all_cr(y_max_loc, y_max)
    call min_all_all_cr(z_min_loc, z_min)
    call max_all_all_cr(z_max_loc, z_max)

    deltax = (x_max - x_min) / NEX_XI / NGLLX
    deltay = (y_max - y_min) / NEX_ETA / NGLLY
    deltaz = (z_max - z_min) / NER / NGLLZ

    call rg%init(x_min, x_max, y_min, y_max, z_min, z_max, deltax, deltay, deltaz)
  
  end block
  call synchronize_all()



end subroutine init_inversion

