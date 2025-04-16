module kernel_io
  use fwat_constants
  use input_params, only: fpar => fwat_par_global
  use fwat_mpi
  use config
  use external_model
  use projection_on_FD_grid_fwat

  implicit none
  integer :: ier
  character(len=MAX_STRING_LEN) :: fprname

contains

  subroutine read_mesh_databases_minimum(is_read_database)
    use specfem_par

    logical, intent(in) :: is_read_database
    character(len=MAX_STRING_LEN) :: database_name

    ! call read_parameter_file(.true.)

    if (is_read_database) then
      call initialize_simulation()

      ! reads in external mesh
      call read_mesh_databases()
    endif
    NSPEC_FWAT = NSPEC_AB
    ! call synchronize_all()

    call check_mesh_distances(worldrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                              x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                              elemsize_min_glob,elemsize_max_glob, &
                              distance_min_glob,distance_max_glob)
    call bcast_all_singlecr(x_min_glob)
    call bcast_all_singlecr(x_max_glob)
    call bcast_all_singlecr(y_min_glob)
    call bcast_all_singlecr(y_max_glob)
    call bcast_all_singlecr(z_min_glob)
    call bcast_all_singlecr(z_max_glob)
    call bcast_all_singlecr(elemsize_min_glob)
    call bcast_all_singlecr(elemsize_max_glob)
    call bcast_all_singlecr(distance_min_glob)
    call bcast_all_singlecr(distance_max_glob)

  end subroutine read_mesh_databases_minimum

  subroutine get_mesh_coord()
    use meshfem3D_subs, only: meshfem3D_fwat
    use generate_databases_subs
    real(kind=cr) :: x_min,x_max
    real(kind=cr) :: y_min,y_max
    real(kind=cr) :: z_min,z_max
    character(len=MAX_STRING_LEN) :: msg

    call meshfem3D_fwat(fpar%sim%MESH_PAR_FILE)
    call generate_databases_fwat(.false.)
    NSPEC_FWAT = NSPEC_AB
    xstore_fwat = xstore
    ystore_fwat = ystore
    zstore_fwat = zstore
  
    x_min = minval(xstore_fwat)
    x_max = maxval(xstore_fwat)

    y_min = minval(ystore_fwat)
    y_max = maxval(ystore_fwat)

    z_min = minval(zstore_fwat)
    z_max = maxval(zstore_fwat)

    ! min and max dimensions of the model
    call min_all_cr(x_min,x_min_glob)
    call max_all_cr(x_max,x_max_glob)

    call min_all_cr(y_min,y_min_glob)
    call max_all_cr(y_max,y_max_glob)

    call min_all_cr(z_min,z_min_glob)
    call max_all_cr(z_max,z_max_glob)

    if (worldrank == 0) then
      if (MEXT_V%x(1) > x_min_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Min x: ',MEXT_V%x(1), &
              ' value of grid larger than x_min_glob:', x_min_glob
        call exit_MPI(0, trim(msg))
      endif
      if (MEXT_V%x(MEXT_V%nx) < x_max_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Max x: ',MEXT_V%x(MEXT_V%nx),&
             ' value of grid smaller than x_max_glob:', x_max_glob
        call exit_MPI(0, trim(msg))
      endif
      if (MEXT_V%y(1) > y_min_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Min y: ',MEXT_V%y(1), &
              ' value of grid larger than y_min_glob:', y_min_glob
        call exit_MPI(0, trim(msg))
      endif
      if (MEXT_V%y(MEXT_V%ny) < y_max_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Max y: ',MEXT_V%y(MEXT_V%ny), &
              ' value of grid smaller than y_max_glob:', y_max_glob
        call exit_MPI(0, trim(msg))
      endif
      if (MEXT_V%z(1) > z_min_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Min z: ',MEXT_V%z(1), &
              ' value of grid larger than z_min_glob:', z_min_glob
        call exit_MPI(0, trim(msg))
      endif
      if (MEXT_V%z(MEXT_V%nz) < z_max_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Max z: ',MEXT_V%z(MEXT_V%nz), &
              ' value of grid smaller than z_max_glob:', z_max_glob
        call exit_MPI(0, trim(msg))
      endif
    endif

    call synchronize_all()
  end subroutine get_mesh_coord

  subroutine read_event_kernel(ievt, dataname, data)
    integer, intent(in) :: ievt
    character(len=*), intent(in) :: dataname
    character(len=MAX_STRING_LEN) :: path
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: data

    path = trim(fpar%acqui%out_fwd_path(ievt))//'/'//trim(EKERNEL_PATH)
    call create_name_database(fprname, worldrank, path)
    ! path = fprname(1:len_trim(fprname))//'/'//trim(dataname)
    open(unit=IIN, file=trim(fprname)//trim(dataname)//'.bin', status='old', action='read', form='unformatted', iostat=ier)
    if (ier /= 0) then
      write(0, *) 'Error could not open database file: ',trim(fprname)//trim(dataname)//'.bin'
      call exit_mpi(myrank,'Error opening database file')
    endif

    allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_FWAT),stat=ier)
    read(IIN) data
    close(IIN)

  end subroutine read_event_kernel

  subroutine remove_event_kernel(ievt, dataname)
    integer, intent(in) :: ievt
    character(len=*), intent(in) :: dataname
    character(len=MAX_STRING_LEN) :: path

    path = trim(fpar%acqui%out_fwd_path(ievt))//'/'//trim(EKERNEL_PATH)
    call create_name_database(fprname, worldrank, path)
    open(IOUT, file=trim(fprname)//trim(dataname)//'.bin', status='old', iostat=ier)
    if(ier == 0) close(IOUT, status='delete')

  end subroutine remove_event_kernel

  subroutine read_kernel(kernel_path, dataname, data)
    character(len=*), intent(in) :: dataname
    character(len=MAX_STRING_LEN) :: kernel_path
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: data

    call create_name_database(fprname, worldrank, kernel_path)
    open(unit=IIN, file=trim(fprname)//trim(dataname)//'.bin', status='old', action='read', form='unformatted', iostat=ier)
    if (ier /= 0) then
      write(0, *) 'Error could not open database file: ',trim(fprname)//trim(dataname)//'.bin'
      call exit_mpi(myrank,'Error opening database file')
    endif

    allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_FWAT),stat=ier)
    read(IIN) data
    close(IIN)

  end subroutine read_kernel

  subroutine write_kernel(path, dataname, data)
    character(len=*), intent(in) :: dataname
    real(kind=cr), dimension(:,:,:,:), intent(in) :: data
    character(len=MAX_STRING_LEN) :: path

    call create_name_database(fprname, worldrank, path)
    ! path = fprname(1:len_trim(fprname))//'/'//trim(dataname)
    open(unit=IOUT, file=trim(fprname)//trim(dataname)//'.bin', status='replace', action='write', form='unformatted', iostat=ier)
    if (ier /= 0) then
      write(0, *) 'Error could not open database file: ',trim(fprname)//trim(dataname)//'.bin'
      call exit_mpi(myrank,'Error opening database file')
    endif

    write(IOUT) data
    close(IOUT)

  end subroutine write_kernel

end module kernel_io