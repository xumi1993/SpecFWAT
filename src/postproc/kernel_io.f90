module kernel_io
  use fwat_constants
  use input_params, only: fpar => fwat_par_global
  use fwat_mpi
  use config
  use specfem_par
  use external_model
  use projection_on_FD_grid_fwat

  implicit none
  integer :: ier
  character(len=MAX_STRING_LEN) :: fprname

contains

  subroutine read_mesh_databases_minimum(is_read_database)
    logical, intent(in) :: is_read_database
    character(len=MAX_STRING_LEN) :: database_name

    ! call read_parameter_file(.true.)

    if (is_read_database) then
      call initialize_simulation()

      ! reads in external mesh
      call read_mesh_databases()
    endif
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

    allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
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

    allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    read(IIN) data
    close(IIN)

  end subroutine read_kernel

  subroutine write_kernel(dataname, data, is_simu_type)
    character(len=*), intent(in) :: dataname
    logical, optional, intent(in) :: is_simu_type
    real(kind=cr), dimension(:,:,:,:), intent(in) :: data
    character(len=MAX_STRING_LEN) :: path

    path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)
    if (present(is_simu_type) .and. is_simu_type) then
      path = path//'_'//trim(simu_type)
    endif

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