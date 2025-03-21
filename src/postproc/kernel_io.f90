module kernel_io
  use fwat_constants
  use input_params, fpar => fwat_par_global
  use fwat_mpi
  use config
  use generate_databases_par

  implicit none
  integer :: ier
  character(len=MAX_STRING_LEN) :: fprname

contains

  subroutine read_mesh_databases_minimum()
    character(len=MAX_STRING_LEN) :: database_name
    integer :: NSPEC_IRREGULAR

    ! read mesh databases
    ! sets file name
    call create_name_database(fprname,worldrank,LOCAL_PATH)
    database_name = fprname(1:len_trim(fprname))//'external_mesh.bin'

    open(unit=IIN,file=trim(database_name),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error could not open database file: ',trim(database_name)
      call exit_mpi(myrank,'Error opening database file')
    endif

    read(IIN) NSPEC_AB
    read(IIN) NGLOB_AB
    read(IIN) NSPEC_IRREGULAR

    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    allocate(xstore(NGLLX,NGLLY,NGLLZ,NGLOB_AB),&
             ystore(NGLLX,NGLLY,NGLLZ,NGLOB_AB),&
             zstore(NGLLX,NGLLY,NGLLZ,NGLOB_AB),stat=ier)

    read(IIN) ibool

    read(IIN) xstore
    read(IIN) ystore
    read(IIN) zstore
    close(IIN)
    call synchronize_all()

    call check_mesh_distances(worldrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                              x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                              elemsize_min_glob,elemsize_max_glob, &
                              distance_min_glob,distance_max_glob)

  end subroutine read_mesh_databases_minimum

  subroutine read_kernel(ievt, dataname, data)
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

  end subroutine read_kernel

  subroutine write_kernel(dataname, data)
    character(len=*), intent(in) :: dataname
    real(kind=cr), dimension(:,:,:,:), intent(in) :: data
    character(len=MAX_STRING_LEN) :: path

    path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)
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