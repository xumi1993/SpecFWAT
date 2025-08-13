module kernel_io
  use fwat_constants
  use input_params, only: fpar => fwat_par_global
  use fwat_mpi
  use config
  use external_model
  use projection_on_FD_grid_fwat

  implicit none
  integer :: ier
  character(len=MAX_STRING_LEN), private :: fprname, msg

contains

  subroutine read_mesh_databases_for_init()
    use specfem_par

    call initialize_simulation_fwat()

    ! reads in external mesh
    call read_mesh_databases_fwat()

    NSPEC_FWAT = NSPEC_AB
    NGLOB_FWAT = NGLOB_AB
    call synchronize_all()

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

    if (worldrank == 0) then
      if (ext_grid%x(1) > x_min_glob) then
        write(msg, '(A,F0.6,A,F20.6)') 'Error: Min x: ',ext_grid%x(1), &
              ' value of grid larger than x_min_glob:', x_min_glob
        call exit_MPI(0, trim(msg))
      endif
      if (ext_grid%x(ext_grid%nx) < x_max_glob) then
        write(msg, '(A,F0.6,A,F20.6)') 'Error: Max x: ',ext_grid%x(ext_grid%nx),&
             ' value of grid smaller than x_max_glob:', x_max_glob
        call exit_MPI(0, trim(msg))
      endif
      if (ext_grid%y(1) > y_min_glob) then
        write(msg, '(A,F0.6,A,F20.6)') 'Error: Min y: ',ext_grid%y(1), &
              ' value of grid larger than y_min_glob:', y_min_glob
        call exit_MPI(0, trim(msg))
      endif
      if (ext_grid%y(ext_grid%ny) < y_max_glob) then
        write(msg, '(A,F0.6,A,F20.6)') 'Error: Max y: ',ext_grid%y(ext_grid%ny), &
              ' value of grid smaller than y_max_glob:', y_max_glob
        call exit_MPI(0, trim(msg))
      endif
      if (ext_grid%z(1) > z_min_glob) then
        write(msg, '(A,F0.6,A,F20.6)') 'Error: Min z: ',ext_grid%z(1), &
              ' value of grid larger than z_min_glob:', z_min_glob
        call exit_MPI(0, trim(msg))
      endif
      if (ext_grid%z(ext_grid%nz) < z_max_glob) then
        write(msg, '(A,F0.6,A,F20.6)') 'Error: Max z: ',ext_grid%z(ext_grid%nz), &
              ' value of grid smaller than z_max_glob:', z_max_glob
        call exit_MPI(0, trim(msg))
      endif
    endif

    call synchronize_all()

  end subroutine read_mesh_databases_for_init

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

  subroutine kernel_cijkl2hti(ievt, ker)
    use specfem_par_elastic, only: c11store, c12store, c13store, c22store, c23store, c24store, &
                                   c33store, c44store, c45store, c55store, c66store, &
                                   rho_vp, rho_vs
    use specfem_par, only: rhostore
    integer, intent(in) :: ievt
    real(kind=cr), dimension(:,:,:,:,:), allocatable, intent(out) :: ker
    real(kind=cr), dimension(:,:,:,:), allocatable :: A, F, C, L, N, vp, vs, gc, gs, gcp, gsp
    real(kind=cr), dimension(:,:,:,:), allocatable :: c11_kl, c44_kl, c55_kl, c45_kl, &
                                                      c12_kl, c13_kl, c22_kl, c23_kl, &
                                                      c33_kl, c66_kl, rho_kl

    ! define model parameters
    L = 0.5 * (c44store + c55store)
    gc = 0.5 * (c44store - c55store)
    gs = c45store
    gcp = gc/L
    gsp = gs/L
    vp = rho_vp / rhostore
    vs = rho_vs / rhostore

    deallocate(L, gc, gs)

    ! read event kernel
    call read_event_kernel(ievt, 'rho_kernel', rho_kl)
    call read_event_kernel(ievt, 'c11_kernel', c11_kl)
    call read_event_kernel(ievt, 'c12_kernel', c12_kl)
    call read_event_kernel(ievt, 'c13_kernel', c13_kl)
    call read_event_kernel(ievt, 'c22_kernel', c22_kl)
    call read_event_kernel(ievt, 'c23_kernel', c23_kl)
    call read_event_kernel(ievt, 'c33_kernel', c33_kl)
    call read_event_kernel(ievt, 'c44_kernel', c44_kl)
    call read_event_kernel(ievt, 'c45_kernel', c45_kl)
    call read_event_kernel(ievt, 'c55_kernel', c55_kl)
    call read_event_kernel(ievt, 'c66_kernel', c66_kl)
    allocate(ker(NGLLX,NGLLY,NGLLZ,NSPEC_FWAT,nkernel),stat=ier)
    ! alpha kernel
    ker(:,:,:,:,1) = c11_kl * (2*rho_vp) + c12_kl * (2*rho_vp)+ c13_kl * (2*rho_vp) &
                   + c22_kl * (2*rho_vp) + c23_kl * (2*rho_vp)+ c33_kl * (2*rho_vp)
    ! beta kernel
    ker(:,:,:,:,2) = c12_kl * (-4*rho_vs) + c13_kl * (-4*rho_vs) + c23_kl * (-4*rho_vs) &
                   + c44_kl * (2*gcp*rho_vs + 2*rho_vs) + c45_kl *(2*gsp*rho_vs) &
                   + c55_kl * (-2*gcp*rho_vs + 2*rho_vs) + c66_kl *(2*rho_vs)
    ! rho kernel
    ker(:,:,:,:,3) = rho_kl + c11_kl * (vp**2) + c12_kl * (vp**2 - 2*vs**2) &
                   + c13_kl * (vp**2 - 2*vs**2) + c22_kl * (vp**2) + c23_kl *(vp**2 - 2*vs**2) &
                   + c33_kl * (vp**2) + c44_kl * (gcp*vs**2 + vs**2) + c45_kl * (gsp*vs**2) &
                   + c55_kl * (-gcp*vs**2 + vs**2) + c66_kl *(vs**2)
    ! gcp kernel
    ker(:,:,:,:,4) = c44_kl * (rho_vs*vs) + c55_kl * (-rho_vs*vs)
    ! gsp kernel
    ker(:,:,:,:,5) = c45_kl * (rho_vs*vs) ! ker_Gs

    ! relative scaling
    ker(:,:,:,:,1) = ker(:,:,:,:,1) * vp
    ker(:,:,:,:,2) = ker(:,:,:,:,2) * vs
    ker(:,:,:,:,3) = ker(:,:,:,:,3) * rhostore
  end subroutine kernel_cijkl2hti

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
    character(len=*), intent(in) :: dataname, kernel_path
    character(len=MAX_STRING_LEN) :: path
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: data

    path = trim(kernel_path)
    call create_name_database(fprname, worldrank, path)
    open(unit=IIN, file=trim(fprname)//trim(dataname)//'.bin', status='old', action='read', form='unformatted', iostat=ier)
    if (ier /= 0) then
      write(0, *) 'Error could not open database file: ',trim(fprname)//trim(dataname)//'.bin'
      call exit_mpi(myrank,'Error opening database file')
    endif

    allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_FWAT),stat=ier)
    read(IIN) data
    close(IIN)

  end subroutine read_kernel

  subroutine write_kernel(kernel_path, dataname, data)
    character(len=*), intent(in) :: dataname, kernel_path
    real(kind=cr), dimension(:,:,:,:), intent(in) :: data
    character(len=MAX_STRING_LEN) :: path

    path = trim(kernel_path)
    call create_name_database(fprname, worldrank, path)
    open(unit=IOUT, file=trim(fprname)//trim(dataname)//'.bin', status='replace', action='write', form='unformatted', iostat=ier)
    if (ier /= 0) then
      write(0, *) 'Error could not open database file: ',trim(fprname)//trim(dataname)//'.bin'
      call exit_mpi(myrank,'Error opening database file')
    endif

    write(IOUT) data
    close(IOUT)

  end subroutine write_kernel

end module kernel_io