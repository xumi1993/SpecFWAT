module model_grid_data
  use config
  use fwat_constants
  use fwat_mpi
  use input_params, only: fpar => fwat_par_global
  use specfem_par
  use projection_on_FD_grid_fwat
  use hdf5_interface
  use utils
  use external_model
  use common_lib, only: get_kernel_names

  implicit none

contains

  subroutine create_grid()
    integer :: i
    character(len=MAX_STRING_LEN) :: msg

    hx_fd_proj = fpar%grid%regular_grid_interval(1)
    hy_fd_proj = fpar%grid%regular_grid_interval(2)
    hz_fd_proj = fpar%grid%regular_grid_interval(3)
    ox_fd_proj = fpar%grid%regular_grid_min_coord(1)
    oy_fd_proj = fpar%grid%regular_grid_min_coord(2)
    oz_fd_proj = fpar%grid%regular_grid_min_coord(3)
    nx_fd_proj = fpar%grid%regular_grid_size(1)
    ny_fd_proj = fpar%grid%regular_grid_size(2)
    nz_fd_proj = fpar%grid%regular_grid_size(3)

    allocate(MEXT_V%x(nx_fd_proj), MEXT_V%y(ny_fd_proj), MEXT_V%z(nz_fd_proj))
    MEXT_V%x = [(ox_fd_proj + (i-1)*hx_fd_proj, i=1,nx_fd_proj)]
    MEXT_V%y = [(oy_fd_proj + (i-1)*hy_fd_proj, i=1,ny_fd_proj)]
    MEXT_V%z = [(oz_fd_proj + (i-1)*hz_fd_proj, i=1,nz_fd_proj)]
    MEXT_V%nx = nx_fd_proj
    MEXT_V%ny = ny_fd_proj
    MEXT_V%nz = nz_fd_proj

    if (worldrank == 0) then
      if (MEXT_V%x(1) > x_min_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Min x: ',MEXT_V%x(1), &
              ' value of grid larger than x_min_glob:', x_min_glob
        call exit_MPI(0, msg)
      endif
      if (MEXT_V%x(MEXT_V%nx) < x_max_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Max x: ',MEXT_V%x(MEXT_V%nx),&
             ' value of grid smaller than x_max_glob:', x_max_glob
        call exit_MPI(0, msg)
      endif
      if (MEXT_V%y(1) > y_min_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Min y: ',MEXT_V%y(1), &
              ' value of grid larger than y_min_glob:', y_min_glob
        call exit_MPI(0, msg)
      endif
      if (MEXT_V%y(MEXT_V%ny) < y_max_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Max y: ',MEXT_V%y(MEXT_V%ny), &
              ' value of grid smaller than y_max_glob:', y_max_glob
        call exit_MPI(0, msg)
      endif
      if (MEXT_V%z(1) > z_min_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Min z: ',MEXT_V%z(1), &
              ' value of grid larger than z_min_glob:', z_min_glob
        call exit_MPI(0, msg)
      endif
      if (MEXT_V%z(MEXT_V%nz) < z_max_glob) then
        write(msg, '(A,F20.6,A,F20.6)') 'Error: Max z: ',MEXT_V%z(MEXT_V%nz), &
              ' value of grid smaller than z_max_glob:', z_max_glob
        call exit_MPI(0, msg)
      endif
    endif

    call synchronize_all()
  end subroutine create_grid

  subroutine database2model(model_data)
    integer :: igll, jgll, kgll, ispec, iker
    real(kind=cr), dimension(:,:,:,:,:), allocatable, intent(out) :: model_data

    if (fpar%update%model_type == 1) then
      allocate(model_data(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nkernel))
      model_data(:,:,:,:,1) = sqrt((kappastore+mustore*FOUR_THIRDS)/rhostore)
      model_data(:,:,:,:,2) = sqrt(mustore/rhostore)
      model_data(:,:,:,:,3) = rhostore
    else
      ! TODO: implement anisotropic model
      call exit_MPI(0, 'Error: Anisotropic model not implemented yet')
    endif

  end subroutine database2model

  subroutine model_gll2grid(grid_model)
    real(kind=cr), dimension(:,:,:,:,:), allocatable :: gll_data
    real(kind=cr), dimension(:,:,:), allocatable :: gm
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: grid_model
    integer :: ier, i
    
    allocate(grid_model(MEXT_V%nx, MEXT_V%ny, MEXT_V%nz, nkernel), stat=ier)
    grid_model = 0.0_cr

    call database2model(gll_data)

    do i = 1, nkernel
      call gll2grid(gll_data(:,:,:,:,i), gm, .true.)
      grid_model(:,:,:,i) = gm
    end do
  end subroutine model_gll2grid

  subroutine write_grid_model(filename, grid_model)
    real(kind=cr), dimension(:,:,:,:), intent(in) :: grid_model
    character(len=*), intent(in) :: filename
    integer :: iker
    type(hdf5_file) :: h5file
    
    if (worldrank == 0) then
      call h5file%open(filename, status='new', action='write')
      call h5file%add('/x', MEXT_V%x)
      call h5file%add('/y', MEXT_V%y)
      call h5file%add('/z', MEXT_V%z)
      do iker = 1, nkernel
        call h5file%add('/'//trim(parameter_names(iker)), transpose_3(grid_model(:,:,:,iker)))
      end do
      call h5file%close(finalize=.true.)
    endif

  end subroutine write_grid_model

  subroutine write_grid_kernel_smooth(grid_kernel, filename)
    real(kind=cr), dimension(:,:,:,:), intent(in) :: grid_kernel
    character(len=*), intent(in) :: filename
    integer :: iker
    type(hdf5_file) :: h5file
    
    if (worldrank == 0) then
      call h5file%open(filename, status='new', action='write')
      call h5file%add('/x', MEXT_V%x)
      call h5file%add('/y', MEXT_V%y)
      call h5file%add('/z', MEXT_V%z)
      do iker = 1, nkernel
        call h5file%add('/'//trim(kernel_names(iker))//'_kernel_smooth',transpose_3(grid_kernel(:,:,:,iker)))
      end do
      call h5file%close(finalize=.true.)
    endif
  end subroutine write_grid_kernel_smooth

  subroutine read_grid_kernel_smooth(filename, grid_model)
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: grid_model
    character(len=*), intent(in) :: filename
    real(kind=cr), dimension(:,:,:), allocatable :: gm  
    integer :: iker

    type(hdf5_file) :: h5file
    if (worldrank == 0) then
      call h5file%open(filename, status='old', action='read')
      do iker = 1, nkernel
        call h5file%get('/'//trim(kernel_names(iker)), gm)
        grid_model(:,:,:,iker) = transpose_3(gm)
      end do
      call h5file%close(finalize=.true.)
    endif
  end subroutine read_grid_kernel_smooth

  subroutine read_grid_model(filename, grid_model)
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: grid_model
    character(len=*), intent(in) :: filename
    real(kind=cr), dimension(:,:,:), allocatable :: gm  
    integer :: iker

    type(hdf5_file) :: h5file
    if (worldrank == 0) then
      call h5file%open(filename, status='old', action='read')
      do iker = 1, nkernel
        call h5file%get('/'//trim(parameter_names(iker)), gm)
        grid_model(:,:,:,iker) = transpose_3(gm)
      end do
      call h5file%close(finalize=.true.)
    endif
  end subroutine read_grid_model
  
  subroutine gll2grid(data, grid_data, is_replace_zero)
    real(kind=cr), dimension(:,:,:,:), intent(in) :: data
    real(kind=cr), dimension(:,:,:,:), allocatable :: dat
    real(kind=cr), dimension(:,:,:), allocatable, intent(out) :: grid_data
    logical, intent(in) :: is_replace_zero
    integer :: i, j, k, ier
    type(profd)  :: projection_fd

    dat = data

    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
    call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
 
    call compute_interpolation_coeff_FD_SEM(projection_fd, worldrank)
    call Project_model_SEM2FD_grid(dat, grid_data, projection_fd, worldrank)
    
    if (is_replace_zero .and. worldrank == 0) then
      do i=1,projection_fd%nx
        do j=1,projection_fd%ny
          do k=1,projection_fd%nz
            if (grid_data(i,j,k) == 0) then
              if (MEXT_V%z(k)>=0) then
                call find_nearestZ_nonzero(grid_data,i,j,k,&
                      projection_fd%nx,projection_fd%ny,projection_fd%nz)
              else
                call find_nearestXY_nonzero(grid_data,i,j,k,&
                      projection_fd%nx,projection_fd%ny,projection_fd%nz)
              endif
            endif
          enddo
        enddo
      enddo
    endif
    call synchronize_all()

  end subroutine gll2grid



end module model_grid_data