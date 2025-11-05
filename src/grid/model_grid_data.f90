module model_grid_data
  use config
  use fwat_constants
  use fwat_mpi
  use input_params, only: fpar => fwat_par_global
  use specfem_par
  use hdf5_interface
  use utils
  use external_model
  use common_lib, only: get_kernel_names
  use projection_on_FD_grid_fwat


  implicit none

contains

  subroutine create_grid()
    integer :: i
    real(kind=cr) :: hx, hy, hz, ox, oy, oz
    integer :: nx, ny, nz

    hx = fpar%grid%regular_grid_interval(1)
    hy = fpar%grid%regular_grid_interval(2)
    hz = fpar%grid%regular_grid_interval(3)
    ox = fpar%grid%regular_grid_min_coord(1)
    oy = fpar%grid%regular_grid_min_coord(2)
    oz = fpar%grid%regular_grid_min_coord(3)
    nx = fpar%grid%regular_grid_size(1)
    ny = fpar%grid%regular_grid_size(2)
    nz = fpar%grid%regular_grid_size(3)

    allocate(ext_grid%x(nx), ext_grid%y(ny), ext_grid%z(nz))
    ext_grid%x = [(ox + (i-1)*hx, i=1,nx)]
    ext_grid%y = [(oy + (i-1)*hy, i=1,ny)]
    ext_grid%z = [(oz + (i-1)*hz, i=1,nz)]
    ext_grid%nx = nx
    ext_grid%ny = ny
    ext_grid%nz = nz
    ext_grid%dx = hx
    ext_grid%dy = hy
    ext_grid%dz = hz

    hx_fd_proj = hx
    hy_fd_proj = hy
    hz_fd_proj = hz
    ox_fd_proj = ox
    oy_fd_proj = oy
    oz_fd_proj = oz
    nx_fd_proj = nx
    ny_fd_proj = ny
    nz_fd_proj = nz

    call synchronize_all()
  end subroutine create_grid

  subroutine database2model(model_data)
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

  ! subroutine model_gll2grid(grid_model)
  !   real(kind=cr), dimension(:,:,:,:,:), allocatable :: gll_data
  !   real(kind=cr), dimension(:,:,:), allocatable :: gm
  !   real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: grid_model
  !   integer :: ier, i
    
  !   allocate(grid_model(ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel), stat=ier)
  !   grid_model = 0.0_cr

  !   call database2model(gll_data)

  !   do i = 1, nkernel
  !     call gll2grid(gll_data(:,:,:,:,i), gm, .true.)
  !     grid_model(:,:,:,i) = gm
  !   end do
  ! end subroutine model_gll2grid

  subroutine write_grid_model(filename, grid_model)
    real(kind=cr), dimension(:,:,:,:), intent(in) :: grid_model
    character(len=*), intent(in) :: filename
    integer :: iker
    type(hdf5_file) :: h5file
    
    if (worldrank == 0) then
      call h5file%open(filename, status='new', action='write')
      call h5file%add('/x', ext_grid%x)
      call h5file%add('/y', ext_grid%y)
      call h5file%add('/z', ext_grid%z)
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
      call h5file%add('/x', ext_grid%x)
      call h5file%add('/y', ext_grid%y)
      call h5file%add('/z', ext_grid%z)
      do iker = 1, nkernel
        call h5file%add('/'//trim(kernel_names(iker))//'_kernel_smooth',transpose_3(grid_kernel(:,:,:,iker)))
      end do
      call h5file%close(finalize=.true.)
    endif
  end subroutine write_grid_kernel_smooth

  subroutine write_grid(filename, key_name, grid_data)
    real(kind=cr), dimension(:,:,:), intent(in) :: grid_data
    character(len=*), intent(in) :: filename, key_name
    type(hdf5_file) :: h5file
    
    if (worldrank == 0) then
      call h5file%open(filename, status='new', action='write')
      call h5file%add('/x', ext_grid%x)
      call h5file%add('/y', ext_grid%y)
      call h5file%add('/z', ext_grid%z)
      call h5file%add('/'//trim(key_name), transpose_3(grid_data(:,:,:)))
      call h5file%close(finalize=.true.)
    endif

  end subroutine write_grid

  subroutine read_grid(filename, key_name, grid_data)
    real(kind=cr), dimension(:,:,:), allocatable, intent(out) :: grid_data
    character(len=*), intent(in) :: filename, key_name
    real(kind=cr), dimension(:,:,:), allocatable :: gm  
    type(hdf5_file) :: h5file

    if (worldrank == 0) then
      call h5file%open(filename, status='old', action='read')
      call h5file%get('/'//trim(key_name), gm)
      grid_data = transpose_3(gm)
      call h5file%close(finalize=.true.)
    endif
  end subroutine read_grid

  subroutine read_grid_kernel_smooth(filename, grid_model)
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: grid_model
    character(len=*), intent(in) :: filename
    real(kind=cr), dimension(:,:,:), allocatable :: gm  
    integer :: iker

    type(hdf5_file) :: h5file
    if (worldrank == 0) then
      allocate(grid_model(ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel))
      call h5file%open(filename, status='old', action='read')
      do iker = 1, nkernel
        call h5file%get('/'//trim(kernel_names(iker))//'_kernel_smooth', gm)
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
      allocate(grid_model(ext_grid%nx, ext_grid%ny, ext_grid%nz, nkernel))
      call h5file%open(filename, status='old', action='read')
      do iker = 1, nkernel
        call h5file%get('/'//trim(parameter_names(iker)), gm)
        grid_model(:,:,:,iker) = transpose_3(gm)
      end do
      call h5file%close(finalize=.true.)
    endif
  end subroutine read_grid_model

  subroutine gll2grid(gll_data, grid_data)
    use shared_parameters

    real(kind=cr), dimension(:,:,:,:), allocatable, intent(in) :: gll_data
    real(kind=cr), dimension(:,:,:), allocatable, intent(out) :: grid_data
    type(profd)  :: projection_fd

    ! Get regular grid properties, and precomputes all interpolation coefficients
    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
    call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

    call compute_interpolation_coeff_FD_SEM(projection_fd, worldrank)

    ! put data from SEM mesh to a regular grid
    call Project_model_SEM2FD_grid(gll_data, grid_data, projection_fd, worldrank)
  
    call synchronize_all()
  
  end subroutine gll2grid

  subroutine grid2gll(grid_data, gll_data)
    use shared_parameters

    real(kind=cr), dimension(:,:,:), allocatable, intent(in) :: grid_data
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: gll_data

    gll_data = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_AB)
    call Project_model_FD_grid2SEM(gll_data, grid_data, worldrank)

    call synchronize_all()
  end subroutine grid2gll

end module model_grid_data