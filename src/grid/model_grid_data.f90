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
    character(len=MAX_STRING_LEN) :: msg
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

  subroutine write_grid_model_iso(filename, grid_model_iso)
    real(kind=cr), dimension(:,:,:,:), intent(in) :: grid_model_iso
    character(len=*), intent(in) :: filename
    integer :: iker
    type(hdf5_file) :: h5file
    
    if (worldrank == 0) then
      call h5file%open(filename, status='old', action='rw')
      do iker = 1, 3
        call h5file%add('/'//trim(MODEL_ISO(iker)), transpose_3(grid_model_iso(:,:,:,iker)))
      end do
      call h5file%close(finalize=.true.)
    endif

  end subroutine write_grid_model_iso

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
    integer :: iker
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
    integer :: iker
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
  
  subroutine gll2grid_simple(data, grid_data)
    real(kind=cr), dimension(:,:,:,:), intent(in) :: data
    real(kind=cr), dimension(:,:,:), allocatable, intent(out) :: grid_data
    real(kind=cr), dimension(:), allocatable :: tmp, tmp_weight, dat_sum, weight_sum
    real(kind=cr) :: wx, wy, wz, wt
    integer :: i, j, k, ier, ispec, idx, idy, idz, n, m, iglob

    tmp = zeros(ext_grid%nx * ext_grid%ny * ext_grid%nz)
    tmp_weight = zeros(ext_grid%nx * ext_grid%ny * ext_grid%nz)
    dat_sum = zeros(ext_grid%nx * ext_grid%ny * ext_grid%nz)
    weight_sum = zeros(ext_grid%nx * ext_grid%ny * ext_grid%nz)

    do ispec = 1, NSPEC_FWAT
      do i = 1, NGLLX
        do j = 1, NGLLY
          do k = 1, NGLLZ
            iglob = ibool_fwat(i, j, k, ispec)
            call locate_bissection(dble(ext_grid%x), ext_grid%nx, dble(xstore_fwat(iglob)), idx)
            if (idx == -1) call exit_mpi(worldrank, 'ERROR GLL2GRID: x is out of boundary')
            wx = (xstore_fwat(iglob) - ext_grid%x(idx)) / (ext_grid%x(idx+1) - ext_grid%x(idx))
            call locate_bissection(dble(ext_grid%y), ext_grid%ny, dble(ystore_fwat(iglob)), idy)
            if (idy == -1) call exit_mpi(worldrank, 'ERROR GLL2GRID: y is out of boundary')
            wy = (ystore_fwat(iglob) - ext_grid%y(idy)) / (ext_grid%y(idy+1) - ext_grid%y(idy))
            call locate_bissection(dble(ext_grid%z), ext_grid%nz, dble(zstore_fwat(iglob)), idz)
            if (idz == -1) call exit_mpi(worldrank, 'ERROR GLL2GRID: z is out of boundary')
            wz = (zstore_fwat(iglob) - ext_grid%z(idz)) / (ext_grid%z(idz+1) - ext_grid%z(idz))
            do n = 1, 8
              select case(n)
                case (1)
                  m = ext_grid%nx * ext_grid%ny * (idz-1) + ext_grid%nx * (idy-1) + idx
                  wt = (1.0_cr - wx) * (1.0_cr - wy) * (1.0_cr - wz)
                case (2)
                  m = ext_grid%nx * ext_grid%ny * (idz-1) + ext_grid%nx * idy + idx
                  wt = (1.0_cr - wx) * wy * (1.0_cr - wz)
                case (3)
                  m = ext_grid%nx * ext_grid%ny * (idz-1) + ext_grid%nx * idy + idx + 1
                  wt = wx * wy * (1.0_cr - wz)
                case (4)
                  m = ext_grid%nx * ext_grid%ny * (idz-1) + ext_grid%nx * (idy-1) + idx + 1
                  wt = wx * (1.0_cr - wy) * (1.0_cr - wz)
                case (5)
                  m = ext_grid%nx * ext_grid%ny * idz + ext_grid%nx * (idy-1) + idx
                  wt = (1.0_cr - wx) * (1.0_cr - wy) * wz
                case (6)
                  m = ext_grid%nx * ext_grid%ny * idz + ext_grid%nx * idy + idx
                  wt = (1.0_cr - wx) * wy * wz
                case (7)
                  m = ext_grid%nx * ext_grid%ny * idz + ext_grid%nx * idy + idx + 1
                  wt = wx * wy * wz
                case (8)
                  m = ext_grid%nx * ext_grid%ny * idz + ext_grid%nx * (idy-1) + idx + 1
                  wt = wx * (1.0_cr - wy) * wz
              end select
              tmp(m) = tmp(m) + wt * data(i, j, k, ispec)
              tmp_weight(m) = tmp_weight(m) + wt
            enddo
          enddo
        enddo
      enddo
    enddo
    call synchronize_all()
    call sum_all_1Darray_cr(tmp, dat_sum, ext_grid%nx * ext_grid%ny * ext_grid%nz)
    call sum_all_1Darray_cr(tmp_weight, weight_sum, ext_grid%nx * ext_grid%ny * ext_grid%nz)

    if (worldrank == 0) then
      ! Normalize the grid data by the weights
      do i = 1, ext_grid%nx * ext_grid%ny * ext_grid%nz
        if (weight_sum(i) > 0.0_cr) then
          dat_sum(i) = dat_sum(i) / weight_sum(i)
        else
          dat_sum(i) = 0.0_cr
        end if
      end do

      ! Reshape the grid data to 3D
      grid_data = reshape(dat_sum, [ext_grid%nx, ext_grid%ny, ext_grid%nz])
    endif
    call synchronize_all()

  end subroutine gll2grid_simple

  subroutine gll2grid(data, grid_data)
    use shared_parameters

    real(kind=cr), dimension(:,:,:,:), allocatable, intent(in) :: data
    real(kind=cr), dimension(:,:,:), allocatable, intent(out) :: grid_data
    real(kind=cr), dimension(:,:,:), allocatable :: model_on_FD_grid
    type(profd)  :: projection_fd
    integer :: ier

    ! Get regular grid properties, and precomputes all interpolation coefficients
    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
    call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

    call compute_interpolation_coeff_FD_SEM(projection_fd, worldrank)

    ! put data from SEM mesh to a regular grid
    call Project_model_SEM2FD_grid(data, grid_data, projection_fd, worldrank)
    ! if (worldrank == 0) then
    !   grid_data = model_on_FD_grid
    ! else
    !   allocate(grid_data(ext_grid%nx, ext_grid%ny, ext_grid%nz))
    ! end if
    ! call bcast_all_cr(grid_data, ext_grid%nx * ext_grid%ny * ext_grid%nz)
  
    call synchronize_all()
  
  end subroutine gll2grid

end module model_grid_data