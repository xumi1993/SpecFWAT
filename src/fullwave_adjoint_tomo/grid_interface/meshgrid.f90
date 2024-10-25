module MeshGrid
  use constants
  use shared_parameters
  use specfem_par
  use utils
  use projection_on_FD_grid_fwat

  implicit none
  
  type, public :: ReglGrid
    real(kind=CUSTOM_REAL) :: xstart, xend, ystart, yend, zstart, zend, dx, dy, dz
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: xfd, yfd, zfd
    integer :: nx, ny, nz
    contains
    procedure :: init, griddata, read_mesh, gridmodel
  end type

contains

subroutine init(this, xstart, xend, ystart, yend, zstart, zend, dx, dy, dz)
  class (ReglGrid), intent(inout) :: this
  real(kind=CUSTOM_REAL), intent(in) :: xstart, xend, ystart, yend, zstart, zend, dx, dy, dz
  integer :: i

  this%xstart = xstart
  this%xend = xend
  this%ystart = ystart
  this%yend = yend
  this%zstart = zstart
  this%zend = zend
  this%dx = dx
  this%dy = dy
  this%dz = dz

  hx_fd_proj = dx
  hy_fd_proj = dy
  hz_fd_proj = dz

  ox_fd_proj = xstart
  oy_fd_proj = ystart
  oz_fd_proj = zstart
  nx_fd_proj = nint((xend - xstart) / hx_fd_proj) + 1
  ny_fd_proj = nint((yend - ystart) / hy_fd_proj) + 1
  nz_fd_proj = nint((zend - zstart) / hz_fd_proj) + 1

  this%nx = nx_fd_proj
  this%ny = ny_fd_proj
  this%nz = nz_fd_proj

  allocate(this%xfd(nx_fd_proj))
  allocate(this%yfd(ny_fd_proj))
  allocate(this%zfd(nz_fd_proj))

  do i=1,nx_fd_proj
    this%xfd(i) = ox_fd_proj + (i-1) * hx_fd_proj
  enddo
  do i=1,ny_fd_proj
    this%yfd(i) = oy_fd_proj + (i-1) * hy_fd_proj
  enddo
  do i=1,nz_fd_proj
    this%zfd(i) = oz_fd_proj + (i-1) * hz_fd_proj
  enddo
end subroutine init

subroutine read_mesh(this)
  class (ReglGrid), intent(inout) :: this
  character(len=MAX_STRING_LEN) :: prname, prname_lp, local_data_file
  integer :: NSPEC_IRREGULAR, ier

 if (myrank==0) print *, 'Reading GLL mesh...'
  ! Get dimensions of current model, stored in proc******_external_mesh.bin
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'
  open(unit=27,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  read(27) NSPEC_AB
  read(27) NGLOB_AB
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1102')
  if (ier /= 0) stop 'error allocating array ibool'
  allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1103')
  if (ier /= 0) stop 'error allocating array xstore etc.'

  read(27) NSPEC_IRREGULAR
  read(27) ibool
  read(27) xstore
  read(27) ystore
  read(27) zstore
  close(27)
  call synchronize_all()
end subroutine read_mesh


subroutine griddata(this, indir, dataname, model_on_FD_grid)
  class (ReglGrid), intent(inout) :: this

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable, intent(out) :: model_on_FD_grid
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: data_sp
  character(len=*), intent(in) :: indir, dataname
  character(len=MAX_STRING_LEN) :: prname, prname_lp, local_data_file
  integer :: ier
  type(profd)  :: projection_fd

  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  call compute_interpolation_coeff_FD_SEM(projection_fd, myrank)
  allocate(model_on_FD_grid(projection_fd%nx, projection_fd%ny, projection_fd%nz),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1104')

  ! Get data to project
  allocate(data_sp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1105')
  if (ier /= 0) stop 'error allocating single precision data array'

  ! data file
  write(prname,'(a,i6.6,a)') trim(indir)//'/proc',myrank,'_'
  local_data_file = trim(prname) // trim(dataname) // '.bin'
  open(unit = 28,file = trim(local_data_file),status='old', &
        action='read',form ='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening ',trim(local_data_file)
    stop
  endif
  read(28) data_sp

  ! put data from SEM mesh to a regular grid
  call Project_model_SEM2FD_grid(data_sp, model_on_FD_grid, projection_fd, myrank)
  call synchronize_all()

end subroutine griddata

subroutine gridmodel(this, indir, dataname, model_on_FD_grid)
  class (ReglGrid), intent(inout) :: this

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable, intent(out) :: model_on_FD_grid
  character(len=*), intent(in) :: indir, dataname
  integer :: i, j, k

  call this%griddata(indir, dataname, model_on_FD_grid)

  do i=1, rg%nx
    do j=1, rg%ny
      do k=1, rg%nz
        if (model_on_FD_grid(i,j,k) == 0) then
          if (rg%zfd(k)>0) then
            call find_nearestZ_nonzero(model_on_FD_grid,i,j,k,rg%nx,rg%ny,rg%nz) 
          else
            call find_nearestXY_nonzero(model_on_FD_grid,i,j,k,rg%nx,rg%ny,rg%nz)
          endif
        endif
      enddo
    enddo
  enddo

end subroutine gridmodel

subroutine togllmodel()

end subroutine togllmodel



end module MeshGrid