module MeshGrid
  use constants
  use shared_parameters
  use specfem_par
  use utils
  use projection_on_FD_grid_fwat
  use hdf5_interface

  implicit none
  
  type, public :: ReglGrid
    real(kind=CUSTOM_REAL) :: xstart, xend, ystart, yend, zstart, zend, dx, dy, dz
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: xfd, yfd, zfd
    integer :: nx, ny, nz
    contains
    procedure :: init, griddata, read_mesh, gridmodel, semdata
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

  do i=1, this%nx
    do j=1, this%ny
      do k=1, this%nz
        if (model_on_FD_grid(i,j,k) == 0) then
          if (this%zfd(k)>0) then
            call find_nearestZ_nonzero(model_on_FD_grid,i,j,k,this%nx,this%ny,this%nz) 
          else
            call find_nearestXY_nonzero(model_on_FD_grid,i,j,k,this%nx,this%ny,this%nz)
          endif
        endif
      enddo
    enddo
  enddo

end subroutine gridmodel

subroutine semdata(this, indir, dataname, model_on_SEM_mesh)
  class (ReglGrid), intent(inout) :: this

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: model_on_FD_grid
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable, intent(out) :: model_on_SEM_mesh
  character(len=*), intent(in) :: indir, dataname

  if (myrank == 0) then
    call h5read(trim(indir)//'/'//trim(dataname)//'.h5', '/vp', model_on_FD_grid)
  endif
  if (myrank /= 0) then
    allocate(model_on_FD_grid(nx_fd_proj, ny_fd_proj, nz_fd_proj))
  endif
  call bcast_all_cr(model_on_FD_grid, nx_fd_proj*ny_fd_proj*nz_fd_proj)

  allocate(model_on_SEM_mesh(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  call Project_model_FD_grid2SEM(model_on_SEM_mesh, model_on_FD_grid, myrank)

end subroutine semdata

subroutine togllmodel(nspec)
  use constants
  use utils
  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IMAIN,MAX_STRING_LEN,ATTENUATION,&
                                    LOCAL_PATH, xstore_db=>xstore, ystore_db=>ystore, zstore_db=>zstore,&
                                    ibool_db=>ibool
  use specfem_par, only: myrank
  use create_regions_mesh_ext_par, only: rhostore,kappastore,mustore,rho_vp,rho_vs,qkappa_attenuation_store,qmu_attenuation_store
  use projection_on_FD_grid_fwat
  use fullwave_adjoint_tomo_par, only: GRID_PATH

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: vp_read,vs_read,rho_read, qmu_read, qkappa_read
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: vp_gll, vs_gll, rho_gll
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: x, y, z
  real(kind=CUSTOM_REAL) :: ox_back, hx_back, oy_back, hy_back, oz_back, hz_back
  integer :: nx_back, ny_back, nz_back
  character(len=MAX_STRING_LEN) :: prname_lp
  integer :: ier, FID=28, nspec, i, j, k, ispec
  logical :: exist


  ! backup the original grid
  ox_back = ox_fd_proj; hx_back = hx_fd_proj; nx_back = nx_fd_proj
  oy_back = oy_fd_proj; hy_back = hy_fd_proj; ny_back = ny_fd_proj
  oz_back = oz_fd_proj; hz_back = hz_fd_proj; nz_back = nz_fd_proj

  ! read h5file
  if (myrank == 0) then
    write(IMAIN,*) '     reading in: ',trim(GRID_PATH)
    inquire(file=trim(GRID_PATH), exist=exist)
    if (.not. exist) then
      print *, 'ERROR: No found mesh file of ', trim(GRID_PATH)
      stop
    endif
    call h5read(GRID_PATH, '/vp', vp_read)
    vp_read = transpose_3(vp_read)
    call h5read(GRID_PATH, '/vs', vs_read)
    vs_read = transpose_3(vs_read)
    call h5read(GRID_PATH, '/rho', rho_read)
    rho_read = transpose_3(rho_read)
    if (ATTENUATION) then
      call h5read(GRID_PATH, '/qmu', qmu_read)
      qmu_read = transpose_3(qmu_read)
      call h5read(GRID_PATH, '/qkappa', qkappa_read)
      qkappa_read = transpose_3(qkappa_read)
    endif
    call h5read(GRID_PATH, '/x', x)
    call h5read(GRID_PATH, '/y', y)
    call h5read(GRID_PATH, '/z', z)
    ox_fd_proj = x(1); hx_fd_proj = x(2) - x(1); nx_fd_proj = size(x)
    oy_fd_proj = y(1); hy_fd_proj = y(2) - y(1); ny_fd_proj = size(y)
    oz_fd_proj = z(1); hz_fd_proj = z(2) - z(1); nz_fd_proj = size(z)
  endif
  call synchronize_all()
  call bcast_all_singlei(nx_fd_proj)
  call bcast_all_singlei(ny_fd_proj)
  call bcast_all_singlei(nz_fd_proj)
  call bcast_all_singlecr(ox_fd_proj)
  call bcast_all_singlecr(hx_fd_proj)
  call bcast_all_singlecr(oy_fd_proj)
  call bcast_all_singlecr(hy_fd_proj)
  call bcast_all_singlecr(oz_fd_proj)
  call bcast_all_singlecr(hz_fd_proj)
  if (myrank /= 0) then
    allocate(vp_read(nx_fd_proj, ny_fd_proj, nz_fd_proj))
    allocate(vs_read(nx_fd_proj, ny_fd_proj, nz_fd_proj))
    allocate(rho_read(nx_fd_proj, ny_fd_proj, nz_fd_proj))
  endif
  call bcast_all_cr(vp_read, nx_fd_proj*ny_fd_proj*nz_fd_proj)
  call bcast_all_cr(vs_read, nx_fd_proj*ny_fd_proj*nz_fd_proj)
  call bcast_all_cr(rho_read, nx_fd_proj*ny_fd_proj*nz_fd_proj)
  if (ATTENUATION) then
    if (myrank /= 0) then
      allocate(qmu_read(nx_fd_proj, ny_fd_proj, nz_fd_proj))
      allocate(qkappa_read(nx_fd_proj, ny_fd_proj, nz_fd_proj))
    endif
    call bcast_all_cr(qmu_read, nx_fd_proj*ny_fd_proj*nz_fd_proj)
    call bcast_all_cr(qkappa_read, nx_fd_proj*ny_fd_proj*nz_fd_proj)
  endif

  allocate(vp_gll(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vs_gll(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho_gll(NGLLX,NGLLY,NGLLZ,nspec))

  ibool = ibool_db
  NSPEC_AB = nspec
  xstore = zeros(NGLLX*NGLLY*NGLLZ*nspec)
  ystore = zeros(NGLLX*NGLLY*NGLLZ*nspec)
  zstore = zeros(NGLLX*NGLLY*NGLLZ*nspec)
  do i=1,NGLLX; do j=1,NGLLY; do k=1,NGLLZ; do ispec=1,nspec
    xstore(ibool(i,j,k,ispec)) = xstore_db(i,j,k,ispec)
    ystore(ibool(i,j,k,ispec)) = ystore_db(i,j,k,ispec)
    zstore(ibool(i,j,k,ispec)) = zstore_db(i,j,k,ispec)
  enddo; enddo; enddo; enddo
  call Project_model_FD_grid2SEM(vp_gll, vp_read, myrank)
  call Project_model_FD_grid2SEM(vs_gll, vs_read, myrank)
  call Project_model_FD_grid2SEM(rho_gll, rho_read, myrank)

  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'
  open(FID,file=trim(prname_lp)//'vp.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error writing model file'
  write(FID) vp_gll
  close(FID)

  open(FID,file=trim(prname_lp)//'vs.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error writing model file'
  write(FID) vs_gll
  close(FID)

  open(FID,file=trim(prname_lp)//'rho.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error writing model file'
  write(FID) rho_gll
  close(FID)

  if (ATTENUATION) then
    call Project_model_FD_grid2SEM(qmu_attenuation_store, qmu_read, myrank)
    call Project_model_FD_grid2SEM(qkappa_attenuation_store, qkappa_read, myrank)

    open(FID,file=trim(prname_lp)//'qmu.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error writing model file'
    write(FID) qmu_attenuation_store
    close(FID)

    open(FID,file=trim(prname_lp)//'qkappa.bin',status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error writing model file'
    write(FID) qkappa_attenuation_store
    close(FID)
  endif

  deallocate(vp_gll, vs_gll, rho_gll)
  deallocate(vp_read, vs_read, rho_read)
  if (ATTENUATION) then
    deallocate(qmu_read, qkappa_read)
  endif

  ! restore the original grid
  ox_fd_proj = ox_back; hx_fd_proj = hx_back; nx_fd_proj = nx_back
  oy_fd_proj = oy_back; hy_fd_proj = hy_back; ny_fd_proj = ny_back
  oz_fd_proj = oz_back; hz_fd_proj = hz_back; nz_fd_proj = nz_back

  call synchronize_all()


end subroutine togllmodel



end module MeshGrid