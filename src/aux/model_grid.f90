program sem_model_slice

  use constants
  use shared_parameters
  use m_npy
  use my_mpi             !! module from specfem
  use hdf5_interface
  use shared_input_parameters, only: SUPPRESS_UTM_PROJECTION
  use projection_on_FD_grid_fwat, only: find_nearestXY_nonzero,find_nearestZ_nonzero
  
  implicit none
  include "precision.h"
  
  integer, parameter :: NMAXPTS = 10000000
  integer :: ier,sizeprocs,myrank,ios, i,j, k, ispec,iglob,ipt, npts,NSPEC_IRREGULAR
  character(len=MAX_STRING_LEN) :: xyz_infile,model_dir,data_name,gmt_outfile,&
             local_data_file, prname, out_dir,suffix, isutm, replace_zero
  ! real(kind=CUSTOM_REAL),dimension(NMAXPTS) : v, &
            !  distmin,distall
  !  double precision,   dimension(NMAXPTS) :: dis
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: dist,vall,vmin,v,distmin
  double precision,dimension(:), allocatable :: x, y
  ! integer,dimension(NMAXPTS) :: ispec_min
  real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: in, out
  integer, dimension(:), allocatable :: grid_count, allcount,ispec_min
  !!!
  integer :: NSPEC_AB,NGLOB_AB
  integer, dimension(:,:,:,:),allocatable :: ibool
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: xstore,ystore,zstore
  double precision, dimension(:),allocatable :: xrange,yrange,zrange
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: v3d
  integer, dimension(:,:,:), allocatable :: igrid
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: vstore
  double precision :: xmin, hx, ymin, hy, zmin, hz,la,lo, ymid, xmid,dis,mindis
  integer :: nx, ny, nz, count,ix_min,iy_min,iz_min,ii,jj,kk,iptmin
  integer, dimension(2) :: ix, iy, iz
  logical :: BROADCAST_AFTER_READ
  type(hdf5_file) :: h5file
  !!!

  ! true --> replace "air" points with NaN (vertical cross sections)
  ! false --> take the closest value to the "air" points (horizontal cross section)
  ! integer, dimension(NX_TOPO_FILE,NY_TOPO_FILE) :: itopo_bathy_basin
  ! double precision :: elevation

  ! MPI initialization
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! input arguments
  if (command_argument_count() >= 5) then
    call getarg(1,xyz_infile) !grid file
    !   call getarg(2,topo_dir)
    call getarg(2,data_name)
    call getarg(3,model_dir)
    call getarg(4,out_dir) ! npz file
    call getarg(5,isutm)
  else
    stop 'Not enough arguments'
  endif
  if (command_argument_count() >= 6) then
    call getarg(6,suffix)
  elseif (command_argument_count() == 5) then
    suffix = 'grid'
  endif
  if (command_argument_count() ==7) then
    call getarg(7,replace_zero)
  elseif (command_argument_count() == 6) then
    replace_zero = 'true'
  else
    stop 'Too more arguments'
  endif
  ! read points to be interpolated
  open(11,file=xyz_infile,iostat=ios)
  read(11,*,iostat=ios) xmin, ymin, zmin
  read(11,*,iostat=ios) hx, hy, hz
  read(11,*,iostat=ios) nx, ny, nz
  close(11)
  allocate(xrange(nx),yrange(ny),zrange(nz))
  do i= 1, nx
    xrange(i) = xmin+(i-1)*hx
  enddo
  do i= 1, ny
    yrange(i) = ymin+(i-1)*hy
  enddo
  do i= 1, nz
    zrange(i) = zmin+(i-1)*hz
  enddo
  npts=nx*ny*nz
  allocate(x(nx),y(ny),dist(npts),grid_count(npts),allcount(npts),&
           vmin(npts), vall(npts),distmin(npts),v(npts),ispec_min(npts))
  allocate(igrid(nx,ny,nz))
  if (trim(isutm) == 'true') then
    do i=1,nx
      ymid = (ymin+yrange(ny))/2
      call utm_geo(xrange(i), ymid, x(i), la, 0)
    enddo
    do i=1,ny
      xmid = (xmin+xrange(nx))/2
      call utm_geo(xmid, yrange(i), lo, y(i), 0)
    enddo
  else
    x=xrange
    y=yrange
  endif
  ! z = zrange
  count = 0
  do i=1,nx
    do j=1,ny
      do k=1,nz
        count = count + 1
        igrid(i,j,k) = count
      enddo
    enddo
  enddo
!   if (myrank == 0) then
!     write(*,*) 'Total number of points = ', npts
!     if (npts > NMAXPTS .or. npts <= 0) call exit_mpi(myrank,'Npts error ...')
!   endif

  !Kai added read mesh
  write(prname,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'
  open(unit=27,file=trim(prname)//'external_mesh.bin', &
        status='old',action='read',form='unformatted',iostat=ier)
  read(27) NSPEC_AB
  read(27) NGLOB_AB

  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ios)
  if (ios /= 0) call exit_MPI_without_rank('error allocating array 1102')
  if (ios /= 0) stop 'error allocating array ibool'
  allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ios)
  if (ios /= 0) call exit_MPI_without_rank('error allocating array 1103')
  if (ios /= 0) stop 'error allocating array xstore etc.'

  read(27) NSPEC_IRREGULAR
  read(27) ibool
  read(27) xstore
  read(27) ystore
  read(27) zstore
  close(27)

  ! read data and topo files
  allocate(vstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ios)
  if (ios /= 0) call exit_MPI_without_rank('error allocating array 1105')
  write(prname,'(a,i6.6,a)') trim(model_dir)//'/proc',myrank,'_'
  local_data_file = trim(prname) // trim(data_name) // '.bin'
  open(unit = 27,file = trim(local_data_file),status='old', &
        action='read',form ='unformatted',iostat=ios)
  read(27) vstore
  close(27)
  call synchronize_all()

  ! search for local minimum-distance point
  distmin(:) = HUGEVAL
  ! grid_count = 0
  vmin = 0.
  do ispec=1,NSPEC_AB
    ! print *, ispec, '/', NSPEC_AB

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ix = nearest_two_points(x, dble(xstore(iglob)))
          iy = nearest_two_points(y, dble(ystore(iglob)))
          iz = nearest_two_points(zrange, dble(zstore(iglob)))
          do ii = 1,2
            do jj = 1,2
              do kk = 1,2
                ipt = igrid(ix(ii), iy(jj), iz(kk))
                dis = dsqrt((x(ix(ii))-dble(xstore(iglob)))**2 &
                                 +(y(iy(jj))-dble(ystore(iglob)))**2 &
                                 +(zrange(iz(kk))-dble(zstore(iglob)))**2)
                ! grid_count(ipt) = grid_count(ipt)+1
                ! vmin(ipt) = vmin(ipt) + vstore(i,j,k,ispec)
                if(dis < distmin(ipt)) then
                  distmin(ipt)=dis
                  ispec_min(ipt)=ispec
                  vmin(ipt)=vstore(i,j,k,ispec)
                endif
              enddo
            enddo
          enddo
          ! grid_count(iptmin) = grid_count(iptmin)+1
          ! vmin(iptmin) = vmin(iptmin) + vstore(i,j,k,ispec)
        enddo
      enddo
    enddo

    ! end of loop on all the elements in current slice
  enddo

  ! frees memory
  deallocate(ibool)
  deallocate(xstore,ystore,zstore,vstore)


  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  if (myrank == 0) print *, 'Done looping over global points ...'

  ! call sum_all_1Darray_dp(vmin,vall,npts)
  ! call MPI_REDUCE(grid_count, allcount,npts,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD,ier)


  allocate(in(2,npts), out(2,npts))
  do i=1, npts
    in(1,i) = distmin(i)
    in(2,i) = myrank    ! myrank is coerced to a double
  enddo
  call MPI_REDUCE(in,out,npts,MPI_2REAL,MPI_MINLOC,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(out,2*npts,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  v(1:npts) = 0
  dist(1:npts) = 0.
  vall = 0
  do i = 1, npts
    if (myrank == nint(out(2,i))) then
      v(i) = vmin(i)
      dist(i) = distmin(i)
    endif
  enddo
  call MPI_REDUCE(v,vall,npts,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
  ! call MPI_REDUCE(dist,distall,npts,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

 
  if (myrank == 0) then
    ! do i = 1,npts
    !   if(allcount(i)>0)  vall(i) = vall(i)/allcount(i)
    ! enddo
    ! zrange = zrange / 1000.
    allocate(v3d(nx,ny,nz))
    count=0
    do i=1,nx
      do j=1,ny
        do k=1,nz
          count = count + 1
          v3d(i,j,k) = vall(count)
        enddo
      enddo
    enddo
    if(replace_zero=='true') then
      do i=1,nx
        do j=1,ny
          do k=1,nz
            if (v3d(i,j,k) == 0) then
              if (zrange(k)>0) then
                call find_nearestZ_nonzero(v3d,i,j,k,nx,ny,nz)
              else
                call find_nearestXY_nonzero(v3d,i,j,k,nx,ny,nz)
              endif
            endif
          enddo
        enddo
      enddo
    endif
    gmt_outfile =  trim(out_dir)//'/'//trim(data_name) // '_'//trim(suffix)//'.h5'
    call h5file%open(gmt_outfile, status='new', action='write')
    call h5file%add('/x', xrange)
    call h5file%add('/y', yrange)
    call h5file%add('/z', zrange)
    call h5file%add('/'//trim(data_name), transpose_3(v3d))
    deallocate(v3d)
  endif
  deallocate(xrange,yrange,zrange)
  deallocate(x,y,dist, igrid)
  call finalize_mpi()

contains

function nearest_two_points(xx, newx) result(idx)
  implicit none
  double precision, dimension(:),allocatable :: xx
  double precision :: newx
  integer :: i, ixx(1), nx
  integer :: idx(2)

  nx = size(xx)
  ixx = minloc(abs(xx-newx))
  i = ixx(1)
  if (xx(i)-newx<=0) then
    idx(1) = i
    if (i==nx) then
      idx(2) = i
    else
      idx(2) = i+1
    endif
  else
    if (i==1) then
      idx(1) = i
    else
      idx(1) = i-1
    endif
    idx(2) = i
  endif
  return
end function nearest_two_points

end program sem_model_slice
