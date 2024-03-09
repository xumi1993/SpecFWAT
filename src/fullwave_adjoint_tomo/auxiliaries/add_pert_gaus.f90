program add_pert_gaus

  use constants
  use shared_parameters
  use my_mpi             !! module from specfem
  use shared_input_parameters, only: SUPPRESS_UTM_PROJECTION
  
  implicit none
  include "precision.h"
  
  integer :: ier,sizeprocs,myrank,ios, i,j, k,n, ispec,iglob,ipt, npts,NSPEC_IRREGULAR
  character(len=MAX_STRING_LEN) :: xyz_infile,model_dir,data_name,&
              local_data_file, prname, out_dir, isutm, pert_str
  double precision,dimension(:), allocatable :: x, y, dist
  !!!
  integer :: NSPEC_AB,NGLOB_AB
  integer, dimension(:,:,:,:),allocatable :: ibool
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: xstore,ystore,zstore
  double precision, dimension(:),allocatable :: xrange,yrange,zrange
  double precision, dimension(:,:), allocatable :: xyzgrid
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: vstore, vstore_new, dlnvs_gaus
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: pertx, perty, pertz, pert_grid
  real(kind=CUSTOM_REAL) :: max_dlnv, max_dlnv_all
  double precision :: xmin, hx, ymin, hy, zmin, hz,la,lo, polar, pert, sigma,&
                     gaus, xmid,ymid, sigma_b, sigma_c, theta, alpha
  integer :: nx, ny, nz
  logical :: BROADCAST_AFTER_READ
  double precision, dimension(3,3) :: rot_mat 

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
  call get_command_argument(1,xyz_infile) !grid file
  call get_command_argument(2,data_name)
  call get_command_argument(3,model_dir)
  call get_command_argument(4,out_dir) ! npz file
  call get_command_argument(5,isutm)
  call get_command_argument(6,pert_str)
  read(pert_str, *) pert
  call get_command_argument(7,pert_str)
  read(pert_str, *) sigma
  if (command_argument_count() == 7) then
    sigma_b = sigma
    sigma_c = sigma
    theta = 0.
    alpha = 0.
  elseif (command_argument_count() == 11) then
    call get_command_argument(8, pert_str)
    read(pert_str, *) sigma_b
    call get_command_argument(9, pert_str)
    read(pert_str, *) sigma_c
    call get_command_argument(10, pert_str)
    read(pert_str, *) theta
    call get_command_argument(11, pert_str)
    read(pert_str, *) alpha
  else
    write(*,*) 'xadd_pert_gaus grid_file filename in_dir out_dir isutm perturbation sigma [b, c, inc, azi]'
    stop 'Not enough or too more arguments'
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
  allocate(x(nx),y(ny),dist(npts))
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
  allocate(pertx(nx),perty(ny),pertz(nz))
  polar = 1.0
  do i=1,nx
    pertx(i) = polar
    polar = -polar
  enddo
  polar = 1.0
  do i=1,ny
    perty(i) = polar
    polar = -polar
  enddo
  polar = 1.0
  do i=1,nz
    pertz(i) = polar
    polar = -polar
  enddo
  allocate(pert_grid(npts), xyzgrid(npts, 3))
  n = 1
  do k=1,nz
    do j=1,ny
      do i=1,nx
        pert_grid(n) = pertx(i) * perty(j) * pertz(k)
        xyzgrid(n, 1) = x(i)
        xyzgrid(n, 2) = y(j)
        xyzgrid(n, 3) = zrange(k)
        n = n + 1
      enddo
    enddo
  enddo

  rot_mat(1, 1) = cosd(theta)
  rot_mat(1, 2) = -sind(theta)*sind(alpha)
  rot_mat(1, 3) = -sind(theta)*cosd(alpha)
  rot_mat(2, 1) = sind(theta)
  rot_mat(2, 2) = cosd(theta)*sind(alpha)
  rot_mat(2, 3) = cosd(theta)*cosd(alpha)
  rot_mat(3, 1) = 0.
  rot_mat(3, 2) = -cosd(alpha)
  rot_mat(3, 3) = sind(alpha)

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

  allocate(vstore_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ios)
  allocate(dlnvs_gaus(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ios)
  vstore_new(:,:,:,:)=0.
  dlnvs_gaus(:,:,:,:)=0.
  do ipt=1,npts ! Number of perturbation points
    do ispec=1,NSPEC_AB
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            ! calculate distance
            call calc_dist(xyzgrid(ipt,1), xyzgrid(ipt,2), xyzgrid(ipt,3), &
                          dble(xstore(iglob)), dble(ystore(iglob)), dble(zstore(iglob)),&
                          sigma, sigma_b, sigma_c, rot_mat, dist(ipt))
            gaus=exp(-1*dist(ipt)*dist(ipt)/2)
            ! dist(ipt)=dsqrt((xyzgrid(ipt,1)-dble(xstore(iglob)))**2 &
            !             +(xyzgrid(ipt,2)-dble(ystore(iglob)))**2 &
            !             +(xyzgrid(ipt,3)-dble(zstore(iglob)))**2)
            ! gaus=exp(-1*dist(ipt)*dist(ipt)/2/sigma/sigma)
            ! if (gaus > 0.01) print *, gaus, dist(ipt)
            ! if (myrank==0) print *, xyzgrid(ipt, :), pert_grid(ipt), gaus, dist(ipt)
            ! vstore_new(i,j,k,ispec)=vstore(i,j,k,ispec)*(1+pert_grid(ipt)*gaus)
            dlnvs_gaus(i,j,k,ispec)=dlnvs_gaus(i,j,k,ispec)+pert_grid(ipt)*gaus
          enddo
        enddo
      enddo
    enddo
  enddo

  max_dlnv = maxval(abs(dlnvs_gaus(:,:,:,:)))
  call synchronize_all()
  call max_all_cr(max_dlnv,max_dlnv_all)
  call bcast_all_singlecr(max_dlnv_all)
  dlnvs_gaus = pert*dlnvs_gaus/max_dlnv_all

  do ispec=1,NSPEC_AB
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          vstore_new(i,j,k,ispec)=vstore(i,j,k,ispec)*(1+dlnvs_gaus(i,j,k,ispec))
          enddo
      enddo
    enddo
  enddo ! end of loop on all the elements in current slice

  call system('mkdir -p '//trim(out_dir))
  write(prname,'(a,i6.6,a)') trim(out_dir)//'/proc',myrank,'_'//trim(data_name)//'.bin'
  ! if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'
  open(27,file=trim(prname),form='unformatted',action='write')
  write(27) vstore_new
  close(27)

  write(prname,'(a,i6.6,a)') trim(out_dir)//'/proc',myrank,'_dln'//trim(data_name)//'.bin'
  ! if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'
  open(27,file=trim(prname),form='unformatted',action='write')
  write(27) dlnvs_gaus
  close(27)

  deallocate(xrange,yrange,zrange)
  deallocate(x,y,dist)
  deallocate(ibool, xyzgrid, pert_grid)
  deallocate(pertx, perty, pertz)
  deallocate(xstore,ystore,zstore,vstore)

  call finalize_mpi()
  end program add_pert_gaus

  subroutine calc_dist(x1, y1, z1, x2, y2, z2, a, b, c, rot_mat, dis)
    implicit none
    double precision, intent(in) :: x1, y1, z1, x2, y2, z2, a, b, c
    double precision, intent(out) :: dis
    double precision, dimension(3) :: delta, delta_hat
    double precision, dimension(3,3) :: rot_mat 

    delta =(/x1-x2, y1-y2, z1-z2/)
    delta_hat = matmul(rot_mat, delta)
    dis = dsqrt((delta_hat(1)/a)**2+(delta_hat(2)/b)**2+(delta_hat(3)/c)**2)
    
  end subroutine calc_dist