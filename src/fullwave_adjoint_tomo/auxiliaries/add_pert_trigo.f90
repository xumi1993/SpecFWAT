program add_pert_trigo

  use constants
  use shared_parameters
  use my_mpi             !! module from specfem
  use interpolation_mod, only: trilin_interp
  
  implicit none
  include "precision.h"
  
  integer :: ier,sizeprocs,myrank,ios, i,j,k, ispec,iglob,NSPEC_IRREGULAR
  character(len=MAX_STRING_LEN) :: xyz_infile,model_dir,data_name,&
              local_data_file, prname, out_dir, isutm, pert_str
  real(kind=CUSTOM_REAL),dimension(:,:,:), allocatable :: pert3d
  !!!
  integer :: NSPEC_AB,NGLOB_AB
  integer, dimension(:,:,:,:),allocatable :: ibool
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: xstore,ystore,zstore
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: xrange,yrange,zrange
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: vstore, vstore_new,dlnvs_trigo
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: pertx, perty, pertz
  real(kind=CUSTOM_REAL) :: xmin, hx, ymin, hy, zmin, hz, val,&
                            xmax, ymax, zmax, pert, valz
  double precision :: xmid, ymid, tmpdble, la, lo
  double precision, dimension(:), allocatable :: x, y
  integer :: nx, ny, nz, npertx, nperty, npertz
  logical :: BROADCAST_AFTER_READ
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
  if (command_argument_count() == 6) then
    call get_command_argument(1,xyz_infile) !grid file
    call get_command_argument(2,data_name)
    call get_command_argument(3,model_dir)
    call get_command_argument(4,out_dir) ! npz file
    call get_command_argument(5,isutm)
    call get_command_argument(6,pert_str)
    read(pert_str, *) pert
  else
    write(*,*) 'xadd_pert_gaus grid_file filename in_dir out_dir isutm perturbation'
    stop 'Not enough arguments'
  endif

  ! read points to be interpolated
  open(11,file=xyz_infile,iostat=ios)
  read(11,*,iostat=ios) xmin, ymin, zmin
  read(11,*,iostat=ios) hx, hy, hz
  read(11,*,iostat=ios) npertx, nperty, npertz
  close(11)
  xmax = xmin + hx*npertx
  ymax = ymin + hy*nperty
  zmax = zmin + hz*npertz
  valz = 1000.
  if (trim(isutm) == 'true') then
    val = 0.01
  else
    val = valz
  endif
  nx = (xmax-xmin)/val+1
  ny = (ymax-ymin)/val+1
  nz = (zmax-zmin)/valz+1

  ! Create trigo perturbation mask
  allocate(xrange(nx), yrange(ny), zrange(nz))
  allocate(pertx(nx),perty(ny),pertz(nz))
  allocate(pert3d(nx,ny,nz))

  do i=1, nx
    xrange(i) = xmin+(i-1)*val
    pertx(i) = sin(npertx*PI*(i-1)/(nx-1))
  enddo
  do i= 1, ny
    yrange(i) = ymin+(i-1)*val
    perty(i) = sin(nperty*PI*(i-1)/(ny-1))
  enddo
  do i= 1, nz
    zrange(i) = zmin+(i-1)*valz
    pertz(i) = sin(npertz*PI*(i-1)/(nz-1))
  enddo
  do i=1,nx;do j=1,ny;do k=1,nz
    pert3d(i,j,k) = pert * pertx(i) * perty(j) * pertz(k)
  enddo;enddo;enddo

  ! Convert to UTM coords
  allocate(x(nx),y(ny))
  if (trim(isutm) == 'true') then
    do i=1,nx
      ymid = (ymin+ymax)/2
      tmpdble = dble(xrange(i))
      call utm_geo(tmpdble, ymid, x(i), la, 0)
    enddo
    do i=1,ny
      xmid = (xmin+xmax)/2
      tmpdble = dble(yrange(i))
      call utm_geo(xmid, tmpdble, lo, y(i), 0)
    enddo
  else
    x=dble(xrange)
    y=dble(yrange)
  endif
  call synchronize_all()
  ! read mesh
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

  ! read data files
  allocate(vstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ios)
  if (ios /= 0) call exit_MPI_without_rank('error allocating array 1105')
  write(prname,'(a,i6.6,a)') trim(model_dir)//'/proc',myrank,'_'
  local_data_file = trim(prname) // trim(data_name) // '.bin'
  open(unit = 27,file = trim(local_data_file),status='old', &
          action='read',form ='unformatted',iostat=ios)
  read(27) vstore
  close(27)
  call synchronize_all()

  allocate(vstore_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ios)
  allocate(dlnvs_trigo(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ios)
  do ispec=1,NSPEC_AB
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          if (xstore(iglob)<x(1) .or. xstore(iglob)>x(nx) .or. &
              ystore(iglob)<y(1) .or. ystore(iglob)>y(ny) .or. &
              zstore(iglob)<zmin .or. zstore(iglob)>zmax) then
            ! if(myrank==0) print *, xstore(iglob), x(1), x(nx), ';',&
            !          ystore(iglob), y(1), y(ny), ';', &
            !          zstore(iglob), zmin, zmax 
            dlnvs_trigo(i,j,k,ispec) = 0.
          else
            call trilin_interp(xstore(iglob),ystore(iglob),zstore(iglob),&
                               real(x),real(y),zrange,nx,ny,nz,pert3d,dlnvs_trigo(i,j,k,ispec))
          endif
          vstore_new(i,j,k,ispec)=vstore(i,j,k,ispec)*(1+dlnvs_trigo(i,j,k,ispec))
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
  write(27) dlnvs_trigo
  close(27)

  deallocate(xrange,yrange,zrange)
  deallocate(x,y)
  deallocate(ibool, pert3d)
  deallocate(pertx, perty, pertz)
  deallocate(xstore,ystore,zstore,vstore)

  call finalize_mpi()
  end program add_pert_trigo