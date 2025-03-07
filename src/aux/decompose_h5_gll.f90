program decompose_h5_gll

  use constants
  use shared_parameters
  ! use specfem_par, only: xigll, yigll, zigll, wxgll, wygll, wzgll
  use m_npy
  use projection_on_FD_grid_fwat
  use hdf5_interface
  use utils, only: transpose_3

  implicit none

  integer, parameter :: NARGS = 3, FID=27

  character(len=MAX_STRING_LEN) :: fname,outdir,grid_file,data_filename,key_name,prname_lp,arg(9)
  integer :: ierr, NSPEC_IRREGULAR, i,myrank,sizeprocs
  ! double precision, dimension(:), allocatable :: xfd, yfd, zfd
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: model_on_SEM_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: model_on_FD_grid
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: x, y, z
  ! type(profd)  :: projection_fd
  logical :: BROADCAST_AFTER_READ

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then
    print *
    print *,'Decomposition of regular grid data on volumetric gll points'
    print *
  endif

  ! needs local_path for mesh files
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! parse command line arguments
  if (command_argument_count() /= NARGS) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/decompose_h5_gll key_name input_file output_dir'
      stop 'Please check command line arguments'
    endif
  endif

  ! reads in arguments
  do i = 1, NARGS
    call get_command_argument(i,arg(i))
  enddo

  ! grid_file = arg(1)
  key_name = arg(1)
  fname= arg(2)
  outdir = arg(3)

  ! call compute_interpolation_coeff_FD_SEM(projection_fd, grid_file, myrank)
  ! allocate(model_on_FD_grid(nx_fd_proj, ny_fd_proj, nz_fd_proj))
  if (myrank == 0) then
    ! call run_sys('unzip -o '//trim(fname), ierr)
    ! fname = trim(data_filename)//'.npy'
    ! call read_npy(fname, model_on_FD_grid)
    data_filename = '/'//trim(key_name)
    call h5read(fname, data_filename, model_on_FD_grid)
    model_on_FD_grid = transpose_3(model_on_FD_grid)
    call h5read(fname, '/x', x)
    call h5read(fname, '/y', y)
    call h5read(fname, '/z', z)
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
    allocate(model_on_FD_grid(nx_fd_proj, ny_fd_proj, nz_fd_proj))
  endif
  call bcast_all_cr(model_on_FD_grid, nx_fd_proj*ny_fd_proj*nz_fd_proj)

  ! Get dimensions of current model, stored in proc******_external_mesh.bin
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'
  open(unit=FID,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin', &
       status='old',action='read',form='unformatted',iostat=ierr)
  read(FID) NSPEC_AB
  read(FID) NGLOB_AB
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ierr)
  if (ierr /= 0) call exit_MPI_without_rank('error allocating array 1102')
  if (ierr /= 0) stop 'error allocating array ibool'
  allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ierr)
  if (ierr /= 0) call exit_MPI_without_rank('error allocating array 1103')
  if (ierr /= 0) stop 'error allocating array xstore etc.'

  read(FID) NSPEC_IRREGULAR
  read(FID) ibool
  read(FID) xstore
  read(FID) ystore
  read(FID) zstore
  close(FID)
  call synchronize_all()
  
  allocate(model_on_SEM_mesh(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  if (myrank==0) print *, 'Interpolation from grid to SEM mesh.'
  call Project_model_FD_grid2SEM(model_on_SEM_mesh, model_on_FD_grid, myrank)
  
  write(prname_lp,'(a,i6.6,a)') trim(outdir)//'/proc',myrank,'_'
  open(FID,file=trim(prname_lp)//trim(key_name)//'.bin',status='unknown',form='unformatted',iostat=ierr)
  if (ierr /= 0) stop 'Error writing model file'
  write(FID) model_on_SEM_mesh
  close(FID)
  call synchronize_all()
  if (myrank==0) print *, 'Finish writing SEM mesh to '//trim(prname_lp)
  deallocate(model_on_SEM_mesh, model_on_FD_grid)
  ! Remove temp .npy files
  ! if(myrank == 0) call run_sys('rm -f x.npy y.npy z.npy '//trim(data_filename)//'.npy', ierr)
  call finalize_mpi()

contains
  subroutine read_npy(fname, x)
    implicit none
    character :: head_magic_num, head_major, head_minor 
    character(len=sizeof(magic_str)) :: head_magic_str
    character(len=MAX_STRING_LEN), intent(in) :: fname
    character(len=:), allocatable :: header_str
    integer(4) :: header_len,ier, i
    integer(2) :: header_len_short
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable, intent(inout) :: x

    open(40, file=fname, status='old', access='stream', form='unformatted')
    read(40) head_magic_num
    read(40) head_magic_str
    read(40) head_major
    read(40) head_minor
    read(40) header_len
    allocate(character(header_len) :: header_str)
    read(40, IOSTAT=ier) header_str
    if (ier/=0) then
      ! rewind(40)
      ! do i = 1,7
      !   read(40,*)
      ! enddo
      call fseek(40, 8, 0)
      read(40) header_len_short
      deallocate(header_str)
      allocate(character(header_len_short) :: header_str)
      read(40) header_str
    endif
    read(40) x
    close(40)

  end subroutine read_npy


end program decompose_h5_gll
