module fwat_mpi
  use mpi
  use my_mpi
  use config

  implicit none

  integer, private :: ier
  integer :: my_node_mpi_comm_world

contains
  subroutine init_mpi_tomo()

    integer :: ier
    integer, dimension(:,:), allocatable :: rank_map_loc

  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
    call MPI_INIT(ier)
    if (ier /= 0 ) stop 'Error initializing MPI'

    ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS and BROADCAST_SAME_MESH_AND_MODEL are read before calling world_split()
    ! thus read the parameter file
    call MPI_COMM_RANK(MPI_COMM_WORLD,worldrank,ier)
    ! if (myrank == 0) then
    !   call open_parameter_file_from_master_only(ier)
    !   ! we need to make sure that NUMBER_OF_SIMULTANEOUS_RUNS and BROADCAST_SAME_MESH_AND_MODEL are read
    !   call read_value_integer(NUMBER_OF_SIMULTANEOUS_RUNS, 'NUMBER_OF_SIMULTANEOUS_RUNS', ier)
    !   if (ier /= 0) stop 'Error reading Par_file parameter NUMBER_OF_SIMULTANEOUS_RUNS'
    !   call read_value_logical(BROADCAST_SAME_MESH_AND_MODEL, 'BROADCAST_SAME_MESH_AND_MODEL', ier)
    !   if (ier /= 0) stop 'Error reading Par_file parameter BROADCAST_SAME_MESH_AND_MODEL'
    !   ! close parameter file
    !   call close_parameter_file()
    ! endif

    ! broadcast parameters read from master to all processes
    my_local_mpi_comm_world = MPI_COMM_WORLD
    ! call bcast_all_singlei(NUMBER_OF_SIMULTANEOUS_RUNS)
    ! call bcast_all_singlel(BROADCAST_SAME_MESH_AND_MODEL)
    ! my_local_mpi_comm_for_bcast = MPI_COMM_NULL

  ! create sub-communicators if needed, if running more than one earthquake from the same job
    ! call world_split()
    call world_rank(worldrank)
    call world_size(worldsize)

    call MPI_Comm_split_type(my_local_mpi_comm_world, MPI_COMM_TYPE_SHARED, worldrank, &
                            MPI_INFO_NULL, my_node_mpi_comm_world, ier)
    call MPI_Comm_rank(my_node_mpi_comm_world, noderank, ier)
    call MPI_Comm_size(my_node_mpi_comm_world, nodesize, ier)

    allocate(rank_map_loc(worldsize, 2), rank_map(worldsize, 2))
    rank_map_loc = 0

    rank_map_loc(worldrank+1, 1) = worldrank
    rank_map_loc(worldrank+1, 2) = noderank
    call sum_all_1Darray_i(rank_map_loc, rank_map, worldsize*2)
  
  end subroutine init_mpi_tomo

  subroutine sum_all_1Darray_i(sendbuf, recvbuf, nx)

    integer :: nx
    integer, dimension(nx) :: sendbuf, recvbuf

    ! checks if anything to do
    if (nx == 0) return

    call MPI_REDUCE(sendbuf,recvbuf,nx,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_1Darray_i

  subroutine prepare_shm_array_dp_1d(buffer, n_elem, win)
    USE, INTRINSIC :: ISO_C_BINDING
    real(kind=dp), dimension(:), pointer :: buffer
    integer :: ierr, istat,n_elem,n
    integer(kind=MPI_ADDRESS_KIND) :: size
    integer :: win, real_size
    type(C_PTR) :: c_window_ptr
    
    ! call world_rank(myrank)
    n = n_elem
    if(noderank /= 0) n = 0
    CALL MPI_Type_size(MPI_DOUBLE_PRECISION, real_size, ierr)
    size = n * real_size
    call MPI_Win_allocate_shared(size, real_size, MPI_INFO_NULL, MPI_COMM_WORLD, c_window_ptr, win, ierr)
    if (noderank /= 0) then
      call MPI_Win_shared_query(win, 0, size, real_size, c_window_ptr, ierr)
    endif
    CALL C_F_POINTER(c_window_ptr, buffer, SHAPE = [n_elem])
    call MPI_Win_fence(0, win, ierr)

  end subroutine prepare_shm_array_dp_1d

  subroutine prepare_shm_array_dp_3d(buffer, nx, ny, nz, win)
    USE, INTRINSIC :: ISO_C_BINDING
    double precision, dimension(:,:,:), pointer :: buffer
    integer :: ierr,nx,ny,nz,n
    integer(kind=MPI_ADDRESS_KIND) :: size
    integer :: win, real_size
    type(C_PTR) :: c_window_ptr
    
    n = nx*ny*nz
    if(noderank /= 0) n = 0
    CALL MPI_Type_size(MPI_DOUBLE_PRECISION, real_size, ierr)
    size = n * real_size
    call MPI_Win_allocate_shared(size, real_size, MPI_INFO_NULL, MPI_COMM_WORLD, c_window_ptr, win, ierr)
    if (noderank /= 0) then
      call MPI_Win_shared_query(win, 0, size, real_size, c_window_ptr, ierr)
    endif
    CALL C_F_POINTER(c_window_ptr, buffer, SHAPE=[nx, ny, nz])
    call MPI_Win_fence(0, win, ierr)

  end subroutine prepare_shm_array_dp_3d

  subroutine prepare_shm_array_ch_1d(buffer, n_elem, nlen, win)
    USE, INTRINSIC :: ISO_C_BINDING
    integer :: ierr, istat,n_elem,n,nlen
    character(len=nlen), dimension(:), pointer :: buffer
    integer(kind=MPI_ADDRESS_KIND) :: size
    integer :: win, char_size
    type(C_PTR) :: c_window_ptr
    
    n = n_elem*nlen
    if(noderank /= 0) n = 0
    CALL MPI_Type_size(MPI_CHARACTER, char_size, ierr)
    size = n * char_size
    call MPI_Win_allocate_shared(size, char_size, MPI_INFO_NULL, MPI_COMM_WORLD, c_window_ptr, win, ierr)
    if (noderank /= 0) then
      call MPI_Win_shared_query(win, 0, size, char_size, c_window_ptr, ierr)
    endif
    CALL C_F_POINTER(c_window_ptr, buffer, SHAPE = [n_elem])
    call MPI_Win_fence(0, win, ierr)

  end subroutine prepare_shm_array_ch_1d

end module fwat_mpi