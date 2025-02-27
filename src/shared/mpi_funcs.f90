module fwat_mpi
  use mpi
  use my_mpi
  use config

  implicit none

  integer, private :: ier
  integer :: my_node_mpi_comm_world

contains
  subroutine init_mpi_fwat()

    integer :: ier
    integer, dimension(:,:), allocatable :: rank_map_loc

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
  
  end subroutine init_mpi_fwat

  subroutine sum_all_1Darray_i(sendbuf, recvbuf, nx)

    integer :: nx
    integer, dimension(nx) :: sendbuf, recvbuf

    ! checks if anything to do
    if (nx == 0) return

    call MPI_REDUCE(sendbuf,recvbuf,nx,MPI_INTEGER,MPI_SUM,0,my_local_mpi_comm_world,ier)

  end subroutine sum_all_1Darray_i

  subroutine send_ch_array(sendbuf, sendcount, nlen, dest, sendtag)
    integer :: dest,sendtag,nlen
    integer :: sendcount
    character(len=nlen),dimension(sendcount):: sendbuf

    integer :: ier

    call MPI_SEND(sendbuf,sendcount*nlen,MPI_CHARACTER,dest,sendtag,my_local_mpi_comm_world,ier)

  end subroutine send_ch_array

  subroutine recv_ch_array(recvbuf, recvcount, nlen, dest, recvtag)
    integer :: dest,recvtag,nlen
    integer :: recvcount
    character(len=nlen),dimension(recvcount):: recvbuf

    integer :: ier

    call MPI_RECV(recvbuf,recvcount*nlen,MPI_CHARACTER,dest,recvtag, &
                  my_local_mpi_comm_world,MPI_STATUS_IGNORE,ier)

  end subroutine recv_ch_array
!

  subroutine land_all_all_l(sendbuf, recvbuf)
    logical, intent(in) :: sendbuf
    logical, intent(out) :: recvbuf

    call MPI_ALLREDUCE(sendbuf, recvbuf, 1, MPI_LOGICAL, MPI_LAND, my_local_mpi_comm_world, ier)
    
  end subroutine land_all_all_l

  subroutine prepare_shm_array_cr_1d(buffer, n_elem, win)
    USE, INTRINSIC :: ISO_C_BINDING
    real(kind=cr), dimension(:), pointer :: buffer
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

  end subroutine prepare_shm_array_cr_1d

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

  subroutine sync_from_main_rank_cr_1d(buffer, countval)
    integer, intent(in) :: countval
    real(kind=cr), dimension(:), intent(inout) :: buffer
    integer :: tag = 1000, i

    if (worldrank == 0) then
      do i = 2, worldsize
        if (rank_map(i, 2) == 0) then
          call send_r(buffer, countval, rank_map(i, 1), tag)
        endif
      enddo
    elseif (noderank == 0) then
      call recv_r(buffer, countval, 0, tag)
    endif
    call synchronize_all()
    
  end subroutine sync_from_main_rank_cr_1d

  subroutine sync_from_main_rank_dp_2d(buffer, nx, ny)
    integer, intent(in) :: nx, ny
    double precision, dimension(:,:), intent(inout) :: buffer
    integer :: tag = 1000, i

    if (worldrank == 0) then
      do i = 2, worldsize
        if (rank_map(i, 2) == 0) then
          call send_dp(buffer, nx*ny, rank_map(i, 1), tag)
        endif
      enddo
    elseif (noderank == 0) then
      call recv_dp(buffer, nx*ny, 0, tag)
    endif
    call synchronize_all()
    
  end subroutine sync_from_main_rank_dp_2d

  subroutine sync_from_main_rank_dp_3d(buffer, nx, ny, nz)
    integer, intent(in) :: nx, ny, nz
    double precision, dimension(:,:,:), intent(inout) :: buffer
    integer :: tag = 1000, i

    if (worldrank == 0) then
      do i = 2, worldsize
        if (rank_map(i, 2) == 0) then
          call send_dp(buffer, nx*ny*nz, rank_map(i, 1), tag)
        endif
      enddo
    elseif (noderank == 0) then
      call recv_dp(buffer, nx*ny*nz, 0, tag)
    endif
    call synchronize_all()
    
  end subroutine sync_from_main_rank_dp_3d


  subroutine sync_from_main_rank_ch(buffer, countval, nlen)
    integer, intent(in) :: countval, nlen
    character(len=nlen), dimension(:), intent(inout) :: buffer
    integer :: tag = 1000, i

    if (worldrank == 0) then
      do i = 2, worldsize
        if (rank_map(i, 2) == 0) then
          call send_ch_array(buffer, countval, nlen, rank_map(i, 1), tag)
        endif
      enddo
    elseif (noderank == 0) then
      call recv_ch_array(buffer, countval, nlen, 0, tag)
    endif
    call synchronize_all()
    
  end subroutine sync_from_main_rank_ch

end module fwat_mpi