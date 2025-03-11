module window_chi
  use config
  use input_params, only: fpar => fwat_par_global
  use fwat_constants
  use ma_constants
  use specfem_par, only: nrec, nrec_local, number_receiver_global, islice_selected_rec
  use fwat_mpi

  implicit none

  type, public :: WindowChi
    integer :: nrow, ncol = 12+NCHI
    real(kind=dp), dimension(:,:), pointer :: chi
    real(kind=dp), dimension(:), pointer :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn
    character(len=MAX_STR_CHI), dimension(:), pointer :: evtid,sta,net,chan
    integer :: win_chi, win_tr, win_am, win_dat, win_syn, win_evtid, win_sta, win_net, win_chan
    contains
    procedure :: init => init_window_chi, finalize, get_column, mean, assemble_window_chi
  end type WindowChi

contains
  subroutine init_window_chi(this)
    class(WindowChi), intent(inout) :: this
    
    this%nrow = nrec
    call prepare_shm_array_dp_2d(this%chi, nrec, NCHI, this%win_chi)
    call prepare_shm_array_dp_1d(this%tr_chi, nrec, this%win_tr)
    call prepare_shm_array_dp_1d(this%am_chi, nrec, this%win_am)
    call prepare_shm_array_dp_1d(this%T_pmax_dat, nrec, this%win_dat)
    call prepare_shm_array_dp_1d(this%T_pmax_syn, nrec, this%win_syn)
    call prepare_shm_array_ch_1d(this%evtid, nrec, MAX_STR_CHI, this%win_evtid)
    call prepare_shm_array_ch_1d(this%sta, nrec, MAX_STR_CHI, this%win_sta)
    call prepare_shm_array_ch_1d(this%net, nrec, MAX_STR_CHI, this%win_net)
    call prepare_shm_array_ch_1d(this%chan, nrec, MAX_STR_CHI, this%win_chan)
  
  end subroutine init_window_chi

  function get_column(this, col) result(array)
    class(WindowChi), intent(inout) :: this
    integer :: col, i
    real(kind=dp), dimension(:), allocatable :: array

    allocate(array(this%nrow))

    if (col <= 8 .or. col > this%ncol) then
      stop 'mean value will be calculated when column > 8'
    elseif (col <= 28) then
      i = col - 8
      array(:) = this%chi(:, i)
    elseif (col == 29) then
      array(:) = this%tr_chi(:)
    elseif (col == 30) then
      array(:) = this%am_chi(:)
    elseif (col == 31) then
      array(:) = this%T_pmax_dat(:)
    elseif (col == 32) then
      array(:) = this%T_pmax_syn(:)
    endif

  end function get_column

  function mean(this, col) result(mean_value)
    class(WindowChi) :: this
    integer :: col
    real(kind=dp) :: mean_value
    real(kind=dp), dimension(:), allocatable :: array

    array = this%get_column(col)
    mean_value = sum(array)/this%nrow
  end function

  subroutine assemble_window_chi(this, window_chi_local)
    class(WindowChi), intent(inout) :: this
    real(kind=dp), dimension(:,:), intent(in) :: window_chi_local
    integer :: irec, irec_local, nsta_irank, iproc, i
    real(kind=dp), dimension(:,:), allocatable :: recv_buffer
    integer, dimension(:), allocatable :: recv_indices, send_indices

    if (worldrank == 0) then
      if (nrec_local > 0) then
        do irec_local = 1, nrec_local
          irec = number_receiver_global(irec_local)
          this%chi(irec, :) = window_chi_local(irec_local, :)
        enddo
      endif
      do iproc = 1, worldsize-1
        nsta_irank=0
        do irec = 1, nrec
          if (islice_selected_rec(irec) == iproc) nsta_irank = nsta_irank + 1
        enddo
      
        if (nsta_irank > 0) then
          allocate(recv_buffer(nsta_irank, NCHI))  ! Allocate a buffer to receive data
          allocate(recv_indices(nsta_irank)) 
          ! Receive the indices first
          call recv_i(recv_indices, nsta_irank, iproc, targ)
          ! Receive the data
          call recv_dp(recv_buffer, NCHI*nsta_irank, iproc, targ)
          ! Copy the received data to the correct location
          do i = 1, nsta_irank
            irec = recv_indices(i)
            this%chi(irec, :) = recv_buffer(i, :)
          enddo
          deallocate(recv_buffer)
          deallocate(recv_indices)
        endif
      enddo
    else
      if (nrec_local > 0) then
        allocate(send_indices(nrec_local))
        do irec_local = 1, nrec_local
          send_indices(irec_local) = number_receiver_global(irec_local)
        enddo
        call send_i(send_indices, nrec_local, 0, targ)
        call send_dp(window_chi_local, NCHI*nrec_local, 0, targ)
        deallocate(send_indices)
      endif
    endif
    call sync_from_main_rank_dp_2d(this%chi, this%nrow, NCHI)


  end subroutine assemble_window_chi

  subroutine finalize(this)
    class(WindowChi), intent(inout) :: this

    call free_shm_array(this%win_chi)
    call free_shm_array(this%win_tr)
    call free_shm_array(this%win_am)
    call free_shm_array(this%win_dat)
    call free_shm_array(this%win_syn)
    call free_shm_array(this%win_evtid)
    call free_shm_array(this%win_sta)
    call free_shm_array(this%win_net)
    call free_shm_array(this%win_chan)
    
  end subroutine finalize
end module window_chi 