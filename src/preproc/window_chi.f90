module window_chi
  use config
  use input_params, only: fpar => fwat_par_global
  use fwat_constants
  use ma_constants
  use specfem_par, only: nrec, number_receiver_global, islice_selected_rec
  use fwat_mpi

  implicit none

  type, public :: WindowChi
    integer :: nrow, ncol = 12+NCHI, FID, ncomp
    real(kind=dp), dimension(:,:,:), pointer :: chi
    real(kind=dp), dimension(:,:), pointer :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn
    real(kind=dp), dimension(:), pointer :: tstart, tend
    character(len=MAX_STR_CHI), dimension(:), pointer :: sta,net
    character(len=MAX_STRING_LEN) :: evtid
    integer :: win_chi, win_tr, win_am, win_dat, win_syn, win_sta, win_net,&
                win_tstart, win_tend
    contains
    procedure :: init => init_window_chi, finalize, get_column, mean_chi, sum_chi, assemble_window_chi,&
                  write, read
  end type WindowChi

contains
  subroutine init_window_chi(this, ievt, band_name, model_name_in, action)
    class(WindowChi), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: chifile, status
    character(len=*), intent(in) :: band_name
    character(len=*), optional, intent(in) :: action, model_name_in
    character(len=MAX_STRING_LEN) :: model_name_loc
    integer :: ievt

    model_name_loc = trim(model_name)
    if (present(model_name_in)) then
      model_name_loc = trim(model_name_in)
    endif
    
    if (present(action)) then
      status = trim(action)
    else
      status = 'unknown'
    endif

    this%nrow = nrec
    this%ncomp = fpar%sim%NRCOMP
    this%evtid = fpar%acqui%evtid_names(ievt)
    call prepare_shm_array_dp_3d(this%chi, nrec, NCHI, this%ncomp, this%win_chi)
    call prepare_shm_array_dp_2d(this%tr_chi, nrec, this%ncomp, this%win_tr)
    call prepare_shm_array_dp_2d(this%am_chi, nrec, this%ncomp, this%win_am)
    call prepare_shm_array_dp_2d(this%T_pmax_dat, nrec, this%ncomp, this%win_dat)
    call prepare_shm_array_dp_2d(this%T_pmax_syn, nrec, this%ncomp, this%win_syn)
    call prepare_shm_array_ch_1d(this%sta, nrec, MAX_STR_CHI, this%win_sta)
    call prepare_shm_array_ch_1d(this%net, nrec, MAX_STR_CHI, this%win_net)
    call prepare_shm_array_dp_1d(this%tstart, nrec, this%win_tstart)
    call prepare_shm_array_dp_1d(this%tend, nrec, this%win_tend)
    chifile = trim(MISFITS_DIR)//'/'//trim(model_name_loc)//'.'//&
              trim(fpar%acqui%evtid_names(ievt))//'_'//trim(band_name)//&
              '_window_chi'
    if (worldrank == 0) open(newunit=this%FID, file=chifile, status=status)
    call synchronize_all()

  end subroutine init_window_chi

  function get_column(this, col, icomp) result(array)
    class(WindowChi), intent(inout) :: this
    integer, intent(in) :: col, icomp
    integer :: i
    real(kind=dp), dimension(:), allocatable :: array

    allocate(array(this%nrow))

    if (col <= 8 .or. col > this%ncol) then
      stop 'mean value will be calculated when column > 8'
    elseif (col <= 28) then
      i = col - 8
      array(:) = this%chi(:, i, icomp)
    elseif (col == 29) then
      array(:) = this%tr_chi(:, icomp)
    elseif (col == 30) then
      array(:) = this%am_chi(:, icomp)
    elseif (col == 31) then
      array(:) = this%T_pmax_dat(:, icomp)
    elseif (col == 32) then
      array(:) = this%T_pmax_syn(:, icomp)
    endif

  end function get_column

  function sum_chi(this, col) result(sum_value)
    class(WindowChi), intent(inout) :: this
    integer, intent(in) :: col
    real(kind=dp) :: sum_value
    real(kind=dp), dimension(:), allocatable :: array
    integer :: i

    sum_value = 0.0_dp
    do i = 1, this%ncomp
      array = this%get_column(col, i)
      sum_value = sum_value + sum(array)
    enddo

  end function

  function mean_chi(this, col) result(mean_value)
    class(WindowChi) :: this
    integer, intent(in) :: col
    real(kind=dp) :: mean_value
    real(kind=dp) :: sum_value

    sum_value = this%sum_chi(col)
    mean_value = sum_value/this%nrow/this%ncomp
  end function

  subroutine assemble_window_chi(this, window_chi_local, tr_chi_local, am_chi_local,&
                                 T_pmax_dat_local, T_pmax_syn_local, sta_local, net_local,&
                                 tstart_local, tend_local)
    class(WindowChi), intent(inout) :: this
    real(kind=dp), dimension(:,:,:), intent(in) :: window_chi_local
    real(kind=dp), dimension(:,:), intent(in) :: tr_chi_local, am_chi_local, T_pmax_dat_local, T_pmax_syn_local
    real(kind=dp), dimension(:), intent(in) :: tstart_local, tend_local
    character(len=MAX_STR_CHI), dimension(:), intent(in) :: sta_local, net_local
    integer :: irec, irec_local, nsta_irank, iproc, i, nrec_loc
    real(kind=dp), dimension(:,:,:), allocatable :: recv_win_chi
    real(kind=dp), dimension(:,:), allocatable :: recv_tr_chi, recv_am_chi, recv_T_pmax_dat, recv_T_pmax_syn
    character(len=MAX_STR_CHI), dimension(:), allocatable :: recv_sta, recv_net
    real(kind=dp), dimension(:), allocatable :: recv_tstart, recv_tend
    integer, dimension(:), allocatable :: recv_indices, send_indices
    
    nrec_loc = get_num_recs_per_proc(nrec, worldrank)

    if (worldrank == 0) then
      if (nrec_loc > 0) then
        do irec_local = 1, nrec_loc
          irec = select_global_id_for_rec(irec_local)
          this%chi(irec, :, :) = window_chi_local(irec_local, :, :)
          this%tr_chi(irec, :) = tr_chi_local(irec_local, :)
          this%am_chi(irec, :) = am_chi_local(irec_local, :)
          this%T_pmax_dat(irec, :) = T_pmax_dat_local(irec_local, :)
          this%T_pmax_syn(irec, :) = T_pmax_syn_local(irec_local, :)
          this%sta(irec) = sta_local(irec_local)
          this%net(irec) = net_local(irec_local)
          this%tstart(irec) = tstart_local(irec_local)
          this%tend(irec) = tend_local(irec_local)
        enddo
      endif
      do iproc = 1, worldsize-1
        nsta_irank = get_num_recs_per_proc(nrec, iproc)
        ! do irec = 1, nrec
        !   if (islice_selected_rec(irec) == iproc) nsta_irank = nsta_irank + 1
        ! enddo
        if (nsta_irank > 0) then
          allocate(recv_win_chi(nsta_irank, NCHI, this%ncomp))  ! Allocate a buffer to receive data
          allocate(recv_tr_chi(nsta_irank, this%ncomp))
          allocate(recv_am_chi(nsta_irank, this%ncomp))
          allocate(recv_T_pmax_dat(nsta_irank, this%ncomp))
          allocate(recv_T_pmax_syn(nsta_irank, this%ncomp))
          allocate(recv_sta(nsta_irank))
          allocate(recv_net(nsta_irank))
          allocate(recv_tstart(nsta_irank))
          allocate(recv_tend(nsta_irank))
          allocate(recv_indices(nsta_irank)) 
          ! Receive the indices first
          call recv_i(recv_indices, nsta_irank, iproc, targ)
          ! Receive the data
          call recv_dp(recv_win_chi, NCHI*nsta_irank*this%ncomp, iproc, targ)
          call recv_dp(recv_tr_chi, nsta_irank*this%ncomp, iproc, targ)
          call recv_dp(recv_am_chi, nsta_irank*this%ncomp, iproc, targ)
          call recv_dp(recv_T_pmax_dat, nsta_irank*this%ncomp, iproc, targ)
          call recv_dp(recv_T_pmax_syn, nsta_irank*this%ncomp, iproc, targ)
          call recv_ch_array(recv_sta, nsta_irank, MAX_STR_CHI, iproc, targ)
          call recv_ch_array(recv_net, nsta_irank, MAX_STR_CHI, iproc, targ)
          call recv_dp(recv_tstart, nsta_irank, iproc, targ)
          call recv_dp(recv_tend, nsta_irank, iproc, targ)
          ! Copy the received data to the correct location
          do i = 1, nsta_irank
            irec = recv_indices(i)
            this%chi(irec, :, :) = recv_win_chi(i, :, :)
            this%tr_chi(irec, :) = recv_tr_chi(i, :)
            this%am_chi(irec, :) = recv_am_chi(i, :)
            this%T_pmax_dat(irec, :) = recv_T_pmax_dat(i, :)
            this%T_pmax_syn(irec, :) = recv_T_pmax_syn(i, :)
            this%sta(irec) = recv_sta(i)
            this%net(irec) = recv_net(i)
            this%tstart(irec) = recv_tstart(i)
            this%tend(irec) = recv_tend(i)
          enddo
          deallocate(recv_win_chi)
          deallocate(recv_tr_chi)
          deallocate(recv_am_chi)
          deallocate(recv_T_pmax_dat)
          deallocate(recv_T_pmax_syn)
          deallocate(recv_indices)
          deallocate(recv_sta)
          deallocate(recv_net)
          deallocate(recv_tstart)
          deallocate(recv_tend)
        endif
      enddo
    else
      if (nrec_loc > 0) then
        allocate(send_indices(nrec_loc))
        do irec_local = 1, nrec_loc
          send_indices(irec_local) = select_global_id_for_rec(irec_local)
        enddo
        call send_i(send_indices, nrec_loc, 0, targ)
        call send_dp(window_chi_local, NCHI*nrec_loc*this%ncomp, 0, targ)
        call send_dp(tr_chi_local, nrec_loc*this%ncomp, 0, targ)
        call send_dp(am_chi_local, nrec_loc*this%ncomp, 0, targ)
        call send_dp(T_pmax_dat_local, nrec_loc*this%ncomp, 0, targ)
        call send_dp(T_pmax_syn_local, nrec_loc*this%ncomp, 0, targ)
        call send_ch_array(sta_local, nrec_loc, MAX_STR_CHI, 0, targ)
        call send_ch_array(net_local, nrec_loc, MAX_STR_CHI, 0, targ)
        call send_dp(tstart_local, nrec_loc, 0, targ)
        call send_dp(tend_local, nrec_loc, 0, targ)
        deallocate(send_indices)
      endif
    endif
    call synchronize_all()
    call sync_from_main_rank_dp_3d(this%chi, this%nrow, NCHI, this%ncomp)
    call sync_from_main_rank_dp_2d(this%tr_chi, this%nrow, this%ncomp)
    call sync_from_main_rank_dp_2d(this%am_chi, this%nrow, this%ncomp)
    call sync_from_main_rank_dp_2d(this%T_pmax_dat, this%nrow, this%ncomp)
    call sync_from_main_rank_dp_2d(this%T_pmax_syn, this%nrow, this%ncomp)
    call sync_from_main_rank_ch(this%sta, this%nrow, MAX_STR_CHI)
    call sync_from_main_rank_ch(this%net, this%nrow, MAX_STR_CHI)
    call sync_from_main_rank_dp_1d(this%tstart, this%nrow)
    call sync_from_main_rank_dp_1d(this%tend, this%nrow)

  end subroutine assemble_window_chi

  subroutine write(this)
    class(WindowChi), intent(inout) :: this
    integer :: irec, icomp
    
    if (worldrank == 0) then
      do irec = 1, this%nrow
        do icomp = 1, this%ncomp
          write(this%FID, '(a20,a8,a4,a5,i4,i4,2e14.6,20e14.6,2e14.6,2f14.6)') &
            this%evtid, this%sta(irec), this%net(irec), &
            trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp)), irec, & 
            fpar%sim%IMEAS, this%tstart(irec), this%tend(irec), &
            this%chi(irec, :, icomp), this%tr_chi(irec, icomp), this%am_chi(irec, icomp), &
            this%T_pmax_dat(irec, icomp), this%T_pmax_syn(irec, icomp)
        enddo
      enddo
    endif

  end subroutine write

  subroutine read(this, model_name_in, ievt, band_name)
    class(WindowChi), intent(inout) :: this
    character(len=*), intent(in) :: band_name, model_name_in
    character(len=MAX_STRING_LEN) :: chan
    integer, intent(in) :: ievt
    integer :: icomp, irec, irec_out, imeas_out

    call this%init(ievt, band_name, model_name_in=model_name_in, action='old')

    if (worldrank == 0) then
      do irec = 1, this%nrow
        do icomp = 1, this%ncomp
          read(this%FID, '(a20,a8,a4,a5,i4,i4,2e14.6,20e14.6,2e14.6,2f14.6)') &
            this%evtid, this%sta(irec), this%net(irec), &
            chan, irec_out, imeas_out, & 
            this%tstart(irec), this%tend(irec), &
            this%chi(irec, :, icomp), this%tr_chi(irec, icomp), this%am_chi(irec, icomp), &
            this%T_pmax_dat(irec, icomp), this%T_pmax_syn(irec, icomp)
        enddo
      enddo
      close(this%FID)
    endif
    call synchronize_all()

  end subroutine read

  subroutine finalize(this)
    class(WindowChi), intent(inout) :: this
    integer :: ier
    
    if (worldrank == 0) close(this%FID, iostat=ier)
    call free_shm_array(this%win_chi)
    call free_shm_array(this%win_tr)
    call free_shm_array(this%win_am)
    call free_shm_array(this%win_dat)
    call free_shm_array(this%win_syn)
    call free_shm_array(this%win_sta)
    call free_shm_array(this%win_net)
    call free_shm_array(this%win_tstart)
    call free_shm_array(this%win_tend)
    
  end subroutine finalize
end module window_chi 