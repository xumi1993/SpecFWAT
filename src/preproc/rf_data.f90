module rf_data
  use config
  use ma_constants
  use common_lib, only: get_band_name, rotate_R_to_NE_dp, dwascii
  use signal, only: bandpass_dp, interpolate_syn_dp, detrend, demean, &
                    myconvolution_dp, time_deconv
  use syn_data, only: SynData, average_amp_scale
  use obs_data, only: ObsData
  use input_params, fpar => fwat_par_global
  use decon_mod, only: deconit
  use fk_coupling
  use fwat_mpi
  use utils, only: zeros_dp, zeros, interp1
  use sacio
  use logger, only: log
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, nrec, nrec_local,OUTPUT_FILES, &
                        number_receiver_global, ispec_selected_rec, islice_selected_rec

  implicit none

  integer :: ier
  type, extends(SynData) :: RFData
    real(kind=cr), dimension(:), pointer :: ttp
    real(kind=dp), dimension(:, :, :), pointer :: rf_dat, rf_syn
    real(kind=cr) :: baz, az
    integer :: ttp_win, rf_win, syn_win
    contains
    procedure :: semd2sac, preprocess, finalize
    procedure, private :: calc_times, calc_rf, interp_data, measure_adj
  end type RFData

  contains

  subroutine semd2sac(this, ievt)
    class(RFData), intent(inout) :: this
    integer, intent(in) :: ievt
    real(kind=cr), dimension(:), allocatable :: bazi
    real(kind=dp), dimension(:), allocatable :: uin, win, rfi
    type(sachead) :: header
    integer :: irec_local, irec, igauss
    character(len=MAX_STRING_LEN) :: datafile, bandname

    call this%init(ievt)
    call this%od%read_stations(ievt, .true.)
    call this%calc_times()

    bazi = zeros(nrec)+this%baz
    call this%read(bazi)

    if (worldrank == 0) call system('mkdir -p '//trim(fpar%acqui%in_dat_path(this%ievt)))
    call synchronize_all()

    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        ! calculate rf
        uin = this%data(:, 2, irec)
        win = this%data(:, 1, irec)
        call bandpass_dp(uin ,NSTEP, dble(DT),1/fpar%sim%LONG_P(1),&
                         1/fpar%sim%SHORT_P(1), IORD)
        call bandpass_dp(win ,NSTEP, dble(DT),1/fpar%sim%LONG_P(1),&
                         1/fpar%sim%SHORT_P(1), IORD)
        do igauss = 1, fpar%sim%rf%NGAUSS
          write(bandname,'(a1,f3.1)') 'F',fpar%sim%rf%f0(igauss)
          datafile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%od%netwk(irec))//'.'&
                     //trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//'R.'&
                     //trim(bandname)//'.rf.sac'
          call deconit(uin, win, real(DT), fpar%sim%rf%tshift, fpar%sim%rf%f0(igauss),&
                       fpar%sim%rf%maxit, fpar%sim%rf%minderr, 0, rfi)
          call sacio_newhead(header, real(fpar%sim%dt), fpar%sim%nstep, -fpar%sim%rf%tshift)
          header%az = this%az
          header%baz = this%baz
          header%stla = this%od%stla(irec)
          header%stlo = this%od%stlo(irec)
          header%stel = this%od%stel(irec)
          header%knetwk = this%od%netwk(irec)
          header%kstnm = this%od%stnm(irec)
          header%kcmpnm = trim(fpar%sim%CH_CODE)//'R'
          header%user1 = fpar%sim%rf%f0(igauss)
          header%kuser1 = 'gauss'
          call sacio_writesac(datafile, header, rfi, ier)
        enddo
      enddo
    endif
    call synchronize_all()

  end subroutine semd2sac

  subroutine preprocess(this, ievt)
    class(RFData), intent(inout) :: this
    integer, intent(in) :: ievt    
    integer :: irec_local, irec
    
    fpar%sim%NRCOMP = 1
    fpar%sim%RCOMPS(1) = 'R'

    call this%init(ievt)

    call this%od%read_stations(ievt)

    call this%od%read_obs_data()
    call this%interp_data()

    call this%read(this%od%baz)
    call this%calc_rf()

    allocate(this%wchi(fpar%sim%rf%NGAUSS))

    ! measure adjoint source
    call this%measure_adj()

  end subroutine preprocess

  subroutine measure_adj(this)
    use measure_adj_mod, only: measure_adj_rf
    class(RFData), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: chan, msg
    integer :: irec_local, irec, igaus
    real(kind=dp), dimension(:,:), allocatable :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn
    real(kind=dp), dimension(:), allocatable :: tstart, tend, synz, synr, adj_2, adj_3
    real(kind=dp), dimension(:,:,:), allocatable :: window_chi, adj_src
    real(kind=dp), dimension(:), allocatable :: adj_r_tw, adj_z_tw, total_misfit
    character(len=MAX_STR_CHI), dimension(:), allocatable :: sta, net
    
    if (nrec_local > 0) then
      adj_src = zeros_dp(NSTEP, 2, nrec_local)
    endif
    total_misfit = zeros_dp(fpar%sim%rf%NGAUSS)
    chan = fpar%sim%CH_CODE//'R'
    do igaus = 1, fpar%sim%rf%NGAUSS
      write(this%band_name, '(a1,F3.1)') 'F', fpar%sim%rf%f0(igaus)
      call this%wchi(igaus)%init(this%ievt, this%band_name)
      if (nrec_local > 0) then
        window_chi = zeros_dp(nrec_local, NCHI, 1)
        tr_chi = zeros_dp(nrec_local, 1)
        am_chi = zeros_dp(nrec_local, 1)
        T_pmax_dat = zeros_dp(nrec_local, 1)
        T_pmax_syn = zeros_dp(nrec_local, 1)
        tstart = zeros_dp(nrec_local)
        tend = zeros_dp(nrec_local)
        allocate(sta(nrec_local))
        allocate(net(nrec_local))
        do irec_local = 1, nrec_local
          irec = number_receiver_global(irec_local)
          synz = this%data(:, 1, irec)
          synr = this%data(:, 2, irec)
          adj_z_tw = zeros_dp(NSTEP)
          adj_r_tw = zeros_dp(NSTEP)
          call bandpass_dp(synz, NSTEP, dble(DT), 1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          call bandpass_dp(synr, NSTEP, dble(DT), 1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          call measure_adj_rf(this%rf_dat(:, igaus, irec), this%rf_syn(:, igaus, irec),&
                              synr, synz, dble(fpar%sim%TIME_WIN(1)), dble(fpar%sim%TIME_WIN(2)),&
                              dble(-t0), dble(this%ttp(irec)), dble(DT), NSTEP, fpar%sim%rf%f0(igaus),&
                              fpar%sim%rf%tshift, fpar%sim%rf%maxit, fpar%sim%rf%minderr, &
                              window_chi(irec_local, :, 1), adj_r_tw, adj_z_tw)
          adj_src(:, 1, irec_local) = adj_src(:, 1, irec_local) + adj_z_tw
          adj_src(:, 2, irec_local) = adj_src(:, 2, irec_local) + adj_r_tw
          tr_chi(irec_local, 1) = fpar%acqui%src_weight(this%ievt)*window_chi(irec_local, 15, 1)
          am_chi(irec_local, 1) = fpar%acqui%src_weight(this%ievt)*window_chi(irec_local, 15, 1)
          T_pmax_dat(irec_local, 1) = 0.0_dp
          T_pmax_syn(irec_local, 1) = 0.0_dp
          if (igaus == 1) then
            sta(irec_local) = this%od%stnm(irec)
            net(irec_local) = this%od%netwk(irec)
            tstart(irec_local) = this%ttp(irec) - dble(fpar%sim%TIME_WIN(1)) 
            tend(irec_local) =  this%ttp(irec) + dble(fpar%sim%TIME_WIN(2))
          endif
        enddo
      endif
      call this%wchi(igaus)%assemble_window_chi(window_chi, tr_chi, am_chi,&
                                               T_pmax_dat, T_pmax_syn, sta, net,&
                                               tstart, tend)
      call this%wchi(igaus)%write()
      total_misfit(igaus) = this%wchi(igaus)%sum_chi(29)
      write(msg, '(a,f12.6)') 'Total misfit: of '//trim(this%band_name)//': ', total_misfit
      call log%write(msg, .true.)
    enddo
    call synchronize_all()

    adj_src = adj_src * fpar%acqui%src_weight(this%ievt)
    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        call this%write_adj(adj_src(:, 1, irec_local), this%comp_name(1), irec)
        adj_2 = zeros_dp(NSTEP)
        adj_3 = zeros_dp(NSTEP)
        call rotate_R_to_NE_dp(adj_src(:, 2, irec_local), adj_2, adj_3, this%baz)
        call this%write_adj(adj_2, this%comp_name(2), irec)
        call this%write_adj(adj_3, this%comp_name(3), irec)
      enddo
    endif
    call synchronize_all()
  end subroutine measure_adj

  subroutine interp_data(this)
    class(RFData), intent(inout) :: this
    integer :: irec, igaus

    call prepare_shm_array_dp_3d(this%rf_syn, NSTEP, fpar%sim%rf%NGAUSS, this%nrec, this%syn_win)
    if (worldrank == 0) then
      do igaus = 1, fpar%sim%rf%NGAUSS
        do irec = 1, this%nrec
          call interpolate_syn_dp(this%od%data(:, igaus, irec), dble(this%od%tbeg(irec)),&
                                  dble(this%od%dt), this%od%npts, &
                                  -dble(fpar%sim%rf%tshift), dble(fpar%sim%dt),NSTEP)
        enddo
      enddo
    endif
    call synchronize_all()
    call sync_from_main_rank_dp_3d(this%rf_syn, NSTEP, fpar%sim%rf%NGAUSS, this%nrec)

  end subroutine interp_data

  subroutine calc_rf(this)
    class(RFData), intent(inout) :: this
    real(kind=dp), dimension(:), allocatable :: uin, win, rfi
    real(kind=dp), dimension(:,:,:), allocatable :: rf_data_local, recv_buffer
    integer, dimension(:), allocatable :: recv_indices, send_indices
    integer :: igaus, irec_local, irec, iproc, nsta_irank, i

    call prepare_shm_array_dp_3d(this%rf_dat, NSTEP, fpar%sim%rf%NGAUSS, this%nrec, this%rf_win)

    if (nrec_local > 0) then
      rf_data_local = zeros_dp(NSTEP, fpar%sim%rf%NGAUSS, nrec_local)
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        uin = this%data(:, 2, irec)
        win = this%data(:, 1, irec)
        call bandpass_dp(uin ,NSTEP, dble(DT),1/fpar%sim%LONG_P(1),&
                         1/fpar%sim%SHORT_P(1), IORD)
        call bandpass_dp(win ,NSTEP, dble(DT),1/fpar%sim%LONG_P(1),&
                         1/fpar%sim%SHORT_P(1), IORD)
        do igaus = 1, fpar%sim%rf%NGAUSS
          call deconit(uin, win, real(DT), fpar%sim%rf%tshift, fpar%sim%rf%f0(igaus),&
                       fpar%sim%rf%maxit, fpar%sim%rf%minderr, 0, rfi)
          rf_data_local(:, igaus, irec_local) = rfi
        enddo
      enddo
    endif

    ! collect to main rank
    if (worldrank == 0) then
      if (nrec_local > 0) then
        do irec_local = 1, nrec_local
          irec = number_receiver_global(irec_local)
          this%rf_dat(:, :, irec) = rf_data_local(:, :, irec_local)
        enddo
      endif
      do iproc = 1, worldsize-1
        nsta_irank = 0
        do irec = 1, nrec
          if (ispec_selected_rec(irec) == iproc) nsta_irank = nsta_irank + 1
        enddo
        if (nsta_irank > 0) then
          allocate(recv_buffer(NSTEP, fpar%sim%rf%NGAUSS, nsta_irank))
          allocate(recv_indices(nsta_irank))
          call recv_i(recv_indices, nsta_irank, iproc, targ)
          call recv_dp(recv_buffer, NSTEP*fpar%sim%rf%NGAUSS*nsta_irank, iproc, targ)
          do i = 1, nsta_irank
            irec = recv_indices(i)
            this%rf_dat(:, :, irec) = recv_buffer(:, :, i)
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
        call send_dp(rf_data_local, NSTEP*fpar%sim%rf%NGAUSS*nrec_local, 0, targ)
        deallocate(send_indices)
      endif
    endif
    call sync_from_main_rank_dp_3d(this%rf_dat, NSTEP, fpar%sim%rf%NGAUSS, this%nrec)
  end subroutine calc_rf

  subroutine calc_times(this)
    class(RFData), intent(inout) :: this
    real(kind=cr), dimension(:), allocatable :: ttp_local
    integer :: irec_local, irec

    ttp_local = zeros(this%nrec)
    call prepare_shm_array_cr_1d(this%ttp, this%nrec, this%ttp_win)
    call read_fk_model(fpar%acqui%evtid_names(this%ievt))

    this%baz = -phi_FK - 90.d0
    this%az = 90.d0 - phi_FK
    
    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        ttp_local(irec) = maxloc(seismograms_d(3, irec_local, 1:NSTEP), dim=1) * DT - T0
      end do
    endif
    call sum_all_1Darray_cr(ttp_local, this%ttp, this%nrec)
    call sync_from_main_rank_cr_1d(this%ttp, this%nrec)
    call free_fk_arrays()

  end subroutine calc_times

  subroutine finalize(this)
    class(RFData), intent(inout) :: this
    integer :: igaus

    call this%od%finalize()
    do igaus = 1, fpar%sim%rf%NGAUSS
      call this%wchi(igaus)%finalize()
    enddo
    call free_shm_array(this%ttp_win)
    call free_shm_array(this%rf_win)
    call free_shm_array(this%syn_win)
  end subroutine finalize

end module rf_data