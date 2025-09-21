module rf_data
  use config
  ! use ma_constants
  use rf_misfit
  use misfit_mod
  use common_lib, only: get_band_name, rotate_R_to_NE_dp, dwascii, mkdir
  use signal, only: bandpass_dp, interpolate_syn_dp, detrend, demean, &
                    myconvolution_dp, time_deconv
  use syn_data, only: SynData, average_amp_scale
  use obs_data, only: ObsData
  use input_params, fpar => fwat_par_global
  use decon_mod, only: deconit
  use fk_coupling
  use fwat_mpi
  use utils, only: zeros_dp, zeros
  use sacio
  use logger, only: log
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, OUTPUT_FILES

  implicit none

  integer :: ier
  type, extends(SynData) :: RFData
    real(kind=cr), dimension(:), pointer :: ttp
    real(kind=dp), dimension(:, :, :), pointer :: rf_dat, rf_syn
    real(kind=cr) :: baz, az
    integer :: ttp_win, rf_win, syn_win
    contains
    procedure :: semd2sac, preprocess, finalize
    procedure, private :: calc_times, calc_rf, interp_data, measure_adj, get_baz
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

    call this%get_baz()
    bazi = zeros(this%nrec)+this%baz
    call this%read(bazi)

    call this%calc_times()

    ! if (worldrank == 0) call system('mkdir -p '//trim(fpar%acqui%in_dat_path(this%ievt)))
    call mkdir(fpar%acqui%in_dat_path(this%ievt))
    call synchronize_all()

    if (this%nrec_loc > 0) then
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
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
          header%knetwk = trim(this%od%netwk(irec))
          header%kstnm = trim(this%od%stnm(irec))
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
    real(kind=cr), dimension(:), allocatable :: bazi
    
    fpar%sim%NRCOMP = 1
    fpar%sim%RCOMPS(1) = 'R'

    call this%init(ievt)

    call this%od%read_stations(ievt)

    call this%od%read_obs_data()
    call this%interp_data()

    call this%get_baz()
    bazi = zeros(this%nrec)+this%baz
    call this%read(bazi)
    call this%calc_times()

    call this%calc_rf()

    ! measure adjoint source
    call this%measure_adj()

  end subroutine preprocess

  subroutine measure_adj(this)
    ! use measure_adj_mod, only: measure_adj_rf
    class(RFData), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: msg
    integer :: irec_local, irec, igaus
    real(kind=dp), dimension(:,:), allocatable :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn
    ! real(kind=dp), dimension(:), allocatable :: tstart, tend, synz, synr, adj_2, adj_3
    real(kind=dp), dimension(:), allocatable :: synz, synr, adj_2, adj_3
    real(kind=dp), dimension(:,:,:), allocatable :: window_chi, adj_src
    ! real(kind=dp), dimension(:), allocatable :: adj_r_tw, adj_z_tw
    ! character(len=MAX_STR_CHI), dimension(:), allocatable :: sta, net
    type(RFMisfit) :: rfm
    type(RECMisfit), dimension(:), allocatable :: recm
    type(EVTMisfit) :: evtm

    if (this%nrec_loc > 0) then
      adj_src = zeros_dp(NSTEP, 2, this%nrec_loc)
      allocate(recm(this%nrec_loc))
      do irec_local = 1, this%nrec_loc
        recm(irec_local) = RECMisfit(fpar%sim%NRCOMP, 1, fpar%sim%rf%NGAUSS)
      enddo
    endif
    do igaus = 1, fpar%sim%rf%NGAUSS
      write(this%band_name, '(a1,F3.1)') 'F', fpar%sim%rf%f0(igaus)
      call this%wchi(igaus)%init(this%ievt, this%band_name)
      if (this%nrec_loc > 0) then
        window_chi = zeros_dp(this%nrec_loc, NCHI, 1)
        tr_chi = zeros_dp(this%nrec_loc, 1)
        am_chi = zeros_dp(this%nrec_loc, 1)
        T_pmax_dat = zeros_dp(this%nrec_loc, 1)
        T_pmax_syn = zeros_dp(this%nrec_loc, 1)
        do irec_local = 1, this%nrec_loc
          irec = select_global_id_for_rec(irec_local)
          synz = this%data(:, 1, irec)
          synr = this%data(:, 2, irec)
          ! adj_z_tw = zeros_dp(NSTEP)
          ! adj_r_tw = zeros_dp(NSTEP)
          call bandpass_dp(synz, NSTEP, dble(DT), 1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          call bandpass_dp(synr, NSTEP, dble(DT), 1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          call rfm%calc_adjoint_source(this%rf_dat(:, igaus, irec), this%rf_syn(:, igaus, irec), synr, synz, DT, &
                                       dble(-this%ttp(irec)+T0), dble([ fpar%sim%TIME_WIN(1), fpar%sim%TIME_WIN(2) ]), &
                                       dble(fpar%sim%rf%f0(igaus)), dble(fpar%sim%rf%tshift), &
                                       fpar%sim%rf%maxit, dble(fpar%sim%rf%minderr))
          ! call measure_adj_rf(this%rf_dat(:, igaus, irec), this%rf_syn(:, igaus, irec),&
          !                     synr, synz, -dble(fpar%sim%TIME_WIN(1)), dble(fpar%sim%TIME_WIN(2)),&
          !                     dble(-t0), dble(this%ttp(irec)), dble(DT), NSTEP, fpar%sim%rf%f0(igaus),&
          !                     fpar%sim%rf%tshift, fpar%sim%rf%maxit, fpar%sim%rf%minderr, &
          !                     window_chi(irec_local, :, 1), adj_r_tw, adj_z_tw,&
          !                     this%od%netwk(irec), this%od%stnm(irec) &
          !                     )
          adj_src(:, 1, irec_local) = adj_src(:, 1, irec_local) + rfm%adj_src_r 
          adj_src(:, 2, irec_local) = adj_src(:, 2, irec_local) + rfm%adj_src_z
          recm(irec_local)%misfits(1, 1, igaus) = rfm%misfits(1)
          recm(irec_local)%residuals(1, 1, igaus) = rfm%residuals(1)
          recm(irec_local)%imeas(1, 1, igaus) = rfm%imeas(1)
          recm(irec_local)%total_misfit(igaus) = recm(irec_local)%total_misfit(igaus) + rfm%total_misfit
          recm(irec_local)%band_name(igaus) = trim(this%band_name)
          if (igaus == 1) then
            recm(irec_local)%sta = trim(this%od%stnm(irec))
            recm(irec_local)%net = trim(this%od%netwk(irec))
            recm(irec_local)%chan(1) = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(1))
            recm(irec_local)%tstart(1, igaus) = this%ttp(irec) + dble(fpar%sim%TIME_WIN(1)) 
            recm(irec_local)%tend(1, igaus) = this%ttp(irec) + dble(fpar%sim%TIME_WIN(2))
          endif
        enddo
      endif
      ! call this%wchi(igaus)%assemble_window_chi(window_chi, tr_chi, am_chi,&
      !                                          T_pmax_dat, T_pmax_syn, sta, net,&
      !                                          tstart, tend)
      ! call this%wchi(igaus)%write()
      ! this%total_misfit(igaus) = this%wchi(igaus)%sum_chi(29)
      ! write(msg, '(a,f12.6)') 'Total misfit: of '//trim(this%band_name)//': ', this%total_misfit(igaus)
      ! call log%write(msg, .true.)
    enddo
    call synchronize_all()

    if (worldrank == 0) evtm = EVTMisfit(fpar%acqui%evtid_names(this%ievt), nrec, fpar%sim%rf%NGAUSS)
    call evtm%assemble(recm)
    if (worldrank == 0) then
      do igaus = 1, fpar%sim%rf%NGAUSS
        write(msg, '(a,f12.6)') 'Total misfit: of '//trim(this%band_name)//': ', evtm%total_misfit(igaus)
        call log%write(msg, .true.)
      enddo
    end if
    call evtm%write()

    call synchronize_all()

    adj_src = adj_src * fpar%acqui%src_weight(this%ievt)
    ! if (IS_OUTPUT_ADJ_SRC) then
    !   call sacio_newhead(header, real(DT), NSTEP, -real(T0))
    !   if (this%nrec_loc > 0) then
    !     do irec_local = 1, this%nrec_loc
    !       irec = select_global_id_for_rec(irec_local)
    !       header%kstnm = this%od%stnm(irec)
    !       header%knetwk = this%od%netwk(irec)
    !       do icomp = 1, 2
    !         if (icomp == 1) then
    !           header%kcmpnm = 'Z'
    !         else
    !           header%kcmpnm = 'R'
    !         end if
    !         call sacio_writesac(trim(fpar%acqui%out_fwd_path(this%ievt))//'/'//trim(ADJOINT_PATH)//'/'//&
    !                             trim(this%od%netwk(irec))//'.'//trim(this%od%stnm(irec))//&
    !                             '.'//trim(fpar%sim%CH_CODE)//trim(header%kcmpnm)//'.adj.sac', &
    !                             header, adj_src(:, icomp, irec_local), ier)
    !       enddo
    !     enddo
    !   endif
    ! endif
    call synchronize_all()

    call this%get_comp_name_adj()
    if (this%nrec_loc > 0) then
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
        call this%write_adj(adj_src(:, 1, irec_local), this%comp_name(1), irec)
        adj_2 = zeros_dp(NSTEP)
        adj_3 = zeros_dp(NSTEP)
        call rotate_R_to_NE_dp(adj_src(:, 2, irec_local), adj_3, adj_2, this%baz)
        ! write N component
        call this%write_adj(adj_2, this%comp_name(2), irec)
        ! write E component
        call this%write_adj(adj_3, this%comp_name(3), irec)
      enddo
    endif
    call synchronize_all()
  end subroutine measure_adj

  subroutine interp_data(this)
    class(RFData), intent(inout) :: this
    integer :: irec, igaus

    call prepare_shm_array_dp_3d(this%rf_dat, NSTEP, fpar%sim%rf%NGAUSS, this%nrec, this%syn_win)
    if (noderank == 0) then
      do igaus = 1, fpar%sim%rf%NGAUSS
        do irec = 1, this%nrec
          call interpolate_syn_dp(this%od%data(:, igaus, irec), dble(this%od%tbeg(irec)),&
                                  dble(this%od%dt), this%od%npts, &
                                  -dble(fpar%sim%rf%tshift), dble(fpar%sim%dt),NSTEP)
          this%rf_dat(:, igaus, irec) = this%od%data(:, igaus, irec)
        enddo
      enddo
    endif
    call synchronize_all()

  end subroutine interp_data

  subroutine calc_rf(this)
    class(RFData), intent(inout) :: this
    real(kind=dp), dimension(:), allocatable :: uin, win, rfi
    real(kind=dp), dimension(:,:,:), allocatable :: rf_data_local, recv_buffer
    integer, dimension(:), allocatable :: recv_indices, send_indices
    integer :: igaus, irec_local, irec, iproc, nsta_irank, i
    type(sachead) :: header

    call prepare_shm_array_dp_3d(this%rf_syn, NSTEP, fpar%sim%rf%NGAUSS, this%nrec, this%rf_win)

    if (this%nrec_loc > 0) then
      rf_data_local = zeros_dp(NSTEP, fpar%sim%rf%NGAUSS, this%nrec_loc)
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
        uin = this%data(:, 2, irec)
        win = this%data(:, 1, irec)
        call bandpass_dp(uin ,NSTEP, dble(DT),1/fpar%sim%LONG_P(1),&
                         1/fpar%sim%SHORT_P(1), IORD)
        call bandpass_dp(win ,NSTEP, dble(DT),1/fpar%sim%LONG_P(1),&
                         1/fpar%sim%SHORT_P(1), IORD)
        if (IS_OUTPUT_PREPROC) then
          call sacio_newhead(header, real(DT), NSTEP, -real(T0))
          header%az = this%az
          header%baz = this%baz
          header%kstnm = trim(this%od%stnm(irec))
          header%knetwk = trim(this%od%netwk(irec))
          header%kcmpnm = trim(fpar%sim%CH_CODE)//'Z'
          header%kevnm = trim(fpar%acqui%evtid_names(this%ievt))
          header%t0 = this%ttp(irec)
          call sacio_writesac(trim(fpar%acqui%out_fwd_path(this%ievt))//'/'//trim(OUTPUT_PATH)//'/wsyn.'//&
                              trim(this%od%netwk(irec))//'.'//trim(this%od%stnm(irec))//&
                              '.'//trim(fpar%sim%CH_CODE)//'Z.sac', header, win, ier)
          header%kcmpnm = trim(fpar%sim%CH_CODE)//'R'
          call sacio_writesac(trim(fpar%acqui%out_fwd_path(this%ievt))//'/'//trim(OUTPUT_PATH)//'/wsyn.'//&
                              trim(this%od%netwk(irec))//'.'//trim(this%od%stnm(irec))//&
                              '.'//trim(fpar%sim%CH_CODE)//'R.sac', header, uin, ier)
        endif
        do igaus = 1, fpar%sim%rf%NGAUSS
          call deconit(uin, win, real(DT), fpar%sim%rf%tshift, fpar%sim%rf%f0(igaus),&
                       fpar%sim%rf%maxit, fpar%sim%rf%minderr, 0, rfi)
          rf_data_local(:, igaus, irec_local) = rfi
          if (IS_OUTPUT_PREPROC) then
            call sacio_newhead(header, real(DT), NSTEP, -real(fpar%sim%rf%tshift))
            write(this%band_name, '(a1,F3.1)') 'F', fpar%sim%rf%f0(igaus)
            header%az = this%az
            header%baz = this%baz
            header%kstnm = trim(this%od%stnm(irec))
            header%knetwk = trim(this%od%netwk(irec))
            header%kcmpnm = trim(fpar%sim%CH_CODE)//'R'
            header%user1 = fpar%sim%rf%f0(igaus)
            header%kuser1 = 'gauss'
            call sacio_writesac(trim(fpar%acqui%out_fwd_path(this%ievt))//'/'//trim(OUTPUT_PATH)//'/syn.'//&
                                trim(this%od%netwk(irec))//'.'//trim(this%od%stnm(irec))//&
                                '.'//trim(fpar%sim%CH_CODE)//'R.rf.sac.'//trim(this%band_name), header, rfi, ier)
            call sacio_writesac(trim(fpar%acqui%out_fwd_path(this%ievt))//'/'//trim(OUTPUT_PATH)//'/obs.'//&
                                trim(this%od%netwk(irec))//'.'//trim(this%od%stnm(irec))//&
                                '.'//trim(fpar%sim%CH_CODE)//'R.rf.sac.'//trim(this%band_name), &
                                header, this%rf_dat(:, igaus, irec), ier)
          endif
        enddo
      enddo
    endif
    call synchronize_all()

    ! collect to main rank
    if (worldrank == 0) then
      if (this%nrec_loc > 0) then
        do irec_local = 1, this%nrec_loc
          irec = select_global_id_for_rec(irec_local)
          this%rf_syn(:, :, irec) = rf_data_local(:, :, irec_local)
        enddo
      endif
      do iproc = 1, worldsize-1
        nsta_irank = get_num_recs_per_proc(this%nrec, iproc)
        if (nsta_irank > 0) then
          allocate(recv_buffer(NSTEP, fpar%sim%rf%NGAUSS, nsta_irank))
          allocate(recv_indices(nsta_irank))
          call recv_i(recv_indices, nsta_irank, iproc, targ)
          call recv_dp(recv_buffer, NSTEP*fpar%sim%rf%NGAUSS*nsta_irank, iproc, targ)
          do i = 1, nsta_irank
            irec = recv_indices(i)
            this%rf_syn(:, :, irec) = recv_buffer(:, :, i)
          enddo
          deallocate(recv_buffer)
          deallocate(recv_indices)
        endif
      enddo
    else
      if (this%nrec_loc > 0) then
        allocate(send_indices(this%nrec_loc))
        do irec_local = 1, this%nrec_loc
          send_indices(irec_local) = select_global_id_for_rec(irec_local)
        enddo
        call send_i(send_indices, this%nrec_loc, 0, targ)
        call send_dp(rf_data_local, NSTEP*fpar%sim%rf%NGAUSS*this%nrec_loc, 0, targ)
        deallocate(send_indices)
      endif
    endif
    call synchronize_all()
    call sync_from_main_rank_dp_3d(this%rf_syn, NSTEP, fpar%sim%rf%NGAUSS, this%nrec)
  end subroutine calc_rf

  subroutine get_baz(this)
    class(RFData), intent(inout) :: this

    call read_fk_model(fpar%acqui%evtid_names(this%ievt)) 
    this%baz = -phi_FK - 90.d0
    this%az = 90.d0 - phi_FK
    call free_fk_arrays()
  end subroutine get_baz

  subroutine calc_times(this)
    class(RFData), intent(inout) :: this
    real(kind=cr), dimension(:), allocatable :: ttp_local
    integer :: irec_local, irec

    ttp_local = zeros(this%nrec)
    call prepare_shm_array_cr_1d(this%ttp, this%nrec, this%ttp_win)
    
    if (this%nrec_loc > 0) then
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
        ttp_local(irec) = maxloc(this%data(:, 1, irec), dim=1) * DT - T0
      end do
    endif
    call sum_all_1Darray_cr(ttp_local, this%ttp, this%nrec)
    call sync_from_main_rank_cr_1d(this%ttp, this%nrec)
  end subroutine calc_times

  subroutine finalize(this)
    class(RFData), intent(inout) :: this
    integer :: igaus

    call this%od%finalize()
    do igaus = 1, fpar%sim%rf%NGAUSS
      call this%wchi(igaus)%finalize()
    enddo
    call free_shm_array(this%dat_win)
    call free_shm_array(this%ttp_win)
    call free_shm_array(this%rf_win)
    call free_shm_array(this%syn_win)
  end subroutine finalize

end module rf_data