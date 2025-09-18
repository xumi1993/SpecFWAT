module telecc_data
  use config
  use measure_adj_mod, only: measure_adj_cross_conv
  use ma_constants
  use common_lib, only: get_band_name, rotate_R_to_NE_dp, dwascii, &
                        get_gauss_fac
  use signal, only: bandpass_dp, interpolate_syn_dp, detrend, demean, &
                    myconvolution_dp, interpolate_func_dp, mycorrelation_dp
  use syn_data, only: SynData, average_amp_scale
  use obs_data, only: ObsData
  use input_params, fpar => fwat_par_global
  use fk_coupling
  use fwat_mpi
  use utils, only: zeros_dp, zeros, interp1, find_maxima_dp
  use sacio
  use logger, only: log
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, DT, NSTEP, nrec,OUTPUT_FILES
  use tele_data, only: TeleData

  implicit none

  integer, private :: ier, ncomp

  type, extends(TeleData) :: TeleCCData
  contains
    procedure :: preprocess, finalize
    procedure, private :: measure_adj, pre_proc, calc_times
  end type TeleCCData

contains

  subroutine preprocess(this, ievt)
    class(TeleCCData), intent(inout) :: this
    integer, intent(in) :: ievt
    character(len=MAX_STRING_LEN) :: msg

    fpar%sim%NRCOMP = 1
    ncomp = size(fpar%sim%RCOMPS)

    call this%init(ievt)

     ! read receiver information
    call this%od%read_stations(this%ievt)

    ! read observed data
    call this%od%read_obs_data()

    ! read syn data
    call this%read(this%od%baz)

    call get_band_name(fpar%sim%SHORT_P(1), fpar%sim%LONG_P(1), this%band_name)

    ! initialize misfits
    call this%wchi(1)%init(ievt, this%band_name)

    call this%get_comp_name_adj()

    ! call this%calc_fktimes()
    call this%calc_times()

    call synchronize_all()

    call this%pre_proc()

    call this%measure_adj()

    call this%wchi(1)%assemble_window_chi(this%window_chi, this%tr_chi, this%am_chi,&
                                          this%T_pmax_dat, this%T_pmax_syn, this%sta, this%net,&
                                          this%tstart, this%tend)

    if (worldrank == 0) then
      call this%wchi(1)%write()
      this%total_misfit(1) = this%wchi(1)%sum_chi(29)
      write(msg, '(a,f12.6)') 'Total misfit: of '//trim(this%band_name)//': ', this%total_misfit(1)
      call log%write(msg, .true.)
    endif
    call synchronize_all()

  end subroutine preprocess

  subroutine pre_proc(this)
    class(TeleCCData), intent(inout) :: this
    integer :: irec_local, irec, icomp
    real(kind=dp), dimension(:), allocatable :: seismo_inp
    real(kind=dp) :: tstart, tend, t01, max_amp

    if (worldrank == 0) then
      max_amp = average_amp_scale(this%od%data, 1)
    endif
    call synchronize_all()
    call bcast_all_singledp(max_amp)

    if (this%nrec_loc > 0) then
      this%seismo_dat = zeros_dp(NSTEP, ncomp, this%nrec_loc)
      this%seismo_syn = zeros_dp(NSTEP, ncomp, this%nrec_loc)
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
        tstart = this%ttp(irec) + fpar%sim%time_win(1)
        tend = this%ttp(irec) + fpar%sim%time_win(2)
        if (tend > NSTEP*DT-T0) then
          call log%write('Error: Time windows length of '//trim(this%od%netwk(irec))&
                          //'.'//trim(this%od%stnm(irec))//' larger than NSTEP', .false.)
          call exit_mpi(worldrank, 'Error: Time windows length larger than NSTEP')
        endif
        do icomp = 1, ncomp
          t01 = this%od%tbeg(irec) - (this%od%tarr(irec) - this%ttp(irec))
          seismo_inp = interpolate_func_dp(this%od%data(:, icomp, irec), t01,&
                                           dble(this%od%dt), this%od%npts, &
                                           -dble(T0), dble(DT), NSTEP)/max_amp
          call detrend(seismo_inp)
          call demean(seismo_inp)
          call bandpass_dp(seismo_inp, NSTEP, dble(DT),&
                           1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          ! call interpolate_syn_dp(seismo_inp, -dble(T0), dble(DT), NSTEP, &
          !                         tstart, dble(dt), NSTEP)
          this%seismo_dat(:, icomp, irec_local) = seismo_inp(1:NSTEP)

          seismo_inp = this%data(:, icomp, irec)
          ! call detrend(seismo_inp)
          ! call demean(seismo_inp)
          call bandpass_dp(seismo_inp, NSTEP, dble(DT),&
                           1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          ! call interpolate_syn_dp(seismo_inp, -dble(T0), dble(DT), NSTEP, &
          !                         tstart, dble(dt), NSTEP)
          this%seismo_syn(:, icomp, irec_local) = seismo_inp(1:NSTEP)
        enddo
      enddo
    endif
    call synchronize_all()
  
  end subroutine pre_proc

  subroutine measure_adj(this)
    class(TeleCCData), intent(inout) :: this
    integer :: irec_local, irec
    real(kind=dp), dimension(:), allocatable :: adj_r_tw, adj_z_tw
    real(kind=dp), dimension(:,:), allocatable :: adj_src
    real(kind=dp), dimension(NCHI) :: window_chi_local
    type(sachead) :: header


    if (this%nrec_loc > 0) then
      allocate(this%window_chi(this%nrec_loc, NCHI, fpar%sim%NRCOMP))
      allocate(this%tr_chi(this%nrec_loc, fpar%sim%NRCOMP))
      allocate(this%am_chi(this%nrec_loc, fpar%sim%NRCOMP))
      allocate(this%T_pmax_dat(this%nrec_loc, fpar%sim%NRCOMP))
      allocate(this%T_pmax_syn(this%nrec_loc, fpar%sim%NRCOMP))
      allocate(this%sta(this%nrec_loc))
      allocate(this%net(this%nrec_loc))
      allocate(this%tstart(this%nrec_loc))
      allocate(this%tend(this%nrec_loc))
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
        if (IS_OUTPUT_PREPROC) then
          block
            character(len=MAX_STRING_LEN) :: sacfile
            call sacio_newhead(header, real(DT), NSTEP, -real(T0))
            header%knetwk = trim(this%od%netwk(irec))
            header%kstnm = trim(this%od%stnm(irec))
            header%kcmpnm = trim(fpar%sim%CH_CODE)//'Z'
            header%kevnm = trim(fpar%acqui%evtid_names(this%ievt))
            header%t0 = this%ttp(irec)
            header%baz = this%od%baz(irec)
            header%dist = this%od%dist(irec)
            sacfile = trim(OUTPUT_FILES)//'/dat.'//trim(this%od%netwk(irec))//&
                      '.'//trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
                      'Z.sac.'//trim(this%band_name)
            call sacio_writesac(sacfile, header, this%seismo_dat(:, 1, irec_local), ier)

            sacfile = trim(OUTPUT_FILES)//'/syn.'//trim(this%od%netwk(irec))//&
                      '.'//trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
                      'Z.sac.'//trim(this%band_name)
            call sacio_writesac(sacfile, header, this%seismo_syn(:, 1, irec_local), ier)
            
            header%kcmpnm = trim(fpar%sim%CH_CODE)//'R'
            sacfile = trim(OUTPUT_FILES)//'/dat.'//trim(this%od%netwk(irec))//&
                      '.'//trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
                      'R.sac.'//trim(this%band_name)
            call sacio_writesac(sacfile, header, this%seismo_dat(:, 2, irec_local), ier)
            sacfile = trim(OUTPUT_FILES)//'/syn.'//trim(this%od%netwk(irec))//&
                      '.'//trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
                      'R.sac.'//trim(this%band_name)
            call sacio_writesac(sacfile, header, this%seismo_syn(:, 2, irec_local), ier)
          end block
        endif

        call measure_adj_cross_conv(this%seismo_dat(:, 2, irec_local), this%seismo_dat(:, 1, irec_local), &
                                    this%seismo_syn(:, 2, irec_local), this%seismo_syn(:, 1, irec_local), &
                                    dble(fpar%sim%TIME_WIN(1)), dble(fpar%sim%TIME_WIN(2)), dble(this%ttp(irec)), &
                                    window_chi_local, adj_r_tw, adj_z_tw)

        adj_r_tw = adj_r_tw * fpar%acqui%src_weight(this%ievt)
        adj_z_tw = adj_z_tw * fpar%acqui%src_weight(this%ievt)
        call this%write_adj(adj_z_tw(1:NSTEP), trim(this%comp_name(1)), irec)
        adj_src = zeros_dp(NSTEP, 2)
        call rotate_R_to_NE_dp(adj_r_tw(1:NSTEP), adj_src(:, 2), adj_src(:, 1), this%baz)
        call this%write_adj(adj_src(:, 1), trim(this%comp_name(2)), irec)
        call this%write_adj(adj_src(:, 2), trim(this%comp_name(3)), irec)

        this%window_chi(irec_local, :, 1) = window_chi_local
        this%tr_chi(irec_local, 1) = fpar%acqui%src_weight(this%ievt)*window_chi_local(15)
        this%am_chi(irec_local, 1) = fpar%acqui%src_weight(this%ievt)*window_chi_local(15)
        this%T_pmax_dat(irec_local, 1) = 0.
        this%T_pmax_syn(irec_local, 1) = 0.
        this%sta(irec_local) = trim(this%od%stnm(irec))
        this%net(irec_local) = trim(this%od%netwk(irec))
        this%tstart(irec_local) = this%ttp(irec) + fpar%sim%time_win(1)
        this%tend(irec_local) = this%ttp(irec) + fpar%sim%time_win(2)
      enddo
    endif
    call synchronize_all()
  end subroutine measure_adj

  subroutine calc_times(this)
    class(TeleCCData), intent(inout) :: this
    real(kind=cr), dimension(:), allocatable :: ttp_local
    integer :: irec_local, irec

    call read_fk_model(fpar%acqui%evtid_names(this%ievt))
    
    this%baz = -phi_FK - 90.d0
    this%az = 90.d0 - phi_FK

    call free_fk_arrays()

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
    class(TeleCCData), intent(inout) :: this

    call this%od%finalize()
    call this%wchi(1)%finalize()
    call free_shm_array(this%ttp_win)
    call free_shm_array(this%dat_win)
    deallocate(this%seismo_dat, stat=ier)
    deallocate(this%seismo_syn, stat=ier)
    deallocate(this%window_chi, stat=ier)
    deallocate(this%tr_chi, stat=ier)
    deallocate(this%am_chi, stat=ier)
    deallocate(this%T_pmax_dat, stat=ier)
    deallocate(this%T_pmax_syn, stat=ier)
    deallocate(this%sta, stat=ier)
    deallocate(this%net, stat=ier)
    deallocate(this%tstart, stat=ier)
    deallocate(this%tend, stat=ier)
  end subroutine finalize

end module telecc_data