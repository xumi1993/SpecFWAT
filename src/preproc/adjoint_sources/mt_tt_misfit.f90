module mt_tt_misfit
  use config
  use signal
  use adj_config, cfg => adj_config_global
  use fftpack
  use cc_tt_misfit, only: CCTTMisfit
  use cross_correlate
  use dpss
  use, intrinsic :: ieee_arithmetic

  implicit none

  logical, save :: is_verbose = .false.
  type, extends(CCTTMisfit) :: MTTTMisfit
    real(kind=dp) :: dt, tlen
    integer :: nlen_f, nlen
    logical :: is_mtm
  contains
    procedure :: calc_adjoint_source
    procedure, private :: initialize, check_time_series_acceptability, &
                          prepare_data_for_mtm, calculate_freq_limits, &
                          calculate_multitaper, calculate_mt_adjsrc, &
                          calculate_mt_error, check_mtm_time_shift_acceptability, &
                          calculate_freq_domain_taper
  end type MTTTMisfit

contains

  subroutine calc_adjoint_source(this, dat, syn, dt, windows)
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(in) :: windows
    real(kind=dp), dimension(:,:), allocatable :: tapers
    real(kind=dp), dimension(:), allocatable :: s, d, adj_tw_p, adj_tw_q, freq, wvec, &
                                                eig, phi_mtm, abs_mtm, dtau_mtm, dlna_mtm,&
                                                err_phi, err_abs, err_dtau, err_dlna, &
                                                wp_w, wq_w
    real(kind=dp) :: df
    integer :: iwin, nb, ne, nlen_win, nfreq_min, nfreq_max, i
    logical :: is_mtm
    type(fft_cls) :: fftins

    ! num of measurements
    if (size(windows, 2) /= 2) then
      write(*,*) 'Error: windows must have two columns (start and end times)'
      error stop
    end if

    this%nwin = size(windows, 1)
    
    ! allocate measurement arrays
    call this%initialize()

    this%nlen = size(dat)
    this%nlen_f = 2 ** exponent(real(this%nlen))
    this%dt = dt
    this%tlen = this%nlen * dt
    if (allocated(this%adj_src)) deallocate(this%adj_src)
    allocate(this%adj_src(this%nlen))
    this%adj_src = 0.0_dp

    !loop over windows
    do iwin = 1, this%nwin
      is_mtm = .true.

      ! trim the windows to the length of the data
      call get_window_info(windows(iwin,:), dt, nb, ne, nlen_win)
      s = syn(nb:ne)
      d = dat(nb:ne)

      ! taper the windows
      call window_taper(s, cfg%taper_percentage, cfg%itaper_type)
      call window_taper(d, cfg%taper_percentage, cfg%itaper_type)

      ! calculate cross-correlation shift
      call calc_cc_shift(d, s, dt, cfg%dt_sigma_min, cfg%dlna_sigma_min, &
                         this%tshift(iwin), this%dlna(iwin), &
                         this%sigma_dt(iwin), this%sigma_dlna(iwin))
      
      this%cc_max(iwin) = cc_max_coef(d, s)
      call this%cc_measure_select(iwin)
      if (.not. this%select_meas(iwin)) cycle
      do while (is_mtm)
        ! prepare data for MTM
        is_mtm = this%check_time_series_acceptability(iwin, nlen_win)
        if (.not. is_mtm) then
          if (is_verbose) write(*,*) 'Warning: Window ', iwin, ' does not meet time series criteria for MTM.'
          exit
        end if

        call this%prepare_data_for_mtm(dat, iwin, windows(iwin,:), d, is_mtm)
        if (.not. is_mtm) then
          if (is_verbose) write(*,*) 'Warning: Window ', iwin, ' could not be prepared for MTM.'
          exit
        end if

        freq = fftins%fftfreq(this%nlen_f, dt)
        df = freq(2) - freq(1)
        wvec = 2.0_dp * PI * freq

        ! calculate frequency limits
        call this%calculate_freq_limits(syn, df, nfreq_min, nfreq_max, is_mtm)
        if (.not. is_mtm) then
          if (is_verbose) write(*,*) 'Warning: Window ', iwin, ' does not meet frequency criteria for MTM.'
          exit
        end if

        ! calculate multitaper transfer function and measurements
        call dpss_windows(nlen_win, cfg%mt_nw, cfg%num_taper, tapers, eig)
        do i = 1, cfg%num_taper
          if (eig(i) > HUGEVAL) then
            if (is_verbose) write(*,*) 'Warning: DPSS taper ', i, ' has infinite eigenvalue:', eig(i),'. Skipping MTM for window ', iwin
            is_mtm = .false.
            exit
          end if
        end do
        if (.not. is_mtm) exit

        tapers = transpose(tapers) * dsqrt(dble(nlen_win))
        call this%calculate_multitaper(d, s, tapers, wvec, nfreq_min, nfreq_max, &
                                      this%tshift(iwin), this%dlna(iwin), &
                                      phi_mtm, abs_mtm, dtau_mtm, dlna_mtm)
        
        if (cfg%use_mt_error) then 
          call this%calculate_mt_error(d, s, tapers, wvec, nfreq_min, nfreq_max, &
                                      this%tshift(iwin), this%dlna(iwin), phi_mtm, abs_mtm, dtau_mtm, dlna_mtm, &
                                      err_phi, err_abs, err_dtau, err_dlna)
        else
          allocate(err_dtau(this%nlen_f), err_dlna(this%nlen_f))
          err_dtau = 0.0_dp
          err_dlna = 0.0_dp
        end if
        ! check MTM measurement quality before calculating adjoint source
        is_mtm = this%check_mtm_time_shift_acceptability(nfreq_min, nfreq_max, df, &
                                                         this%tshift(iwin), dtau_mtm, err_dtau)
        if (.not. is_mtm) exit

        call this%calculate_freq_domain_taper(nfreq_min, nfreq_max, df, &
                                              dtau_mtm, dlna_mtm, &
                                              this%sigma_dt(iwin), this%sigma_dlna(iwin), &
                                              err_dtau, err_dlna, wp_w, wq_w)
        ! Store the misfits for this window
        this%misfit_p(iwin) = simpson(dtau_mtm ** 2 * wp_w, df)
        this%misfit_q(iwin) = simpson(dlna_mtm ** 2 * wq_w, df)

        ! calculate multitaper phase shift misfit and adjoint source
        call this%calculate_mt_adjsrc(s, tapers, nfreq_min, nfreq_max, &
                                      dtau_mtm, dlna_mtm, wp_w, wq_w, adj_tw_p, adj_tw_q)
        this%imeas(iwin) = cfg%imeasure_type  
        exit

      end do ! end of is_mtm loop

      if (.not. is_mtm) then
        ! calculate adjoint source
        call calc_cc_adjsrc(s, this%tshift(iwin), this%dlna(iwin), dt, &
                            this%sigma_dt(iwin), this%sigma_dlna(iwin), &
                            this%misfit_p(iwin), this%misfit_q(iwin), &
                            adj_tw_p, adj_tw_q)
        select case (cfg%imeasure_type)
          case(IMEAS_CC_TT_MT) ! MT-CC-TT (phase shift)
            this%imeas(iwin) = IMEAS_CC_TT
          case(IMEAS_CC_DLNA_MT) ! MT-CC-DLNA (amplitude)
            this%imeas(iwin) = IMEAS_CC_DLNA
        end select
      endif

      call window_taper(adj_tw_p, cfg%taper_percentage, cfg%itaper_type)
      call window_taper(adj_tw_q, cfg%taper_percentage, cfg%itaper_type)

      ! select measurements based on measurement type and add to total adjoint source
      select case (cfg%imeasure_type)
        case(IMEAS_CC_TT_MT) ! MT-CC-TT (phase shift)
          this%total_misfit = this%total_misfit + this%misfit_p(iwin)
          this%misfits(iwin) = this%misfit_p(iwin)
          this%residuals(iwin) = this%tshift(iwin)
          this%errors(iwin) = this%sigma_dt(iwin)
          this%adj_src(nb:ne) = this%adj_src(nb:ne) + adj_tw_p(:)
        case(IMEAS_CC_DLNA_MT) ! MT-CC-DLNA (amplitude)
          this%total_misfit = this%total_misfit + this%misfit_q(iwin)
          this%misfits(iwin) = this%misfit_q(iwin)
          this%residuals(iwin) = this%dlna(iwin)
          this%errors(iwin) = this%sigma_dlna(iwin)
          this%adj_src(nb:ne) = this%adj_src(nb:ne) + adj_tw_q(:)
      end select
      
      ! deallocate tmp variables
      if (allocated(freq)) deallocate(freq)
      if (allocated(s)) deallocate(s)
      if (allocated(d)) deallocate(d)
      if (allocated(adj_tw_p)) deallocate(adj_tw_p)
      if (allocated(adj_tw_q)) deallocate(adj_tw_q)
      if (allocated(wvec)) deallocate(wvec)
      if (allocated(eig)) deallocate(eig)
      if (allocated(phi_mtm)) deallocate(phi_mtm)
      if (allocated(abs_mtm)) deallocate(abs_mtm)
      if (allocated(dtau_mtm)) deallocate(dtau_mtm)
      if (allocated(dlna_mtm)) deallocate(dlna_mtm)
      if (allocated(err_phi)) deallocate(err_phi)
      if (allocated(err_abs)) deallocate(err_abs)
      if (allocated(err_dtau)) deallocate(err_dtau)
      if (allocated(err_dlna)) deallocate(err_dlna)
      if (allocated(wp_w)) deallocate(wp_w)
      if (allocated(wq_w)) deallocate(wq_w)
      if (allocated(tapers)) deallocate(tapers)
    end do
  end subroutine calc_adjoint_source

  subroutine calculate_multitaper(this, d, s, tapers, wvec, nfreq_min, nfreq_max, cc_tshift, cc_dlna,&
                                  phi_w, abs_w, dtau_w, dlna_w)
    ! Calculate the multitaper transfer function and measurements
    ! Arguments:
    ! d, s: time series of data and synthetics
    ! tapers: array of tapers (nlen_t, ntaper)
    ! wvec: angular frequency vector (nlen_f)
    ! nfreq_min, nfreq_max: frequency indices for analysis
    ! cc_tshift, cc_dlna: cross-correlation time shift and log amplitude ratio
    ! Outputs:
    ! phi_w: phase of transfer function
    ! abs_w: amplitude of transfer function
    ! dtau_w: derivative of phase with respect to frequency
    ! dlna_w: derivative of log amplitude with respect to frequency
    
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: d, s, wvec
    real(kind=dp), dimension(:,:), intent(in) :: tapers
    integer, intent(in) :: nfreq_min, nfreq_max
    real(kind=dp), intent(in) :: cc_tshift, cc_dlna
    real(kind=dp), dimension(:), allocatable, intent(out) :: phi_w, abs_w, dtau_w, dlna_w
    complex(kind=dp), dimension(:), allocatable :: top_tf, bot_tf, d_tw, s_tw, trans_func
    real(kind=dp), dimension(:), allocatable :: taper, d_t, s_t
    real(kind=dp) :: wlevel
    integer :: nlen_t, ntaper, fnum, itaper, i
    type(fft_cls) :: fftins
    
    nlen_t = size(d)
    ntaper = size(tapers, 2)
    fnum = int(this%nlen_f / 2) + 1

    allocate(top_tf(this%nlen_f), bot_tf(this%nlen_f), trans_func(this%nlen_f))
    allocate(d_t(this%nlen_f), s_t(this%nlen_f))
    top_tf = (0.0_dp, 0.0_dp)
    bot_tf = (0.0_dp, 0.0_dp)
    trans_func = (0.0_dp, 0.0_dp)  ! Initialize trans_func

    do itaper = 1, ntaper
      taper = tapers(1:nlen_t, itaper)

      ! Apply taper to data and synthetics
      d_t = 0.0_dp
      s_t = 0.0_dp
      d_t(1:nlen_t) = d * taper
      s_t(1:nlen_t) = s * taper

      ! Prepare zero-padded arrays for FFT

      d_tw = fftins%fft_dp(d_t, this%nlen_f) * this%dt
      s_tw = fftins%fft_dp(s_t, this%nlen_f) * this%dt

      ! Accumulate the top and bottom terms for multitaper
      top_tf = top_tf + d_tw * conjg(s_tw)
      bot_tf = bot_tf + s_tw * conjg(s_tw)
    enddo

    ! Calculate the transfer function with water level stabilization
    wlevel = maxval(abs(bot_tf(1:fnum))) * cfg%transfunc_waterlevel**2

    ! Calculate the transfer function
    do i = nfreq_min, nfreq_max
      if (abs(bot_tf(i)) >= wlevel) then
        trans_func(i) = top_tf(i) / bot_tf(i)
      else
        trans_func(i) = top_tf(i) / (bot_tf(i) + wlevel)
      end if
    enddo

    allocate(phi_w(this%nlen_f), abs_w(this%nlen_f), dtau_w(this%nlen_f), dlna_w(this%nlen_f))
    phi_w = 0.0_dp
    abs_w = 0.0_dp
    dtau_w = 0.0_dp
    dlna_w = 0.0_dp

    ! Calculate phase, amplitude, and their derivatives
    phi_w(nfreq_min:nfreq_max) = atan2(aimag(trans_func(nfreq_min:nfreq_max)), real(trans_func(nfreq_min:nfreq_max)))
    call process_cycle_skipping(phi_w, nfreq_min, nfreq_max, wvec(1:this%nlen_f), this%nlen_f)

    ! Calculate amplitude
    abs_w(nfreq_min:nfreq_max) = abs(trans_func(nfreq_min:nfreq_max))

    ! Add CC measurement to transfer function
    dtau_w(1) = cc_tshift
    dtau_w(max(2, nfreq_min):nfreq_max) = &
      - 1.0_dp / wvec(max(2, nfreq_min):nfreq_max) * &
      phi_w(max(2, nfreq_min):nfreq_max) + cc_tshift

    dlna_w(nfreq_min:nfreq_max) = log(abs_w(nfreq_min:nfreq_max)) + cc_dlna

  end subroutine calculate_multitaper

  subroutine calculate_mt_error(this, d, s, tapers, wvec, nfreq_min, nfreq_max, &
                               cc_tshift, cc_dlna, phi_mtm, abs_mtm, dtau_mtm, dlna_mtm, &
                               err_phi, err_abs, err_dtau, err_dlna)
    ! Calculate multitaper error with Jackknife MT estimates
    ! Uses jackknife method by systematically leaving out each observation
    
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: d, s, wvec
    real(kind=dp), dimension(:,:), intent(in) :: tapers
    integer, intent(in) :: nfreq_min, nfreq_max
    real(kind=dp), intent(in) :: cc_tshift, cc_dlna
    real(kind=dp), dimension(:), intent(in) :: phi_mtm, abs_mtm, dtau_mtm, dlna_mtm
    real(kind=dp), dimension(:), allocatable, intent(out) :: err_phi, err_abs, err_dtau, err_dlna
    
    real(kind=dp), dimension(:,:), allocatable :: phi_mul, abs_mul, dtau_mul, dlna_mul, tapers_om
    real(kind=dp), dimension(:), allocatable :: ephi_ave, eabs_ave, edtau_ave, edlna_ave
    real(kind=dp), dimension(:), allocatable :: phi_om, abs_om, dtau_om, dlna_om
    integer :: nlen_t, ntaper, itaper, j, icol
    
    nlen_t = size(d)
    ntaper = size(tapers, 2)
    
    ! Initialize arrays
    allocate(phi_mul(this%nlen_f, ntaper), abs_mul(this%nlen_f, ntaper))
    allocate(dtau_mul(this%nlen_f, ntaper), dlna_mul(this%nlen_f, ntaper))
    allocate(ephi_ave(this%nlen_f), eabs_ave(this%nlen_f))
    allocate(edtau_ave(this%nlen_f), edlna_ave(this%nlen_f))
    allocate(err_phi(this%nlen_f), err_abs(this%nlen_f))
    allocate(err_dtau(this%nlen_f), err_dlna(this%nlen_f))
    allocate(tapers_om(nlen_t, ntaper-1))
    
    phi_mul = 0.0_dp
    abs_mul = 0.0_dp
    dtau_mul = 0.0_dp
    dlna_mul = 0.0_dp
    ephi_ave = 0.0_dp
    eabs_ave = 0.0_dp
    edtau_ave = 0.0_dp
    edlna_ave = 0.0_dp
    err_phi = 0.0_dp
    err_abs = 0.0_dp
    err_dtau = 0.0_dp
    err_dlna = 0.0_dp
    
    ! Jackknife loop: leave out one taper at a time
    do itaper = 1, ntaper
      ! Create taper array with one taper removed
      tapers_om = 0.0_dp
      icol = 1
      do j = 1, ntaper
        if (j /= itaper) then
          tapers_om(1:this%nlen_f, icol) = tapers(1:this%nlen_f, j)
          icol = icol + 1
        end if
      end do
      
      ! Recalculate MT measurements with reduced taper set
      call this%calculate_multitaper(d, s, tapers_om, wvec, nfreq_min, nfreq_max, &
                                    cc_tshift, cc_dlna, phi_om, abs_om, dtau_om, dlna_om)
      ! Store jackknife estimates
      phi_mul(:, itaper) = phi_om(:)
      abs_mul(:, itaper) = abs_om(:)
      dtau_mul(:, itaper) = dtau_om(:)
      dlna_mul(:, itaper) = dlna_om(:)
      
      ! Error estimation using jackknife formula
      ephi_ave(nfreq_min:nfreq_max) = ephi_ave(nfreq_min:nfreq_max) + &
        ntaper * phi_mtm(nfreq_min:nfreq_max) - &
        (ntaper - 1) * phi_mul(nfreq_min:nfreq_max, itaper)
        
      eabs_ave(nfreq_min:nfreq_max) = eabs_ave(nfreq_min:nfreq_max) + &
        ntaper * abs_mtm(nfreq_min:nfreq_max) - &
        (ntaper - 1) * abs_mul(nfreq_min:nfreq_max, itaper)
        
      edtau_ave(nfreq_min:nfreq_max) = edtau_ave(nfreq_min:nfreq_max) + &
        ntaper * dtau_mtm(nfreq_min:nfreq_max) - &
        (ntaper - 1) * dtau_mul(nfreq_min:nfreq_max, itaper)
        
      edlna_ave(nfreq_min:nfreq_max) = edlna_ave(nfreq_min:nfreq_max) + &
        ntaper * dlna_mtm(nfreq_min:nfreq_max) - &
        (ntaper - 1) * dlna_mul(nfreq_min:nfreq_max, itaper)
    end do
    
    ! Take average over each taper band
    ephi_ave = ephi_ave / ntaper
    eabs_ave = eabs_ave / ntaper
    edtau_ave = edtau_ave / ntaper
    edlna_ave = edlna_ave / ntaper
    
    ! Calculate deviation
    do itaper = 1, ntaper
      err_phi(nfreq_min:nfreq_max) = err_phi(nfreq_min:nfreq_max) + &
        (phi_mul(nfreq_min:nfreq_max, itaper) - ephi_ave(nfreq_min:nfreq_max))**2
      err_abs(nfreq_min:nfreq_max) = err_abs(nfreq_min:nfreq_max) + &
        (abs_mul(nfreq_min:nfreq_max, itaper) - eabs_ave(nfreq_min:nfreq_max))**2
      err_dtau(nfreq_min:nfreq_max) = err_dtau(nfreq_min:nfreq_max) + &
        (dtau_mul(nfreq_min:nfreq_max, itaper) - edtau_ave(nfreq_min:nfreq_max))**2
      err_dlna(nfreq_min:nfreq_max) = err_dlna(nfreq_min:nfreq_max) + &
        (dlna_mul(nfreq_min:nfreq_max, itaper) - edlna_ave(nfreq_min:nfreq_max))**2
    end do
    
    ! Calculate standard deviation
    err_phi(nfreq_min:nfreq_max) = sqrt(err_phi(nfreq_min:nfreq_max) / (ntaper * (ntaper - 1)))
    err_abs(nfreq_min:nfreq_max) = sqrt(err_abs(nfreq_min:nfreq_max) / (ntaper * (ntaper - 1)))
    err_dtau(nfreq_min:nfreq_max) = sqrt(err_dtau(nfreq_min:nfreq_max) / (ntaper * (ntaper - 1)))
    err_dlna(nfreq_min:nfreq_max) = sqrt(err_dlna(nfreq_min:nfreq_max) / (ntaper * (ntaper - 1)))
    
  end subroutine calculate_mt_error

  subroutine calculate_mt_adjsrc(this, s, tapers, nfreq_min, nfreq_max, &
                                dtau_mtm, dlna_mtm, wp_w, wq_w, adj_p, adj_q)
    ! Calculate the adjoint source for a multitaper measurement, which
    ! tapers synthetics in various windowed frequency-dependent tapers and
    ! scales them by phase dependent travel time measurements (which
    ! incorporate the observed data).
    
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: s, dtau_mtm, dlna_mtm, wp_w, wq_w
    real(kind=dp), dimension(:,:), intent(in) :: tapers
    integer, intent(in) :: nfreq_min, nfreq_max
    real(kind=dp), dimension(:), allocatable, intent(out) :: adj_p, adj_q
    
    real(kind=dp), dimension(:), allocatable :: s_t, s_tv, fp_t, fq_t
    real(kind=dp), dimension(:), allocatable :: p_wt, q_wt, taper_full
    complex(kind=dp), dimension(:), allocatable :: bottom_p, bottom_q
    complex(kind=dp), dimension(:,:), allocatable :: s_tw, s_tvw
    complex(kind=dp), dimension(:), allocatable :: p_w, q_w
    complex, dimension(:), allocatable :: p_wt_c, q_wt_c
    integer :: nlen_t, ntaper, itaper, i
    type(fft_cls) :: fftins
    
    nlen_t = size(s)
    ntaper = size(tapers, 2)
    
    ! Allocate arrays
    allocate(fp_t(nlen_t), fq_t(nlen_t))
    allocate(bottom_p(this%nlen_f), bottom_q(this%nlen_f))
    allocate(s_tw(this%nlen_f, ntaper), s_tvw(this%nlen_f, ntaper))
    allocate(p_w(this%nlen_f), q_w(this%nlen_f))
    allocate(p_wt_c(this%nlen_f), q_wt_c(this%nlen_f))
    allocate(p_wt(this%nlen_f), q_wt(this%nlen_f))
    allocate(taper_full(this%nlen_f))
    allocate(adj_p(nlen_t), adj_q(nlen_t))
    allocate(s_t(this%nlen_f), s_tv(this%nlen_f))
    
    ! Start piecing together transfer functions that will be applied to synthetics
    bottom_p = (0.0_dp, 0.0_dp)
    bottom_q = (0.0_dp, 0.0_dp)
    s_tw = (0.0_dp, 0.0_dp)
    s_tvw = (0.0_dp, 0.0_dp)
    
    ! Construct the bottom term of the adjoint formula which requires
    ! summed contributions from each of the taper bands
    do itaper = 1, ntaper
      taper_full = 0.0_dp
      taper_full(1:nlen_t) = tapers(1:nlen_t, itaper)
      
      ! Taper synthetics (s_t) and take the derivative (s_tv)
      s_t = 0.0_dp
      s_tv = 0.0_dp
      s_t(1:nlen_t) = s * taper_full(1:nlen_t)
      s_tv(1:nlen_t) = gradient(s_t(1:nlen_t), this%dt)
      
      ! Apply FFT to tapered measurements to get to freq. domain.
      s_tw(:, itaper) = fftins%fft_dp(s_t, this%nlen_f) * this%dt
      s_tvw(:, itaper) = fftins%fft_dp(s_tv, this%nlen_f) * this%dt
      
      ! Calculate bottom term of the adjoint equation
      bottom_p = bottom_p + s_tvw(:, itaper) * conjg(s_tvw(:, itaper))
      bottom_q = bottom_q + s_tw(:, itaper) * conjg(s_tw(:, itaper))
    end do
    
    ! Now we generate the adjoint sources using each of the tapers
    fp_t = 0.0_dp
    fq_t = 0.0_dp
    
    do itaper = 1, ntaper
      taper_full = 0.0_dp
      taper_full(1:nlen_t) = tapers(1:nlen_t, itaper)
      
      ! Calculate the full adjoint terms p_w, q_w
      p_w = (0.0_dp, 0.0_dp)
      q_w = (0.0_dp, 0.0_dp)
      
      do i = nfreq_min, nfreq_max
        if (abs(bottom_p(i)) > 0.0_dp) then
          p_w(i) = s_tvw(i, itaper) / bottom_p(i)
        end if
        if (abs(bottom_q(i)) > 0.0_dp) then
          q_w(i) = -s_tw(i, itaper) / bottom_q(i)
        end if
      end do
      
      ! Weight the adjoint terms by the phase + amplitude measurements
      do i = 1, this%nlen_f
        p_w(i) = p_w(i) * dtau_mtm(i) * wp_w(i)  ! phase
        q_w(i) = q_w(i) * dlna_mtm(i) * wq_w(i)  ! amplitude
      end do
      
      ! Inverse FFT of weighted adjoint to get back to the time domain
      p_wt_c = cmplx(p_w)
      q_wt_c = cmplx(q_w)
      
      p_wt = real(fftins%ifft(p_wt_c, this%nlen_f)) * 2.0_dp / this%dt
      q_wt = real(fftins%ifft(q_wt_c, this%nlen_f)) * 2.0_dp / this%dt

      ! Taper adjoint term before adding it back to full adj source
      fp_t = fp_t + p_wt(1:nlen_t) * taper_full(1:nlen_t)
      fq_t = fq_t + q_wt(1:nlen_t) * taper_full(1:nlen_t)
    end do

    
    ! Return the adjoint sources
    adj_p = fp_t
    adj_q = fq_t
    
  end subroutine calculate_mt_adjsrc

  function check_mtm_time_shift_acceptability(this, nfreq_min, nfreq_max, df, &
                                             cc_tshift, dtau_mtm, err_dtau) result(is_mtm)
    ! Check MTM time shift measurements to see if they are within allowable
    ! bounds set by the config. If any of the phases used in MTM do not
    ! meet criteria, we will fall back to CC measurement.
    
    class(MTTTMisfit), intent(inout) :: this
    integer, intent(in) :: nfreq_min, nfreq_max
    real(kind=dp), intent(in) :: df, cc_tshift
    real(kind=dp), dimension(:), intent(in) :: dtau_mtm, err_dtau
    logical :: is_mtm
    
    real(kind=dp) :: wave_period, max_dt_allowed, max_err_allowed
    integer :: j
    
    ! True unless set False
    is_mtm = .true.
    
    ! If any MTM measurements is out of the reasonable range, switch to CC
    do j = nfreq_min, nfreq_max
      if (j * df > 0.0_dp) then
        wave_period = 1.0_dp / (j * df)
        
        ! dt larger than 1/dt_fac of the wave period
        max_dt_allowed = wave_period / cfg%dt_fac
        if (abs(dtau_mtm(j)) > max_dt_allowed) then
          if (is_verbose) write(*,*) 'Warning: reject MTM: dt measurement is too large at frequency', j
          is_mtm = .false.
          exit
        end if
        
        ! Error larger than 1/err_fac of wave period
        max_err_allowed = wave_period / cfg%err_fac
        if (err_dtau(j) > max_err_allowed) then
          if (is_verbose) write(*,*) 'Warning: reject MTM: dt error is too large at frequency', j
          is_mtm = .false.
          exit
        end if
        
        ! dt larger than the maximum allowable time shift
        if (abs(dtau_mtm(j)) > cfg%dt_max_scale * abs(cc_tshift)) then
          if (is_verbose) write(*,*) 'Warning: reject MTM: dt is larger than maximum '//&
                                      'allowable time shift at frequency', j
          is_mtm = .false.
          exit
        end if
      end if
    end do
    
  end function check_mtm_time_shift_acceptability
    
  function check_time_series_acceptability(this, iwin, nlen_w) result(is_acceptable)
    class(MTTTMisfit), intent(in) :: this
    integer, intent(in) :: nlen_w, iwin
    logical :: is_acceptable

    ! Check if the time shift is within acceptable limits
    if (abs(this%tshift(iwin)) <= this%dt) then
      is_acceptable = .false.
    elseif (cfg%min_cycle_in_window * cfg%min_period > nlen_w) then
      is_acceptable = .false.
    else
      is_acceptable = .true.
    end if
  end function check_time_series_acceptability

  subroutine prepare_data_for_mtm(this, dat, iwin, window, d, is_acceptable)
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(out) :: d
    real(kind=dp), dimension(:), intent(in) :: dat
    real(kind=dp), dimension(:), intent(in) :: window
    integer, intent(in) :: iwin
    integer :: nb, ne, ishift, nb_d, ne_d, nlen_d, nlen_w
    logical, intent(out) :: is_acceptable

    ! Prepare data for MTM analysis
    call get_window_info(window, this%dt, nb, ne, nlen_w)
    ishift = int(this%tshift(iwin) / this%dt)

    nb_d = max(1, nb + ishift)
    ne_d = min(this%nlen, ne + ishift)
    nlen_d = ne_d - nb_d + 1

    if (nlen_d == nlen_w) then
      ! If the data length matches the window length, use it directly
      d(1:nlen_w) = dat(nb_d:ne_d)
      d = d * exp(-this%dlna(iwin))
      call window_taper(d, cfg%taper_percentage, cfg%itaper_type)
      is_acceptable = .true.
    else
      d = dat(nb:ne)
      is_acceptable = .false.
    end if
  end subroutine prepare_data_for_mtm

  subroutine calculate_freq_limits(this, syn, df, nfreq_min, nfreq_max, is_acceptable)
    class(MTTTMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: syn
    real(kind=dp), intent(in) :: df
    integer, intent(out) :: nfreq_min, nfreq_max
    logical, intent(out) :: is_acceptable
    complex(kind=dp), dimension(:), allocatable :: s_spec
    real(kind=dp) :: ampmax, scaled_wl, half_taper_bandwidth, &
                     chosen_bandwidth
    real(kind=dp), dimension(:), allocatable :: syn_expand
    integer :: fnum, i_ampmax, ifreq_min, ifreq_max, iw
    logical :: is_search
    type(fft_cls) :: fftins

    ! calculate frequency limits for MTM analysis using synthetic data
    fnum = int(this%nlen_f / 2) + 1
    
    ! Use synthetic data (not adj_src) to match Python implementation
    allocate(syn_expand(this%nlen_f))
    syn_expand = 0.0_dp
    syn_expand(1:size(syn)) = syn(:)
    s_spec = fftins%fft_dp(syn_expand, this%nlen_f) * this%dt

    ! Calculate the frequency limits based on the sampling rate and window length
    ! Only consider the positive frequencies (1:fnum)
    ampmax = maxval(abs(s_spec(1:fnum)))
    i_ampmax = maxloc(abs(s_spec(1:fnum)), dim=1)

    ! Scale the maximum amplitude to the water threshold
    scaled_wl = cfg%water_threshold * ampmax

    ! Find the frequency index corresponding to the maximum amplitude
    ifreq_min = max(1, int(1.0_dp / (cfg%max_period * df)))  ! Ensure minimum is 1
    ifreq_max = min(fnum, int(1.0_dp / (cfg%min_period * df)))  ! Ensure maximum doesn't exceed fnum

    ! write(*,*) 'DEBUG: fnum=', fnum, ', i_ampmax=', i_ampmax
    ! write(*,*) 'DEBUG: ifreq_min=', ifreq_min, ', ifreq_max=', ifreq_max
    ! write(*,*) 'DEBUG: ampmax=', ampmax, ', scaled_wl=', scaled_wl

    ! get the maximum frequency limits
    nfreq_max = fnum  ! Start from fnum (will be reduced by search)
    is_search = .true.
    if (i_ampmax < fnum) then  ! Only search if there's space to search
      do iw = i_ampmax + 1, fnum
        call search_frequency_limit(is_search, iw, nfreq_max, s_spec, scaled_wl, 10)
      end do
    end if
    ! Make sure `nfreq_max` does not go beyond reasonable limits
    nfreq_max = min(nfreq_max, ifreq_max, int(1.0_dp / (2.0_dp * this%dt) / df))

    ! get the minimum frequency limits  
    nfreq_min = 1  ! Start from 1 (will be increased by search)
    is_search = .true.
    if (i_ampmax > 1) then  ! Only search if there's space to search
      do iw = i_ampmax - 1, 1, -1
        call search_frequency_limit(is_search, iw, nfreq_min, s_spec, scaled_wl, 10)
      end do
    end if
    ! Make sure `nfreq_min` does not go below reasonable limits
    nfreq_min = max(nfreq_min, ifreq_min, max(1, int(cfg%min_cycle_in_window / this%tlen / df)))

    half_taper_bandwidth = cfg%mt_nw / (4.0_dp * this%tlen)
    chosen_bandwidth = (nfreq_max - nfreq_min) * df
    if (chosen_bandwidth < half_taper_bandwidth) then
      is_acceptable = .false.
    else
      is_acceptable = .true.
    end if

  end subroutine calculate_freq_limits

  subroutine search_frequency_limit(is_search, index, nfreq_limit, spectra, water_threshold, c)
    integer, intent(in) :: index          ! index of the frequency to check
    integer, intent(in) :: c              ! scaling factor for water threshold
    integer, intent(inout) :: nfreq_limit ! frequency limit index
    real(kind=dp), intent(in) :: water_threshold ! water threshold value
    complex(kind=dp), intent(in) :: spectra(:)     ! spectrum data
    logical, intent(inout) :: is_search   ! search flag

    ! First condition: stop search if amplitude is below water threshold
    if (abs(spectra(index)) < water_threshold .and. is_search) then
      is_search   = .false.
      nfreq_limit = index
    end if

    ! Second condition: restart search if amplitude is above c * water_threshold  
    ! This works like a "heating thermostat" as described in Python docstring
    if (abs(spectra(index)) > c * water_threshold .and. .not. is_search) then
      is_search   = .true.
      nfreq_limit = index
    end if

  end subroutine search_frequency_limit

  subroutine calculate_freq_domain_taper(this, nfreq_min, nfreq_max, df, &
                                        dtau_mtm, dlna_mtm, err_dt_cc, err_dlna_cc, &
                                        err_dtau_mt, err_dlna_mt, wp_w, wq_w)
    ! Calculate frequency domain taper weighted by misfit (either CC or MTM)
    ! Frequency-domain tapers are based on adjusted frequency band and error estimation.
    ! They are not one of the filtering processes that needs to be applied to the 
    ! adjoint source but rather a frequency domain weighting function for adjoint 
    ! source and misfit function.
    
    class(MTTTMisfit), intent(inout) :: this
    integer, intent(in) :: nfreq_min, nfreq_max
    real(kind=dp), intent(in) :: df, err_dt_cc, err_dlna_cc
    real(kind=dp), dimension(:), intent(in) :: dtau_mtm, dlna_mtm, err_dtau_mt, err_dlna_mt
    real(kind=dp), dimension(:), allocatable, intent(out) :: wp_w, wq_w
    
    real(kind=dp), dimension(:), allocatable :: w_taper, win_taper
    real(kind=dp), dimension(:), allocatable :: err_dtau_mt_work, err_dlna_mt_work
    real(kind=dp) :: ffac, dtau_wtr, dlna_wtr
    integer :: win_taper_len, i
    
    ! Initialize arrays
    allocate(w_taper(this%nlen_f))
    allocate(err_dtau_mt_work(this%nlen_f))
    allocate(err_dlna_mt_work(this%nlen_f))
    w_taper = 0.0_dp
    err_dtau_mt_work = err_dtau_mt
    err_dlna_mt_work = err_dlna_mt
    
    ! Create frequency domain taper
    win_taper_len = nfreq_max - nfreq_min
    allocate(win_taper(win_taper_len))
    win_taper = 1.0_dp
    
    ! Create cosine taper over frequency range
    call window_taper(win_taper, 1.0_dp, 3) ! cos taper with 100% tapering
    w_taper(nfreq_min:nfreq_max-1) = win_taper(1:win_taper_len)
    
    ! Normalization factor, factor 2 is needed for integration -inf to inf
    ffac = 2.0_dp * df * sum(w_taper(nfreq_min:nfreq_max-1))
    
    ! write(*,*) 'DEBUG: frequency bound (idx): [', nfreq_min, ',', nfreq_max-1, ']', &
    !            ' (Hz) [', df * (nfreq_min - 1), ',', df * nfreq_max, ']'
    ! write(*,*) 'DEBUG: frequency domain taper normalization coeff:', ffac
    ! write(*,*) 'DEBUG: frequency domain sampling length df=', df
    
    ! if (ffac <= 0.0_dp) then
      ! write(*,*) 'WARNING: frequency band too narrow:'
      ! write(*,*) 'fmin=', nfreq_min, ', fmax=', nfreq_max, ', ffac=', ffac
    ! end if
    
    ! Normalized, tapered window in the frequency domain
    wp_w = w_taper / ffac
    wq_w = w_taper / ffac
    
    ! Choose whether to scale by CC error or by calculated MT errors
    if (cfg%use_cc_error) then
      wp_w = wp_w / (err_dt_cc ** 2)
      wq_w = wq_w / (err_dlna_cc ** 2)
    else if (cfg%use_mt_error) then
      ! Calculate water threshold for MT measurements
      dtau_wtr = cfg%water_threshold * sum(abs(dtau_mtm(nfreq_min:nfreq_max-1))) / &
                 (nfreq_max - nfreq_min)
      dlna_wtr = cfg%water_threshold * sum(abs(dlna_mtm(nfreq_min:nfreq_max-1))) / &
                 (nfreq_max - nfreq_min)
      
      ! Apply water threshold to error estimates where needed
      do i = nfreq_min, nfreq_max
        if (err_dtau_mt_work(i) < dtau_wtr) then
          err_dtau_mt_work(i) = err_dtau_mt_work(i) + dtau_wtr
        end if
        if (err_dlna_mt_work(i) < dlna_wtr) then
          err_dlna_mt_work(i) = err_dlna_mt_work(i) + dlna_wtr
        end if
      end do
      
      ! Scale by multitaper error estimates
      do i = nfreq_min, nfreq_max
        wp_w(i) = wp_w(i) / (err_dtau_mt_work(i) ** 2)
        wq_w(i) = wq_w(i) / (err_dlna_mt_work(i) ** 2)
      end do
    end if
    
  end subroutine calculate_freq_domain_taper

  subroutine initialize(this)
    class(MTTTMisfit), intent(inout) :: this

    ! deallocate arrays
    if (allocated(this%adj_src)) deallocate(this%adj_src)
    if (allocated(this%misfits)) deallocate(this%misfits)
    if (allocated(this%tshift)) deallocate(this%tshift)
    if (allocated(this%dlna)) deallocate(this%dlna)
    if (allocated(this%sigma_dt)) deallocate(this%sigma_dt)
    if (allocated(this%sigma_dlna)) deallocate(this%sigma_dlna)
    if (allocated(this%misfit_p)) deallocate(this%misfit_p)
    if (allocated(this%misfit_q)) deallocate(this%misfit_q)
    if (allocated(this%cc_max)) deallocate(this%cc_max)
    if (allocated(this%select_meas)) deallocate(this%select_meas)
    if (allocated(this%imeas)) deallocate(this%imeas)
    if (allocated(this%residuals)) deallocate(this%residuals)
    if (allocated(this%errors)) deallocate(this%errors)

    allocate(this%misfits(this%nwin))
    allocate(this%tshift(this%nwin))
    allocate(this%dlna(this%nwin))
    allocate(this%sigma_dt(this%nwin))
    allocate(this%sigma_dlna(this%nwin))
    allocate(this%misfit_p(this%nwin))
    allocate(this%misfit_q(this%nwin))
    allocate(this%cc_max(this%nwin))
    allocate(this%select_meas(this%nwin))
    allocate(this%imeas(this%nwin))
    allocate(this%residuals(this%nwin))
    allocate(this%errors(this%nwin))
    this%select_meas = .true.
    this%misfits = 0.0_dp
    this%total_misfit = 0.0_dp
    
    ! ! Set measurement type based on imeasure_type
    ! select case (imeasure_type)
    !   case() ! MT-CC-TT
    !     this%imeas = IMEAS_CC_TT_MT
    !   case(5) ! MT-CC-DLNA  
    !     this%imeas = IMEAS_CC_DLNA_MT
    !   case default
    !     this%imeas = IMEAS_CC_TT_MT ! default to phase shift
    ! end select

  end subroutine initialize
end module mt_tt_misfit