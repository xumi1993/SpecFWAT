module flexwin
  use config
  use signal
  use cross_correlate
  use sta_lta_mod
  use fftpack
  use adj_config
  use interval_scheduling
  implicit none

  integer, private, parameter :: max_good_win = 500000

  type :: window_candidate_type
     integer :: left_idx
     integer :: right_idx
     integer :: peak_idx
     real(kind=dp) :: cc
     real(kind=dp) :: dlnA
     real(kind=dp) :: dt_shift
     real(kind=dp) :: weight
     real(kind=dp) :: snr_amp
     real(kind=dp) :: snr_eng
     logical :: rejected
  end type window_candidate_type

  type :: flexwin_type
    integer :: n_win = 0
    real(kind=dp), allocatable :: dat(:)  ! data array for windowing
    real(kind=dp), allocatable :: syn(:)  ! synthetic array for windowing
    real(kind=dp), allocatable :: twin(:,:)   ! window time (s)
    real(kind=dp), allocatable :: wt(:)   ! window weight
    real(kind=dp), allocatable :: cc_coe(:)   ! cross-correlation coefficient
    real(kind=dp), allocatable :: time_shift(:)   ! time shift
    real(kind=dp), allocatable :: dlnA(:)     ! amplitude ratio
    real(kind=dp) :: noise_level, tstart, tend, min_period, max_period, dt, t0
    integer :: nstart, nend, npts
    integer, allocatable :: win_samp(:,:)
  contains
    procedure :: select_windows
    procedure, private :: initialize, calculate_criteria, reject_on_minima_water_level, &
                          reject_on_prominence_of_central_peak, reject_on_phase_separation, &
                          curtail_length_of_windows, reject_on_traveltimes, remove_duplicates, merge_windows
    final :: finalize
  end type flexwin_type

  interface flexwin_type
     module procedure new_flexwin
  end interface

contains
  function new_flexwin(dat, syn, dt, t0, tp, dis, min_period, max_period) result(this)
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt, t0, tp, dis, min_period, max_period
    type(flexwin_type) :: this
    call this%initialize(dat, syn, dt, t0, tp, dis, min_period, max_period)
  end function new_flexwin

  subroutine finalize(this)
    type(flexwin_type), intent(inout) :: this
    if (allocated(this%dat)) deallocate(this%dat)
    if (allocated(this%syn)) deallocate(this%syn)
    if (allocated(this%twin)) deallocate(this%twin)
    if (allocated(this%wt)) deallocate(this%wt)
    if (allocated(this%cc_coe)) deallocate(this%cc_coe)
    if (allocated(this%time_shift)) deallocate(this%time_shift)
    if (allocated(this%dlnA)) deallocate(this%dlnA)
    if (allocated(this%win_samp)) deallocate(this%win_samp)
  end subroutine finalize

  subroutine initialize(this, dat, syn, dt, t0, tp, dis, min_period, max_period)
    class(flexwin_type), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt, t0, tp, dis, min_period, max_period

    this%dat = dat / maxval(abs(syn))
    this%syn = syn / maxval(abs(syn))
    this%dt = dt
    this%t0 = t0 ! t0 is the same as specfem t0
    this%npts = size(dat)
    this%min_period = min_period
    this%max_period = max_period
    ! Match Pyflex calculate_preliminiaries limits
    ! Store tstart and tend relative to trace start (subtract t0)
    this%tstart = tp - win_config_global%max_time_before_first_arrival + t0
    this%tend = dis / win_config_global%min_velocity + max_period + t0
    this%nstart = int(this%tstart / dt) + 1
    this%nend = min(size(dat), int(this%tend / dt) + 1)
    ! Use Max Amplitude for noise level to match Pyflex default
    ! In initialize: noise_level = sqrt(sum(dat(1:nstart)**2) / nstart)
    ! So noise_level^2 is mean energy of noise.
    ! Energy of window = sum(dat_win**2) / nlen
    if (this%nstart > 0) then
       this%noise_level = maxval(abs(dat(1:this%nstart)))
    else
       this%noise_level = 1.0_dp ! Avoid division by zero
    end if

  end subroutine initialize

  subroutine calculate_criteria(this, cand)
    class(flexwin_type), intent(in) :: this
    type(window_candidate_type), intent(inout) :: cand
    
    real(kind=dp), allocatable :: dat_win(:), syn_win(:)
    integer :: nlen
    real(kind=dp) :: sigma_dt, sigma_dlna, dt_sigma_min, dlna_sigma_min
    
    nlen = cand%right_idx - cand%left_idx + 1
    allocate(dat_win(nlen), syn_win(nlen))
    
    dat_win = this%dat(cand%left_idx:cand%right_idx)
    syn_win = this%syn(cand%left_idx:cand%right_idx)
    
    ! Taper
    ! Pyflex does NOT taper windows before CC calculation.
    ! call window_taper(dat_win, 0.1_dp, 1) ! Using 0.1 taper fraction as default
    ! call window_taper(syn_win, 0.1_dp, 1)
    
    ! Calculate CC and shift
    ! Using calc_cc_shift from cross_correlate
    ! We need dummy values for sigmas if we just want shift and dlna
    dt_sigma_min = 1.0_dp ! Dummy
    dlna_sigma_min = 1.0_dp ! Dummy
    
    call calc_cc_shift(dat_win, syn_win, this%dt, dt_sigma_min, dlna_sigma_min, &
                       cand%dt_shift, cand%dlnA, sigma_dt, sigma_dlna)
                       
    ! Calculate max CC coefficient
    cand%cc = cc_max_coef(dat_win, syn_win)
    
    ! Calculate SNR (Amplitude and Energy)
    ! Assuming noise level is pre-calculated in this%noise_level (amplitude)
    if (this%noise_level > 0.0_dp) then
       cand%snr_amp = maxval(abs(dat_win)) / this%noise_level
    else
       cand%snr_amp = 1.0e5_dp ! High value if no noise
    end if
    
    ! Energy SNR
    ! Need noise energy. this%noise_level is RMS amplitude?
    ! In initialize: noise_level = sqrt(sum(dat(1:nstart)**2) / nstart)
    ! So noise_level^2 is mean energy of noise.
    ! Energy of window = sum(dat_win**2) / nlen
    if (this%noise_level > 0.0_dp) then
       cand%snr_eng = (sum(dat_win**2) / real(nlen, kind=dp)) / (this%noise_level**2)
    else
       cand%snr_eng = 1.0e5_dp
    end if
    
    ! Calculate weight (example: based on CC and shift)
    ! Python: weight_function = lambda win: win.max_cc_value * win.energy_ratio * win.len_seconds
    ! But here we can use a simpler one or configurable one.
    ! For now, let's use CC * length
    ! cand%weight = cand%cc * abs(cand%dlnA) ! Just a placeholder, need to check config
    ! Actually pyflex default weight function is often just length or something similar if not specified.
    ! Let's use length * CC for now.
    ! Use (right - left) instead of nlen (right - left + 1) to match Pyflex and avoid splitting bias
    ! Also divide by min_period to match Pyflex exactly
    cand%weight = real(cand%right_idx - cand%left_idx, kind=dp) * this%dt / this%min_period * cand%cc
    ! cand%weight = real(cand%right_idx - cand%left_idx, kind=dp)

  end subroutine calculate_criteria

  subroutine select_windows(this)
    class(flexwin_type), intent(inout) :: this
    
    real(kind=dp), allocatable :: env_data(:), stalta(:)
    complex(kind=dp), allocatable :: hilb(:)
    type(fft_cls) :: fft_ins
    integer, allocatable :: peaks(:), troughs(:)
    integer :: n_peaks, n_troughs, i, j, k, n_cand
    type(window_candidate_type), allocatable :: candidates(:)
    type(interval_type), allocatable :: intervals(:)
    integer, allocatable :: selected_indices(:)
    integer :: n_selected
    
    ! 1. Calculate STA/LTA
    hilb = fft_ins%hilbert(this%syn)
    env_data = abs(hilb)
    
    stalta = STA_LTA(env_data, this%dt, this%min_period)
    
    deallocate(hilb, env_data)
    
    ! 2. Find peaks and troughs
    peaks = find_maxima(stalta)
    troughs = find_maxima(-stalta) ! Minima are maxima of negative
    
    n_peaks = size(peaks)
    n_troughs = size(troughs)
    
    if (n_peaks == 0 .or. n_troughs < 2) return
    
    ! 3. Generate candidates
    allocate(candidates(max_good_win))
    n_cand = 0
    
    do i = 1, n_peaks
       if (stalta(peaks(i)) < win_config_global%stalta_water_level) cycle
       
       ! Iterate over ALL troughs to the left
       do j = 1, n_troughs
          if (troughs(j) >= peaks(i)) exit ! No more left troughs
          
          ! Iterate over ALL troughs to the right
          do k = 1, n_troughs
             if (troughs(k) <= peaks(i)) cycle
             
             ! Create candidate
             if (n_cand >= max_good_win) exit
             
             n_cand = n_cand + 1
             candidates(n_cand)%left_idx = troughs(j)
             candidates(n_cand)%right_idx = troughs(k)
             candidates(n_cand)%peak_idx = peaks(i)
             candidates(n_cand)%rejected = .false.
             
             ! Initial rejection based on travel times (if available)
             ! In Fortran, tstart/tend are set in initialize.
             if ((real(candidates(n_cand)%left_idx - 1, kind=dp) * this%dt < this%tstart) .or. &
                 (real(candidates(n_cand)%right_idx - 1, kind=dp) * this%dt > this%tend)) then
                candidates(n_cand)%rejected = .true.
             end if
             
          end do
          if (n_cand >= max_good_win) exit
       end do
       if (n_cand >= max_good_win) exit
    end do
    
    ! 4. Apply Rejection Steps (Order matters!)
    
    ! 4.0 Reject on travel times (Pyflex does this after initial selection)
    call this%reject_on_traveltimes(candidates, n_cand, 0.0_dp, 0.0_dp) ! tp and dis are not needed as we use stored tstart/tend/max_period
    
    ! 4.1 Reject on minimum length (first pass)
    do i = 1, n_cand
       if (candidates(i)%rejected) cycle
       if ((candidates(i)%right_idx - candidates(i)%left_idx) * this%dt < &
            win_config_global%min_win_len_fac * this%min_period) then
          candidates(i)%rejected = .true.
       end if
    end do

    ! 4.2 Reject on minima water level
    call this%reject_on_minima_water_level(candidates, n_cand, stalta, troughs)

    ! 4.3 Reject on prominence of central peak
    call this%reject_on_prominence_of_central_peak(candidates, n_cand, stalta, troughs)

    ! 4.4 Reject on phase separation
    call this%reject_on_phase_separation(candidates, n_cand, stalta, peaks, troughs)

    ! 4.5 Curtail length of windows
    call this%curtail_length_of_windows(candidates, n_cand, peaks)
    
    ! 4.6 Remove duplicates (and recenter)
    call this%remove_duplicates(candidates, n_cand)
    
    ! 4.7 Reject on minimum length (second pass)
    do i = 1, n_cand
       if (candidates(i)%rejected) cycle
       if ((candidates(i)%right_idx - candidates(i)%left_idx) * this%dt < &
            win_config_global%min_win_len_fac * this%min_period) then
          candidates(i)%rejected = .true.
       end if
    end do

    ! 4.8 Calculate Criteria and Reject based on SNR and Data Fit
    do i = 1, n_cand
      if (candidates(i)%rejected) cycle
      
      call this%calculate_criteria(candidates(i))
      
      ! Check CC
      if (candidates(i)%cc < win_config_global%threshold_corr) then
        candidates(i)%rejected = .true.
      end if
      
      ! Check Shift
      if (abs(candidates(i)%dt_shift) > win_config_global%threshold_shift_fac * this%min_period) then
        candidates(i)%rejected = .true.
      end if
      
      ! Reject on dlnA
      if (abs(candidates(i)%dlnA) > win_config_global%threshold_dlna) then
        candidates(i)%rejected = .true.
      end if
      
      ! Reject on SNR
      if (candidates(i)%snr_amp < win_config_global%min_snr_window) then
        candidates(i)%rejected = .true.
      end if
      
    end do
    
    ! 5. Resolution Strategy
    if (.not. win_config_global%is_split_phases) then
      call this%merge_windows(candidates, n_cand)
      
      ! All remaining candidates are selected
      n_selected = n_cand
      if (n_selected > 0) then
        if (allocated(selected_indices)) deallocate(selected_indices)
        allocate(selected_indices(n_selected))
        do i = 1, n_selected
          selected_indices(i) = i
        end do
      end if
       
    else
    
      ! Count valid candidates
      k = 0
      do i = 1, n_cand
        if (.not. candidates(i)%rejected) k = k + 1
      end do
      
      if (k == 0) return
      
      allocate(intervals(k))
      k = 0
      do i = 1, n_cand
        if (.not. candidates(i)%rejected) then
          k = k + 1
          intervals(k)%original_index = i
          intervals(k)%left = real(candidates(i)%left_idx, kind=dp) * this%dt
          intervals(k)%right = real(candidates(i)%right_idx, kind=dp) * this%dt
          intervals(k)%weight = candidates(i)%weight
        end if
      end do
      
      ! 6. Schedule
      call schedule_weighted_intervals(intervals, selected_indices, n_selected)
      deallocate(intervals)
    
    end if
    
    ! 7. Store results
    this%n_win = n_selected
    if (allocated(this%twin)) deallocate(this%twin)
    if (allocated(this%win_samp)) deallocate(this%win_samp)
    if (allocated(this%wt)) deallocate(this%wt)
    if (allocated(this%cc_coe)) deallocate(this%cc_coe)
    if (allocated(this%time_shift)) deallocate(this%time_shift)
    if (allocated(this%dlnA)) deallocate(this%dlnA)
    
    allocate(this%twin(n_selected, 2))
    allocate(this%win_samp(n_selected, 2))
    allocate(this%wt(n_selected))
    allocate(this%cc_coe(n_selected))
    allocate(this%time_shift(n_selected))
    allocate(this%dlnA(n_selected))
    
    do i = 1, n_selected
      j = selected_indices(i) ! Index in candidates array
      this%win_samp(i, 1) = candidates(j)%left_idx
      this%win_samp(i, 2) = candidates(j)%right_idx
      this%twin(i, 1) = real(candidates(j)%left_idx - 1, kind=dp) * this%dt - this%t0
      this%twin(i, 2) = real(candidates(j)%right_idx - 1, kind=dp) * this%dt - this%t0
      this%wt(i) = candidates(j)%weight
      this%cc_coe(i) = candidates(j)%cc
      this%time_shift(i) = candidates(j)%dt_shift
      this%dlnA(i) = candidates(j)%dlnA
    end do
    
  end subroutine select_windows

  subroutine reject_on_traveltimes(this, candidates, n_cand, tp, dis)
    class(flexwin_type), intent(in) :: this
    type(window_candidate_type), intent(inout) :: candidates(:)
    integer, intent(in) :: n_cand
    real(kind=dp), intent(in) :: tp, dis ! Unused, using stored values
    
    real(kind=dp) :: min_time, max_time
    real(kind=dp) :: win_start, win_end
    integer :: i
    
    min_time = this%tstart + win_config_global%max_time_before_first_arrival - this%min_period
    max_time = this%tend - this%max_period
    
    do i = 1, n_cand
      if (candidates(i)%rejected) cycle
      
      win_start = real(candidates(i)%left_idx - 1, kind=dp) * this%dt
      win_end = real(candidates(i)%right_idx - 1, kind=dp) * this%dt
      
      if (.not. (win_end >= min_time .and. win_start <= max_time)) then
        candidates(i)%rejected = .true.
      end if
    end do
  end subroutine reject_on_traveltimes

  subroutine reject_on_minima_water_level(this, candidates, n_cand, stalta, troughs)
    class(flexwin_type), intent(in) :: this
    type(window_candidate_type), intent(inout) :: candidates(:)
    integer, intent(in) :: n_cand
    real(kind=dp), intent(in) :: stalta(:)
    integer, intent(in) :: troughs(:)
    
    integer :: i, j
    real(kind=dp) :: waterlevel_midpoint
    logical :: reject
    
    do i = 1, n_cand
      if (candidates(i)%rejected) cycle
      
      waterlevel_midpoint = win_config_global%c_0 * win_config_global%stalta_water_level
      
      reject = .false.
      ! Check internal minima
      do j = 1, size(troughs)
        if (troughs(j) > candidates(i)%left_idx .and. troughs(j) < candidates(i)%right_idx) then
          if (stalta(troughs(j)) <= waterlevel_midpoint) then
            reject = .true.
            exit
          end if
        end if
      end do
      
      if (reject) candidates(i)%rejected = .true.
    end do
  end subroutine reject_on_minima_water_level

  subroutine reject_on_prominence_of_central_peak(this, candidates, n_cand, stalta, troughs)
    class(flexwin_type), intent(in) :: this
    type(window_candidate_type), intent(inout) :: candidates(:)
    integer, intent(in) :: n_cand
    real(kind=dp), intent(in) :: stalta(:)
    integer, intent(in) :: troughs(:)
    
    integer :: i, j
    integer :: left_trough_idx, right_trough_idx
    real(kind=dp) :: left_val, right_val, center_val, delta_left, delta_right
    
    if (win_config_global%c_2 == 0.0_dp) return
    
    do i = 1, n_cand
      if (candidates(i)%rejected) cycle
      
      left_trough_idx = -1
      right_trough_idx = -1
      
      ! Find left trough (largest index < peak_idx)
      do j = size(troughs), 1, -1
        if (troughs(j) < candidates(i)%peak_idx) then
          left_trough_idx = troughs(j)
          exit
        end if
      end do
      
      ! Find right trough (smallest index > peak_idx)
      do j = 1, size(troughs)
        if (troughs(j) > candidates(i)%peak_idx) then
          right_trough_idx = troughs(j)
          exit
        end if
      end do
      
      if (left_trough_idx == -1 .or. right_trough_idx == -1) then
        candidates(i)%rejected = .true.
        cycle
      end if
      
      left_val = stalta(left_trough_idx)
      right_val = stalta(right_trough_idx)
      center_val = stalta(candidates(i)%peak_idx)
      
      delta_left = center_val - left_val
      delta_right = center_val - right_val
      
      if (delta_left < win_config_global%c_2 * center_val .or. &
          delta_right < win_config_global%c_2 * center_val) then
        candidates(i)%rejected = .true.
      end if
      
    end do
  end subroutine reject_on_prominence_of_central_peak

  subroutine reject_on_phase_separation(this, candidates, n_cand, stalta, peaks, troughs)
    class(flexwin_type), intent(in) :: this
    type(window_candidate_type), intent(inout) :: candidates(:)
    integer, intent(in) :: n_cand
    real(kind=dp), intent(in) :: stalta(:)
    integer, intent(in) :: peaks(:), troughs(:)
    
    integer :: i, j
    real(kind=dp) :: stalta_min, d_stalta_center, d_stalta, d_time, f_time
    integer :: max_idx
    logical :: reject
    
    do i = 1, n_cand
      if (candidates(i)%rejected) cycle
      
      ! Find lowest minimum within window (inclusive of boundaries)
      stalta_min = huge(stalta_min)
      do j = 1, size(troughs)
        if (troughs(j) >= candidates(i)%left_idx .and. troughs(j) <= candidates(i)%right_idx) then
          if (stalta(troughs(j)) < stalta_min) stalta_min = stalta(troughs(j))
        end if
      end do
      
      d_stalta_center = stalta(candidates(i)%peak_idx) - stalta_min
      
      ! Find all internal maxima (excluding center)
      reject = .false.
      do j = 1, size(peaks)
        if (peaks(j) >= candidates(i)%left_idx .and. peaks(j) <= candidates(i)%right_idx .and. &
          peaks(j) /= candidates(i)%peak_idx) then
          
          max_idx = peaks(j)
          d_stalta = stalta(max_idx) - stalta_min
          
          d_time = abs(candidates(i)%peak_idx - max_idx) * this%dt / this%min_period
          
          if (d_time >= win_config_global%c_3b) then
            f_time = exp(-((d_time - win_config_global%c_3b) / win_config_global%c_3b)**2)
          else
            f_time = 1.0_dp
          end if
          
          if (d_stalta > (win_config_global%c_3a * d_stalta_center * f_time)) then
            reject = .true.
            exit
          end if
        end if
      end do
      
      if (reject) candidates(i)%rejected = .true.
    end do
  end subroutine reject_on_phase_separation

  subroutine curtail_length_of_windows(this, candidates, n_cand, peaks)
    class(flexwin_type), intent(in) :: this
    type(window_candidate_type), intent(inout) :: candidates(:)
    integer, intent(in) :: n_cand
    integer, intent(in) :: peaks(:)
    
    integer :: i, j
    real(kind=dp) :: time_decay_left, time_decay_right
    integer :: i_left, i_right
    real(kind=dp) :: delta_left, delta_right
    integer, allocatable :: internal_maxima(:)
    integer :: n_internal
    
    time_decay_left = this%min_period * win_config_global%c_4a / this%dt
    time_decay_right = this%min_period * win_config_global%c_4b / this%dt
    
    do i = 1, n_cand
       if (candidates(i)%rejected) cycle
       
       ! Find internal maxima (excluding center)
       n_internal = 0
       do j = 1, size(peaks)
          if (peaks(j) >= candidates(i)%left_idx .and. peaks(j) <= candidates(i)%right_idx .and. &
              peaks(j) /= candidates(i)%peak_idx) then
             n_internal = n_internal + 1
          end if
       end do
       
       if (n_internal < 2) cycle
       
       allocate(internal_maxima(n_internal))
       n_internal = 0
       do j = 1, size(peaks)
          if (peaks(j) >= candidates(i)%left_idx .and. peaks(j) <= candidates(i)%right_idx .and. &
              peaks(j) /= candidates(i)%peak_idx) then
             n_internal = n_internal + 1
             internal_maxima(n_internal) = peaks(j)
          end if
       end do
       
       i_left = internal_maxima(1)
       i_right = internal_maxima(n_internal)
       
       delta_left = real(i_left - candidates(i)%left_idx, kind=dp)
       delta_right = real(candidates(i)%right_idx - i_right, kind=dp)
       
       if (delta_left > time_decay_left) then
          candidates(i)%left_idx = int(real(i_left, kind=dp) - time_decay_left)
       end if
       
       if (delta_right > time_decay_right) then
          candidates(i)%right_idx = int(real(i_right, kind=dp) + time_decay_right)
       end if
       
       deallocate(internal_maxima)
    end do
  end subroutine curtail_length_of_windows

  subroutine remove_duplicates(this, candidates, n_cand)
    class(flexwin_type), intent(in) :: this
    type(window_candidate_type), intent(inout) :: candidates(:)
    integer, intent(in) :: n_cand
    
    integer :: i, j
    
    do i = 1, n_cand
       if (candidates(i)%rejected) cycle
       
       ! Recenter (match Pyflex)
       candidates(i)%peak_idx = (candidates(i)%left_idx + candidates(i)%right_idx) / 2
       
       do j = i + 1, n_cand
          if (candidates(j)%rejected) cycle
          
          if (candidates(j)%left_idx == candidates(i)%left_idx .and. &
              candidates(j)%right_idx == candidates(i)%right_idx) then
             candidates(j)%rejected = .true.
          end if
       end do
    end do
  end subroutine remove_duplicates

  subroutine merge_windows(this, candidates, n_cand)
    class(flexwin_type), intent(inout) :: this
    type(window_candidate_type), intent(inout) :: candidates(:)
    integer, intent(inout) :: n_cand
    
    type(window_candidate_type), allocatable :: valid_candidates(:), merged_candidates(:)
    integer :: i, k, n_valid, n_merged
    type(window_candidate_type) :: right_win
    
    ! 1. Collect valid candidates
    n_valid = 0
    do i = 1, n_cand
      if (.not. candidates(i)%rejected) n_valid = n_valid + 1
    end do
    
    if (n_valid == 0) then
      n_cand = 0
      return
    end if
    
    allocate(valid_candidates(n_valid))
    k = 0
    do i = 1, n_cand
      if (.not. candidates(i)%rejected) then
        k = k + 1
        valid_candidates(k) = candidates(i)
      end if
    end do
    
    ! 2. Sort by left_idx
    call sort_candidates_by_left(valid_candidates, n_valid)
    
    ! 3. Merge
    allocate(merged_candidates(n_valid)) ! Max possible size
    n_merged = 0
    
    if (n_valid > 0) then
      n_merged = 1
      merged_candidates(1) = valid_candidates(1)
      
      do i = 2, n_valid
        right_win = valid_candidates(i)
        ! Check overlap
        ! Python: if (left_win.right + 1) < right_win.left:
        if (merged_candidates(n_merged)%right_idx + 1 < right_win%left_idx) then
            ! No overlap, add new window
            n_merged = n_merged + 1
            merged_candidates(n_merged) = right_win
        else
            ! Overlap, merge
            merged_candidates(n_merged)%right_idx = max(merged_candidates(n_merged)%right_idx, right_win%right_idx)
        end if
      end do
    end if
    
    ! 4. Recalculate criteria and update candidates
    n_cand = n_merged
    do i = 1, n_merged
      ! Recenter
      merged_candidates(i)%peak_idx = (merged_candidates(i)%left_idx + merged_candidates(i)%right_idx) / 2
      
      ! Recalculate criteria
      call this%calculate_criteria(merged_candidates(i))
      
      ! Copy back to candidates
      candidates(i) = merged_candidates(i)
      candidates(i)%rejected = .false.
    end do
    
    ! Mark the rest as rejected (though n_cand is updated, caller might iterate up to old n_cand)
    ! Actually, caller usually uses n_cand. But let's be safe.
    ! Since candidates is intent(inout) and size is max_good_win, we just update n_cand.
    
    deallocate(valid_candidates, merged_candidates)

  end subroutine merge_windows

  subroutine sort_candidates_by_left(candidates, n)
    type(window_candidate_type), intent(inout) :: candidates(:)
    integer, intent(in) :: n
    
    integer :: i, j
    type(window_candidate_type) :: temp
    
    ! Simple bubble sort for small n
    do i = 1, n-1
       do j = i+1, n
          if (candidates(j)%left_idx < candidates(i)%left_idx) then
             temp = candidates(i)
             candidates(i) = candidates(j)
             candidates(j) = temp
          end if
       end do
    end do
  end subroutine sort_candidates_by_left

end module flexwin