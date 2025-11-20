module win_sel
  use config
  use signal
  use cross_correlate
  use sta_lta_mod
  use fftpack
  use adj_config, only: win_config_global
  implicit none

  integer, private, parameter :: max_good_win = 1000

  type :: win_sel_type
    integer :: n_win = 0
    real(kind=dp), allocatable :: dat(:)  ! data array for windowing
    real(kind=dp), allocatable :: syn(:)  ! synthetic array for windowing
    real(kind=dp), allocatable :: twin(:,:)   ! window time (s)
    real(kind=dp), allocatable :: wt(:)   ! window weight
    real(kind=dp), allocatable :: cc_coe(:)   ! cross-correlation coefficient
    real(kind=dp), allocatable :: time_shift(:)   ! time shift
    real(kind=dp), allocatable :: times_cc(:) ! logical array for good windows
    real(kind=dp) :: noise_level, tstart, tend, min_period, jump_buffer, dt, t0
    integer :: nstart, nend, npts
    integer, allocatable :: win_samp(:,:)
  contains
    procedure :: init => initialize, gen_good_windows
    procedure, private :: sliding_cc, split_window_by_phases
  end type win_sel_type

contains
  subroutine initialize(this, dat, syn, dt, t0, tp, dis, min_period)
    class(win_sel_type), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt, t0, tp, dis, min_period

    if (allocated(this%dat)) deallocate(this%dat)
    if (allocated(this%syn)) deallocate(this%syn)
    if (allocated(this%time_shift)) deallocate(this%time_shift)
    if (allocated(this%cc_coe)) deallocate(this%cc_coe)
    if (allocated(this%times_cc)) deallocate(this%times_cc)
    if (allocated(this%twin)) deallocate(this%twin)
    if (allocated(this%wt)) deallocate(this%wt)
    if (allocated(this%win_samp)) deallocate(this%win_samp)
    this%dat = dat / maxval(abs(syn))
    this%syn = syn / maxval(abs(syn))
    this%dt = dt
    this%t0 = t0 ! t0 is the same as specfem t0
    this%npts = size(dat)
    this%min_period = min_period
    this%tstart = tp - min_period * 1.5 + t0
    this%tend = dis / win_config_global%min_velocity + min_period * 1.5 + t0
    this%jump_buffer = win_config_global%jump_fac * this%min_period
    this%nstart = int(this%tstart / dt) + 1
    this%nend = min(size(dat), int(this%tend / dt) + 1)
    this%noise_level = sqrt(sum(dat(1:this%nstart)**2) / real(this%nstart, kind=dp))

  end subroutine initialize

  subroutine sliding_cc(this)
    class(win_sel_type), intent(inout) :: this
    integer :: nlen, i, nb, ne
    real(kind=dp), allocatable :: dat_win(:), syn_win(:)

    nlen = int(win_config_global%sliding_win_len_fac * this%min_period / this%dt) + 1
    allocate(dat_win(nlen), syn_win(nlen))
    allocate(this%time_shift(this%npts - nlen))
    allocate(this%cc_coe(this%npts - nlen))
    allocate(this%times_cc(this%npts - nlen))
    this%time_shift = 0.0_dp
    this%cc_coe = 0.0_dp
    this%times_cc = 0.0_dp
    do i = 1, this%npts - nlen
      nb = i
      ne = i + nlen - 1
      dat_win(:) = this%dat(nb:ne)
      syn_win(:) = this%syn(nb:ne)
      call window_taper(dat_win, 1.0_dp, 1)
      call window_taper(syn_win, 1.0_dp, 1)
      call xcorr_shift_coe(dat_win, syn_win, this%dt, this%time_shift(i), this%cc_coe(i))
      this%times_cc(i) = (dble(i - 1) + dble(nlen) / 2.0_dp) * this%dt
    enddo
    
  end subroutine sliding_cc

  subroutine gen_good_windows(this)
    class(win_sel_type), intent(inout) :: this
    integer :: i, n_jump, n_groups, group_start, group_end, ib, ie, n_peaks, &
               groups(max_good_win, 2)
    integer, allocatable, dimension(:) :: max_peaks, min_peaks
    logical, allocatable, dimension(:) :: good_windows, good_cc, good_shift, good_time
    real(kind=dp), allocatable :: time_shift_grad(:)
    logical :: in_group, good_groups(max_good_win)
    real(kind=dp) :: max_amp, eng_dat, eng_syn

    ! Variables for STA/LTA and splitting within gen_good_windows
    integer, allocatable :: peaks(:), troughs(:)
    real(kind=dp), allocatable :: stalta(:)

    ! Final window containers (collect passed sub-windows)
    integer, allocatable :: final_win_samp(:,:)
    real(kind=dp), allocatable :: final_twin(:,:)
    integer :: final_count, k, n_sub_windows
    integer, allocatable :: sub_windows(:,:)

    ! First call sliding_cc to compute time_shift and cc_coe
    call this%sliding_cc()

    ! Compute envelope + STA/LTA once and find peaks/troughs for possible splitting
    block
      real(kind=dp), allocatable :: env_data(:)
      complex(kind=dp), allocatable :: hilb(:)
      type(fft_cls) :: fft_ins
      
      if (win_config_global%is_split_phases) then
        hilb = fft_ins%hilbert(this%syn)
        env_data = abs(hilb)
        stalta = STA_LTA(env_data, this%dt, this%min_period)
        peaks = find_maxima(stalta)
        troughs = find_maxima(-stalta)
        deallocate(hilb, env_data)
      end if
    end block

    ! Allocate logical arrays
    allocate(good_cc(size(this%cc_coe)))
    allocate(good_shift(size(this%time_shift)))
    allocate(good_time(size(this%times_cc)))
    allocate(good_windows(size(this%time_shift)))
    
    good_cc = (this%cc_coe > win_config_global%threshold_corr)
    good_shift = (abs(this%time_shift) < win_config_global%threshold_shift_fac * this%min_period)
    good_time = (this%times_cc >= this%tstart) .and. (this%times_cc <= this%tend)
    good_windows = good_cc .and. good_shift .and. good_time

    allocate(time_shift_grad(size(this%time_shift)))
    time_shift_grad = 0.0_dp
    do concurrent (i = 2:size(this%time_shift))
      time_shift_grad(i) = this%time_shift(i) - this%time_shift(i-1)
    end do

    ! Mark windows around large jumps as bad
    do i = 1, size(this%time_shift)
      if (abs(time_shift_grad(i)) > this%jump_buffer) then
        n_jump = nint(this%jump_buffer / this%dt)
        ib = max(1, i - n_jump)
        ie = min(size(this%time_shift), i + n_jump)
        good_windows(ib:ie) = .false.
      end if
    end do

    ! Identify continuous good window segments
    n_groups = 0
    in_group = .false.
    good_groups = .false.
    group_start = 0
    do i = 1, size(this%time_shift)
      if (good_windows(i) .and. .not. in_group) then
        in_group = .true.
        group_start = i
      else if (.not. good_windows(i) .and. in_group) then
        in_group = .false.
        group_end = i - 1
        n_groups = n_groups + 1
        good_groups(n_groups) = .true.
        groups(n_groups, 1) = group_start
        groups(n_groups, 2) = group_end
      end if
    end do

    if (in_group) then
      group_end = size(good_windows)
      n_groups = n_groups + 1
      good_groups(n_groups) = .true.
      groups(n_groups, 1) = group_start
      groups(n_groups, 2) = group_end
    end if

    ! Prepare final containers (worst case: one output window per input sample group)
    final_count = 0
    allocate(final_win_samp(max_good_win, 2))
    allocate(final_twin(max_good_win, 2))

    ! Process each group: split first (if enabled), then apply all quality checks
    do i = 1, n_groups
      ! First split into phases if enabled, otherwise use original group as single window
      if (win_config_global%is_split_phases) then
        call this%split_window_by_phases(groups(i,1), groups(i,2), stalta, peaks, troughs, &
                                         int(win_config_global%min_win_len_fac * this%min_period / this%dt), &
                                         sub_windows, n_sub_windows)
      else
        ! No splitting: treat original group as single sub-window
        allocate(sub_windows(1, 2))
        sub_windows(1, :) = groups(i, :)
        n_sub_windows = 1
      end if

      ! Now apply all quality checks to each sub-window
      do k = 1, n_sub_windows
        ! Check 1: time length too short
        if ((sub_windows(k, 2) - sub_windows(k, 1)) * this%dt < &
             win_config_global%min_win_len_fac * this%min_period) then
          cycle
        end if

        ! Convert times_cc indices to waveform sample indices
        ib = int(this%times_cc(sub_windows(k, 1)) / this%dt) + 1
        ie = int(this%times_cc(sub_windows(k, 2)) / this%dt) + 1
        ib = max(1, min(ib, this%npts))
        ie = max(1, min(ie, this%npts))
        if (ie <= ib) cycle

        ! Check 2: number of peaks too few (check both observed and synthetic)
        max_peaks = find_maxima(this%dat(ib:ie))
        min_peaks = find_maxima(-this%dat(ib:ie))
        n_peaks = size(max_peaks) + size(min_peaks)
        max_peaks = find_maxima(this%syn(ib:ie))
        min_peaks = find_maxima(-this%syn(ib:ie))
        n_peaks = min(n_peaks, size(max_peaks) + size(min_peaks))
        if (n_peaks < win_config_global%min_peaks_troughs) cycle

        ! Check 3: signal to noise ratio too low
        max_amp = maxval(abs(this%dat(ib:ie)))
        if (max_amp / this%noise_level < win_config_global%min_snr_window) cycle

        ! Check 4: Energy ratio check
        eng_dat = sum(this%dat(ib:ie)**2)
        eng_syn = sum(this%syn(ib:ie)**2)
        if (eng_dat / eng_syn > win_config_global%min_energy_ratio .or. &
            eng_dat / eng_syn < 1/win_config_global%min_energy_ratio) cycle

        ! All checks passed - add to final windows
        final_count = final_count + 1
        final_win_samp(final_count, 1) = sub_windows(k, 1)
        final_win_samp(final_count, 2) = sub_windows(k, 2)
        final_twin(final_count, 1) = this%times_cc(sub_windows(k, 1))
        final_twin(final_count, 2) = this%times_cc(sub_windows(k, 2))
      end do ! k = 1, n_sub_windows

      ! Deallocate sub_windows for next iteration
      if (allocated(sub_windows)) deallocate(sub_windows)
    end do ! i = 1, n_groups

    ! Finalize results
    this%n_win = final_count
    if (final_count > 0) then
      if (allocated(this%twin)) deallocate(this%twin)
      if (allocated(this%win_samp)) deallocate(this%win_samp)
      allocate(this%twin(this%n_win, 2))
      allocate(this%win_samp(this%n_win, 2))
      this%twin = final_twin(1:this%n_win, :)
      this%win_samp = final_win_samp(1:this%n_win, :)
    end if

    ! cleanup
    if (allocated(peaks)) deallocate(peaks)
    if (allocated(troughs)) deallocate(troughs)
    if (allocated(final_win_samp)) deallocate(final_win_samp)
    if (allocated(final_twin)) deallocate(final_twin)

  end subroutine gen_good_windows

  subroutine split_window_by_phases(this, start_idx, end_idx, stalta, peaks, troughs, &
                                    min_window_length, sub_windows, n_sub_windows)
    ! Use STA/LTA to detect multiple phases within a window and split at troughs
    ! 
    ! Simplified logic:
    ! - Find all peaks within the window
    ! - If multiple peaks exist, split at the troughs between them
    ! - No significance threshold - pure phase separation based on STA/LTA structure
    class(win_sel_type), intent(inout) :: this
    integer, intent(in) :: start_idx, end_idx
    real(kind=dp), dimension(:), intent(in) :: stalta
    integer, dimension(:), intent(in) :: peaks, troughs
    integer, dimension(:, :), allocatable, intent(out) :: sub_windows
    integer, intent(out) :: n_sub_windows
    integer, intent(in) :: min_window_length
    
    integer :: start_sample, end_sample, i, j, n_window_peaks, n_between_troughs
    integer :: peak1, peak2, min_trough_idx, split_point
    integer, dimension(:), allocatable :: window_peaks, split_points, between_troughs
    real(kind=dp) :: min_trough_value
    logical, dimension(:), allocatable :: mask
        
    ! Convert times_cc indices to stalta indices (sample points)
    start_sample = int(this%times_cc(start_idx) / this%dt) + 1
    end_sample = int(this%times_cc(end_idx) / this%dt) + 1
    
    ! Ensure indices are within valid range
    start_sample = max(1, min(start_sample, size(stalta)))
    end_sample = max(1, min(end_sample, size(stalta)))
    
    ! Return original window if invalid range
    if (start_sample >= end_sample) then
      allocate(sub_windows(1, 2))
      sub_windows(1, 1) = start_idx
      sub_windows(1, 2) = end_idx
      n_sub_windows = 1
      return
    end if
    
    ! Find all peaks within the window
    allocate(mask(size(peaks)))
    mask = (peaks >= start_sample) .and. (peaks <= end_sample)
    n_window_peaks = count(mask)
    
    ! Single peak or no peak, don't split
    if (n_window_peaks <= 1) then
      allocate(sub_windows(1, 2))
      sub_windows(1, 1) = start_idx
      sub_windows(1, 2) = end_idx
      n_sub_windows = 1
      deallocate(mask)
      return
    end if
    
    ! Extract window peaks
    allocate(window_peaks(n_window_peaks))
    j = 0
    do i = 1, size(peaks)
      if (mask(i)) then
        j = j + 1
        window_peaks(j) = peaks(i)
      end if
    end do
    deallocate(mask)
    
    ! Build split points array
    allocate(split_points(n_window_peaks + 1))
    split_points(1) = start_sample  ! Start point
    
    ! Find split points between consecutive peaks
    do i = 1, n_window_peaks - 1
      peak1 = window_peaks(i)
      peak2 = window_peaks(i + 1)
      
      ! Find troughs between two peaks
      allocate(mask(size(troughs)))
      mask = (troughs > peak1) .and. (troughs < peak2)
      n_between_troughs = count(mask)
      
      if (n_between_troughs > 0) then
        ! Extract troughs between peaks
        allocate(between_troughs(n_between_troughs))
        j = 0
        do j = 1, size(troughs)
          if (mask(j)) then
            between_troughs(count(mask(1:j))) = troughs(j)
          end if
        end do
        
        ! Choose the trough with minimum STA/LTA value as split point
        min_trough_value = huge(min_trough_value)
        min_trough_idx = between_troughs(1)
        do j = 1, n_between_troughs
          if (stalta(between_troughs(j)) < min_trough_value) then
            min_trough_value = stalta(between_troughs(j))
            min_trough_idx = between_troughs(j)
          end if
        end do
        split_points(i + 1) = min_trough_idx
        deallocate(between_troughs)
      else
        ! If no trough, use midpoint between two peaks
        split_points(i + 1) = (peak1 + peak2) / 2
      end if
      deallocate(mask)
    end do
    
    split_points(n_window_peaks + 1) = end_sample  ! End point
    
    ! Build sub-windows after splitting
    ! First count valid sub-windows (checking both sample length and times_cc length)
    n_sub_windows = 0
    do i = 1, n_window_peaks
      ! Check if sample length is sufficient
      if (split_points(i + 1) - split_points(i) >= min_window_length) then
        n_sub_windows = n_sub_windows + 1
      end if
    end do
    
    ! If no valid sub-windows, return original window
    if (n_sub_windows == 0) then
      allocate(sub_windows(1, 2))
      sub_windows(1, 1) = start_idx
      sub_windows(1, 2) = end_idx
      n_sub_windows = 1
      deallocate(split_points, window_peaks)
      return
    end if
    
    ! Allocate and fill sub-windows
    allocate(sub_windows(n_sub_windows, 2))
    j = 0
    do i = 1, n_window_peaks
      ! Check if sample length is sufficient
      if (split_points(i + 1) - split_points(i) >= min_window_length) then
        ! Temporarily store the sub-window indices
        j = j + 1
        ! Convert back to times_cc indices
        ! Find the times_cc index closest to split_points(i) * dt
        sub_windows(j, 1) = start_idx
        sub_windows(j, 2) = end_idx
        do split_point = start_idx, end_idx
          if (abs(this%times_cc(split_point) - split_points(i) * this%dt) < &
              abs(this%times_cc(sub_windows(j, 1)) - split_points(i) * this%dt)) then
            sub_windows(j, 1) = split_point
          end if
          if (abs(this%times_cc(split_point) - split_points(i + 1) * this%dt) < &
              abs(this%times_cc(sub_windows(j, 2)) - split_points(i + 1) * this%dt)) then
            sub_windows(j, 2) = split_point
          end if
        end do
        
        ! Ensure within original window range
        sub_windows(j, 1) = max(start_idx, min(sub_windows(j, 1), end_idx))
        sub_windows(j, 2) = max(start_idx, min(sub_windows(j, 2), end_idx))
      end if
    end do
    
  end subroutine split_window_by_phases

  subroutine xcorr_shift_coe(d, s, dt, tshift, cc_max_coef)
    real(kind=dp), dimension(:), intent(in) :: s, d
    real(kind=dp), intent(in) :: dt
    real(kind=dp), intent(out) :: cc_max_coef, tshift
    real(kind=dp), dimension(:), allocatable :: s_shift, d_shift
    integer :: nlen, ishift, index, index_shift

    nlen = size(s)
    allocate(s_shift(nlen))
    s_shift = 0.0_dp
    d_shift = d(:)

    ishift = xcorr_shift(d, s)
    do index = 1, nlen
      index_shift = index - ishift
      if (index_shift >= 1 .and. index_shift <= nlen) then
        s_shift(index) = s(index_shift)
      end if
    end do
    tshift = real(ishift, kind=dp) * dt

    ! taper the 0 edges on d
    if (ishift > 0) then
      d_shift(1:ishift) = 0.0_dp
    else if (ishift < 0) then
      d_shift(nlen+ishift+1:nlen) = 0.0_dp
    end if

    cc_max_coef = sum(d_shift * s_shift) / sqrt(sum(d_shift**2) * sum(s_shift**2))

  end subroutine xcorr_shift_coe

  function merge_sorted_arrays(arr1, arr2) result(merged)
    integer, intent(in) :: arr1(:), arr2(:)
    integer, allocatable :: merged(:)
    integer :: n1, n2, n_total
    integer :: i1, i2, idx
    
    n1 = size(arr1)
    n2 = size(arr2)
    n_total = n1 + n2
    
    allocate(merged(n_total))
    
    i1 = 1
    i2 = 1
    idx = 1
    
    ! merge the two arrays
    do while (i1 <= n1 .and. i2 <= n2)
      if (arr1(i1) <= arr2(i2)) then
        merged(idx) = arr1(i1)
        i1 = i1 + 1
      else
        merged(idx) = arr2(i2)
        i2 = i2 + 1
      end if
      idx = idx + 1
    end do
    
    ! append remaining elements
    do while (i1 <= n1)
      merged(idx) = arr1(i1)
      i1 = i1 + 1
      idx = idx + 1
    end do
    
    do while (i2 <= n2)
      merged(idx) = arr2(i2)
      i2 = i2 + 1
      idx = idx + 1
    end do
    
  end function merge_sorted_arrays

end module win_sel