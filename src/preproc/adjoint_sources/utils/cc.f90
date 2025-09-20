module cross_correlate
  use config
  use signal
  use, intrinsic :: ieee_arithmetic
  use utils, only: simpson, gradient
  implicit none

contains

  subroutine calc_cc_shift(d, s, dt, dt_sigma_min, dlna_sigma_min, &
                           tshift, dlna, sigma_dt, sigma_dlna)
    ! Calculate the cross-correlation shift between two signals

    real(kind=dp), dimension(:), intent(in) :: s, d
    real(kind=dp), intent(in) :: dt_sigma_min, dlna_sigma_min, dt
    real(kind=dp), intent(out) :: tshift, dlna, sigma_dt, sigma_dlna
    integer :: ishift

    ishift = xcorr_shift(d, s)
    tshift = ishift * dt ! Convert index shift to time shift

    ! Amplitude error
    dlna = 0.5_dp * log(sum(d(:)**2) / sum(s(:)**2))

    ! Calculate the cross-correlation error
    call calc_cc_error(d, s, dt, ishift, dlna, dt_sigma_min, dlna_sigma_min, &
                       sigma_dt, sigma_dlna)

  end subroutine calc_cc_shift

  subroutine calc_cc_adjsrc(s, tshift, dlna, dt, sigma_dt, sigma_dlna, &
                            misfit_p, misfit_q, adj_p, adj_q)
    ! Calculate the adjoint source for the cross-correlation traveltime

    real(kind=dp), dimension(:), intent(in) :: s
    real(kind=dp), intent(in) :: tshift, dlna, dt, sigma_dt, sigma_dlna
    real(kind=dp), intent(out) :: misfit_p, misfit_q
    real(kind=dp), dimension(:), allocatable, intent(out) :: adj_p, adj_q
    real(kind=dp), dimension(:), allocatable :: dsdt
    real(kind=dp) :: nnorm
    integer :: nt

    nt = size(s)

    allocate(adj_p(nt), adj_q(nt))
    adj_p = 0.0_dp
    adj_q = 0.0_dp

    ! Calculate misfits for time shift and amplitude anomaly
    misfit_p = 0.5_dp * (tshift / sigma_dt)**2
    misfit_q = 0.5_dp * (dlna / sigma_dlna)**2

    ! Calculate the adjoint source for time shift
    dsdt = gradient(s, dt)
    nnorm = simpson(dsdt**2, dt)
    adj_p(1:nt) = dsdt(1:nt) * tshift / nnorm / sigma_dt**2

    nnorm = simpson(s**2, dt)
    adj_q(1:nt) = -1.0_dp * s(1:nt) * dlna / nnorm / sigma_dlna**2

  end subroutine calc_cc_adjsrc

  subroutine calc_cc_error(d, s, dt, ishift, dlna, dt_sigma_min, dlna_sigma_min, &
                           sigma_dt, sigma_dlna)
    ! Calculate the cross-correlation error between two signals

    real(kind=dp), dimension(:), intent(in) :: s, d
    integer, intent(in) :: ishift
    real(kind=dp), intent(in) :: dt, dlna, dt_sigma_min, dlna_sigma_min
    real(kind=dp), intent(out) :: sigma_dt, sigma_dlna
    real(kind=dp), dimension(:), allocatable :: s_cc_dt, s_cc_dtdlna, s_cc_vel
    real(kind=dp) :: sigma_dt_num, sigma_dt_den, sigma_dlna_num, sigma_dlna_den

    ! Apply a scaling and time shift to the synthetic data
    call cc_correction(s, ishift, dlna, s_cc_dt, s_cc_dtdlna)

    s_cc_vel = gradient(s_cc_dtdlna, dt)

    sigma_dt_num = sum((d - s_cc_dtdlna)**2)
    sigma_dt_den = sum(s_cc_vel**2)
    sigma_dt = sqrt(sigma_dt_num / sigma_dt_den)

    sigma_dlna_num = sigma_dt_num
    sigma_dlna_den = sum(s_cc_dt**2)
    sigma_dlna = sqrt(sigma_dlna_num / sigma_dlna_den)

    ! Check if the calculated errors are within the specified limits
    if (sigma_dt < dt_sigma_min .or. ieee_is_nan(sigma_dt)) then
      sigma_dt = dt_sigma_min
    end if

    if (sigma_dlna < dlna_sigma_min .or. ieee_is_nan(sigma_dlna)) then
      sigma_dlna = dlna_sigma_min
    end if

  end subroutine calc_cc_error

  subroutine cc_correction(s, ishift, dlna, s_cc_dt, s_cc_dtdlna) ! reconstruct_syn_cc
    real(kind=dp), dimension(:), intent(in) :: s
    real(kind=dp), intent(in) :: dlna
    integer, intent(in) :: ishift
    real(kind=dp), dimension(:), allocatable, intent(out) :: s_cc_dt, s_cc_dtdlna
    integer :: index, index_shift, nlen_t

    ! Initialize output arrays
    nlen_t = size(s)
    allocate(s_cc_dt(nlen_t))
    allocate(s_cc_dtdlna(nlen_t))
    s_cc_dt = 0.0_dp
    s_cc_dtdlna = 0.0_dp

    do index = 1, nlen_t
      index_shift = index - ishift
      if (index_shift >= 1 .and. index_shift <= nlen_t) then
        s_cc_dt(index) = s(index_shift)
        s_cc_dtdlna(index) = exp(dlna) * s(index_shift)
      end if
    end do

  end subroutine cc_correction

  function xcorr_shift(d, s)
    real(kind=dp), dimension(:), intent(in) :: s, d
    real(kind=dp), dimension(:), allocatable :: cc
    integer :: xcorr_shift

    call mycorrelation_dp(d, s, cc, 1)
    xcorr_shift = maxloc(cc, dim=1) - size(s)

  end function xcorr_shift

  function cc_max_coef(d, s)
    real(kind=dp), dimension(:), intent(in) :: s, d
    real(kind=dp) :: cc_max_coef
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

    ! taper the 0 edges on d
    if (ishift > 0) then
      d_shift(1:ishift) = 0.0_dp
    else if (ishift < 0) then
      d_shift(nlen+ishift+1:nlen) = 0.0_dp
    end if

    cc_max_coef = sum(d_shift * s_shift) / sqrt(sum(d_shift**2) * sum(s_shift**2))

  end function cc_max_coef


end module cross_correlate