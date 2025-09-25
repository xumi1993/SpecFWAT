module exponentiated_phase_misfit
  use config
  use signal
  use adj_config, cfg => adj_config_global
  use fftpack
  use utils, only: simpson
  implicit none

  type(fft_cls), private :: fft_ins

  type, extends(AdjointMeasurement) :: ExponentiatedPhaseMisfit
    real(kind=dp), dimension(:), allocatable :: misfits_real, misfits_imag, &
                                                residuals_real, residuals_imag
  contains
    procedure :: calc_adjoint_source
    procedure, private :: initialize
  end type ExponentiatedPhaseMisfit

contains

  subroutine calc_adjoint_source(this, dat, syn, dt, windows)
    class(ExponentiatedPhaseMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(in) :: windows
    real(kind=dp), dimension(:), allocatable :: s, d, env_s, env_d, hilb_s, hilb_d,&
                                              env_s_wl, env_d_wl, diff_real, diff_imag, &
                                              env_s_wtr_cubic, adj_real, adj_imag
    complex(kind=dp), dimension(:), allocatable :: h_s, h_d
    real(kind=dp) :: threshold_s, threshold_d
    integer :: iwin, nlen, nb, ne, nlen_win

    ! num of measurements
    if (size(windows, 2) /= 2) then
      write(*,*) 'Error: windows must have two columns (start and end times)'
      error stop
    end if
    this%nwin = size(windows, 1)
    ! allocate measurement arrays
    call this%initialize()
  
    ! number of time samples
    nlen = size(dat)
    if (allocated(this%adj_src)) deallocate(this%adj_src)
    allocate(this%adj_src(nlen))
    this%adj_src = 0.0_dp

    ! loop over windows
    do iwin = 1, this%nwin
      call get_window_info(windows(iwin,:), dt, nb, ne, nlen_win)
      s = syn(nb:ne)
      d = dat(nb:ne)

      ! taper the windows
      call window_taper(s, cfg%taper_percentage, cfg%itaper_type)
      call window_taper(d, cfg%taper_percentage, cfg%itaper_type)

      ! calculate envelopes and hilbert transforms
      h_s = fft_ins%hilbert(s)
      h_d = fft_ins%hilbert(d)
      env_s = abs(h_s)
      env_d = abs(h_d)
      hilb_s = aimag(h_s)
      hilb_d = aimag(h_d)

      ! calculate water level
      threshold_s = cfg%wtr_env * maxval(env_s)
      threshold_d = cfg%wtr_env * maxval(env_d)
      env_s_wl = env_s + threshold_s
      env_d_wl = env_d + threshold_d

      ! calculate residuals
      diff_real = d/env_d_wl - s/env_s_wl
      diff_imag = hilb_d/env_d_wl - hilb_s/env_s_wl

      this%residuals_real(iwin) = sum(diff_real)/nlen_win
      this%residuals_imag(iwin) = sum(diff_imag)/nlen_win
      this%residuals(iwin) = this%residuals_real(iwin) + this%residuals_imag(iwin)

      this%misfits_real(iwin) = 0.5_dp * simpson(diff_real**2, dt)
      this%misfits_imag(iwin) = 0.5_dp * simpson(diff_imag**2, dt)
      this%misfits(iwin) = this%misfits_real(iwin) + this%misfits_imag(iwin)

      this%total_misfit = this%total_misfit + this%misfits(iwin)

      env_s_wtr_cubic = env_s_wl**3
      adj_real = -1.0_dp * (diff_real * hilb_s ** 2 / env_s_wtr_cubic) - &
                  aimag(fft_ins%hilbert(diff_real * s * hilb_s / env_s_wtr_cubic))
      adj_imag = diff_imag * s * hilb_s / env_s_wtr_cubic + &
                  aimag(fft_ins%hilbert(diff_imag * s ** 2 / env_s_wtr_cubic))

      ! taper the adjoint source
      call window_taper(adj_real, cfg%taper_percentage, cfg%itaper_type)
      call window_taper(adj_imag, cfg%taper_percentage, cfg%itaper_type)

      this%adj_src(nb:ne) = this%adj_src(nb:ne) + adj_real(1:nlen_win) + adj_imag(1:nlen_win)

      deallocate(s, d, h_d, h_s, env_s, env_d, hilb_s, hilb_d)
      deallocate(env_s_wl, env_d_wl, env_s_wtr_cubic, diff_real, diff_imag)
      deallocate(adj_imag, adj_real)
    end do

  end subroutine calc_adjoint_source

  subroutine initialize(this)
    class(ExponentiatedPhaseMisfit), intent(inout) :: this

    if (allocated(this%residuals)) deallocate(this%residuals)
    if (allocated(this%misfits)) deallocate(this%misfits)
    if (allocated(this%misfits_real)) deallocate(this%misfits_real)
    if (allocated(this%misfits_imag)) deallocate(this%misfits_imag)
    if (allocated(this%select_meas)) deallocate(this%select_meas)
    if (allocated(this%imeas)) deallocate(this%imeas)
    if (allocated(this%residuals_real)) deallocate(this%residuals_real)
    if (allocated(this%residuals_imag)) deallocate(this%residuals_imag)

    allocate(this%misfits_real(this%nwin))
    allocate(this%misfits_imag(this%nwin))
    allocate(this%select_meas(this%nwin))
    allocate(this%imeas(this%nwin))
    allocate(this%residuals_real(this%nwin))
    allocate(this%residuals_imag(this%nwin))
    allocate(this%misfits(this%nwin))
    allocate(this%residuals(this%nwin))

    this%select_meas = .true.
    this%imeas = IMEAS_EXP_PHASE
    this%misfits_real = 0.0_dp
    this%misfits_imag = 0.0_dp
    this%residuals_real = 0.0_dp
    this%residuals_imag = 0.0_dp
    this%total_misfit = 0.0_dp
    this%misfits = 0.0_dp
    this%residuals = 0.0_dp

  end subroutine initialize

end module exponentiated_phase_misfit