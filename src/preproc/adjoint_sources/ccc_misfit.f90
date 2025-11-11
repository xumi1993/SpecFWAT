module ccc_misfit
  use config
  use signal
  use adj_config,  cfg => adj_config_global
  implicit none

  type, extends(AdjointMeasurement) :: CCCMisfit
    real(kind=dp), dimension(:), allocatable :: tshift, dlna, sigma_dt, sigma_dlna, &
                                            misfit_p, misfit_q, cc_max
  contains
    procedure :: calc_adjoint_source
    procedure, private :: initialize
  end type CCCMisfit

contains

  subroutine calc_adjoint_source(this, dat, syn, dt, windows)
    class(CCCmisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(in) :: windows
    real(kind=dp), dimension(:), allocatable :: s, d, adj_tw
    real(kind=dp) :: p, d_eng, s_eng, coef
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

      p = sum(d * s) * dt
      d_eng = sum(d * d) * dt
      s_eng = sum(s * s) * dt

      ! calculate cross-correlation coefficient
      coef = p / sqrt(d_eng * s_eng)
      if (coef < cfg%CC_MIN) then
        this%select_meas(iwin) = .false.
        this%misfits(iwin) = 0.0_dp
        this%residuals(iwin) = 0.0_dp
        this%imeas(iwin) = 0
        cycle
      end if

      ! calculate cross-correlation adjoint source
      adj_tw = (d - (p/s_eng) * s) / sqrt(d_eng * s_eng)
      call window_taper(adj_tw, cfg%taper_percentage, cfg%itaper_type)

      ! calculate cross-correlation coefficient misfit
      this%misfits(iwin) = 1 - coef
      this%total_misfit = this%total_misfit + this%misfits(iwin)
      this%residuals(iwin) = coef
      this%imeas(iwin) = IMEAS_CCC
      this%adj_src(nb:ne) = adj_tw
    end do

  end subroutine calc_adjoint_source

  subroutine initialize(this)
    class(CCCmisfit), intent(inout) :: this

    ! deallocate if already allocated
    if(allocated(this%misfits)) deallocate(this%misfits)
    if(allocated(this%select_meas)) deallocate(this%select_meas)
    if(allocated(this%residuals)) deallocate(this%residuals)
    if(allocated(this%imeas)) deallocate(this%imeas)

    ! allocate measurement arrays
    allocate(this%misfits(this%nwin))
    allocate(this%select_meas(this%nwin))
    allocate(this%residuals(this%nwin))
    allocate(this%imeas(this%nwin))

    ! initialize values
    this%misfits = 0.0_dp
    this%select_meas = .true.
    this%total_misfit = 0.0_dp
    this%imeas = 0

  end subroutine initialize

end module ccc_misfit