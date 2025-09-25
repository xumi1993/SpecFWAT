module waveform_misfit
  use adj_config, cfg => adj_config_global
  use config
  use signal
  use utils, only: simpson
  implicit none

  type, extends(AdjointMeasurement) :: WaveformMisfit
  contains
    procedure :: calc_adjoint_source
    procedure, private :: initialize
  end type WaveformMisfit

contains
  subroutine calc_adjoint_source(this, dat, syn, dt, windows)
    class(WaveformMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(in) :: windows
    real(kind=dp), dimension(:), allocatable :: s, d, adj_tw
    real(kind=dp) :: misfit
    integer :: iwin, nlen, nb, ne, nlen_win

    ! num of measurements
    if (size(windows, 2) /= 2) then
      write(*,*) 'Error: windows must have two columns (start and end times)'
      error stop
    end if
  
    ! allocate measurement arrays
    this%nwin = size(windows, 1)
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

      adj_tw = s - d

      call window_taper(adj_tw, cfg%taper_percentage, cfg%itaper_type)

      ! calculate waveform misfit and adjoint source
      misfit = 0.5_dp * simpson((s - d)**2, dt)
      this%misfits(iwin) = misfit
      this%residuals(iwin) = sum(s - d) / nlen_win
      this%total_misfit = this%total_misfit + misfit
      this%adj_src(nb:ne) = adj_tw

      deallocate(s, d, adj_tw)
    end do

  end subroutine calc_adjoint_source

  subroutine initialize(this)
    class(WaveformMisfit), intent(inout) :: this

    ! Initialize the adjoint source and measurement arrays

    if (allocated(this%misfits)) deallocate(this%misfits)
    allocate(this%misfits(this%nwin))
    if (allocated(this%imeas)) deallocate(this%imeas)
    allocate(this%imeas(this%nwin))
    if (allocated(this%select_meas)) deallocate(this%select_meas)
    allocate(this%select_meas(this%nwin))
    if (allocated(this%residuals)) deallocate(this%residuals)
    allocate(this%residuals(this%nwin))
    this%misfits = 0.0_dp
    this%imeas = IMEAS_WAVEFORM
    this%select_meas = .true.
    this%total_misfit = 0.0_dp

  end subroutine initialize

end module waveform_misfit