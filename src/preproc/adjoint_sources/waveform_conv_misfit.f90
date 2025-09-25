module waveform_conv_misfit
  use config
  use signal
  use adj_config, cfg => adj_config_global
  use utils, only: simpson
  implicit none

  type, extends(AdjointMeasurement) :: WaveformConvMisfit
    real(kind=dp), dimension(:), allocatable :: adj_src_r, adj_src_z
    real(kind=dp) :: dt, tlen
    contains
      procedure :: calc_adjoint_source
      procedure, private :: initialize
  end type WaveformConvMisfit

contains
  subroutine calc_adjoint_source(this, dat_r, dat_z, syn_r, syn_z, dt, window)
    class(WaveformConvMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat_r, dat_z, syn_r, syn_z
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(2), intent(in) :: window
    real(kind=dp), dimension(:), allocatable :: c1, c2, &
                                                adj_tw_r, adj_tw_z, conv_diff
    integer :: nlen_win, nlen, nb, ne

    this%nwin = 1
    call this%initialize()

    nlen = size(dat_r)
    if (allocated(this%adj_src_r)) deallocate(this%adj_src_r)
    if (allocated(this%adj_src_z)) deallocate(this%adj_src_z)
    allocate(this%adj_src_r(nlen))
    allocate(this%adj_src_z(nlen))
    this%adj_src_r = 0.0_dp
    this%adj_src_z = 0.0_dp

    call get_window_info(window, dt, nb, ne, nlen_win)
   
    ! calculate adjoint sources
    call myconvolution_dp(syn_r, dat_z, c1, 0)
    call myconvolution_dp(dat_r, syn_z, c2, 0)
    conv_diff = (c1 - c2) * dt
    call myconvolution_dp(dat_z(nlen:1:-1), &
                          conv_diff, adj_tw_r, 0)
    adj_tw_r = adj_tw_r * dt
    adj_tw_r = adj_tw_r(nb:ne)
    call myconvolution_dp(-dat_r(nlen:1:-1), &
                          conv_diff, adj_tw_z, 0)
    adj_tw_z = adj_tw_z * dt
    adj_tw_z = adj_tw_z(nb:ne)

    ! add window to adjoint source
    conv_diff = conv_diff(nb:ne)
    call window_taper(conv_diff, cfg%taper_percentage, cfg%itaper_type)
    call window_taper(adj_tw_r, cfg%taper_percentage, cfg%itaper_type)
    call window_taper(adj_tw_z, cfg%taper_percentage, cfg%itaper_type)

    this%misfits(1) = 0.5_dp * simpson((conv_diff)**2, dt)
    this%total_misfit = this%total_misfit + this%misfits(1)
    this%residuals(1) = sum(conv_diff) / nlen_win
    this%adj_src_r(nb:ne) = adj_tw_r
    this%adj_src_z(nb:ne) = adj_tw_z
    
  end subroutine calc_adjoint_source

  subroutine initialize(this)
    class(WaveformConvMisfit), intent(inout) :: this

    ! Initialize the adjoint source and measurement arrays
    if (allocated(this%misfits)) deallocate(this%misfits)
    allocate(this%misfits(this%nwin))
    if (allocated(this%residuals)) deallocate(this%residuals)
    allocate(this%residuals(this%nwin))
    if (allocated(this%imeas)) deallocate(this%imeas)
    allocate(this%imeas(this%nwin))
    this%total_misfit = 0.0_dp
    this%misfits = 0.0_dp
    this%residuals = 0.0_dp
    this%imeas = IMEAS_WAVEFORM_CONV
  end subroutine initialize
end module waveform_conv_misfit