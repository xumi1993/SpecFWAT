module rf_misfit
  use config
  use signal
  use adj_config, cfg => adj_config_global
  use decon_mod, only: deconit
  use utils, only: simpson

  type, extends(AdjointMeasurement) :: RFMisfit
    real(kind=dp), dimension(:), allocatable :: adj_src_r, adj_src_z
  contains
    procedure :: calc_adjoint_source
    procedure, private :: initialize
  end type RFMisfit

contains

  subroutine calc_adjoint_source(this, dat, syn, synr, synz, dt, tp, window, &
                                 f0, shift, maxit, minderr)
    class(RFMisfit), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat, syn, synr, synz
    real(kind=dp), intent(in) :: dt, f0, tp
    real(kind=dp), dimension(2), intent(in) :: window
    real(kind=dp), optional, intent(in) :: shift, minderr
    integer, optional, intent(in) :: maxit
    real(kind=dp), dimension(2) :: win_tr
    real(kind=dp) :: shift_loc, minderr_loc, shift_rev
    integer :: maxit_loc, nlen_win, nlen_win_rf, nb, ne, nlen, nlen_rf, nstart
    real(kind=dp), dimension(:), allocatable :: zrf, dat_norm, syn_norm, r_rev, z_rev, diff, &
                                                adj_r, adj_z, adj_r_tw, adj_z_tw, num, den, &
                                                dat_tw, syn_tw

    if (present(shift)) then
      shift_loc = shift
    else
      shift_loc = 5.0_dp
    end if

    if (present(maxit)) then
      maxit_loc = maxit
    else
      maxit_loc = 200
    end if

    if (present(minderr)) then
      minderr_loc = minderr
    else
      minderr_loc = 1.0e-3_dp
    end if

    this%nwin = 1
    nlen = size(synr)
    nlen_rf = size(dat)

    allocate(dat_norm(nlen), syn_norm(nlen))
    dat_norm = 0.0_dp
    syn_norm = 0.0_dp

    if (allocated(this%adj_src_r)) deallocate(this%adj_src_r)
    if (allocated(this%adj_src_z)) deallocate(this%adj_src_z)
    allocate(this%adj_src_r(nlen))
    allocate(this%adj_src_z(nlen))
    this%adj_src_r = 0.0_dp
    this%adj_src_z = 0.0_dp
    call this%initialize()

    ! calculate z-rf for normalization
    call deconit(synz, synz, real(dt), 10., real(f0), 10, 0.001, 0, zrf)

    if (nlen >= nlen_rf) then
      dat_norm(1:nlen_rf) = dat(1:nlen_rf) / maxval(abs(zrf))
      syn_norm(1:nlen_rf) = syn(1:nlen_rf) / maxval(abs(zrf))
    else
      dat_norm(1:nlen) = dat(1:nlen) / maxval(abs(zrf))
      syn_norm(1:nlen) = syn(1:nlen) / maxval(abs(zrf))
    end if

    ! trim for measurement
    nstart = int((tp - shift_loc)/ dt) + 1
    if (nstart < 1) error stop 'Error: Not enough samples for the beginning of the window'
    allocate(z_rev(nlen), r_rev(nlen))
    z_rev = 0.0_dp
    r_rev = 0.0_dp
    z_rev(1:nlen) = synz(nlen:nstart:-1)
    r_rev(1:nlen) = synr(nlen:nstart:-1)
    call window_taper(z_rev, cfg%taper_percentage, cfg%itaper_type)
    call window_taper(r_rev, cfg%taper_percentage, cfg%itaper_type)

    ! calculate adjoint source
    shift_rev = nlen * dt - shift_loc
    diff = syn_norm - dat_norm
    call deconit(diff, z_rev, real(dt), real(shift_rev), real(f0), maxit_loc, real(minderr_loc), 1, adj_r)
    call myconvolution_dp(-diff, r_rev, num, 1)
    call myconvolution_dp(z_rev, z_rev, den, 1)
    call deconit(num, den, real(dt), real(shift_rev), real(f0), maxit_loc, real(minderr_loc), 1, adj_z)
    deallocate(num, den, z_rev, r_rev)

    ! tapper
    nlen_win_rf = int((window(2) - window(1))/dt) + 1
    nb = nstart
    ne = nb + nlen_win_rf - 1
    adj_r_tw = adj_r(nb:ne)
    adj_z_tw = adj_z(nb:ne)
    call window_taper(adj_r_tw, cfg%taper_percentage, cfg%itaper_type)
    call window_taper(adj_z_tw, cfg%taper_percentage, cfg%itaper_type)
    this%adj_src_r(nb:ne) = adj_r_tw
    this%adj_src_z(nb:ne) = adj_z_tw

    ! trim data and synthetic
    win_tr = window + shift_loc
    if (win_tr(1) < 0.0_dp) win_tr(1) = 0.0_dp
    call get_window_info(win_tr, dt, nb, ne, nlen_win)
    dat_tw = dat_norm(nb:ne)
    syn_tw = syn_norm(nb:ne)
    call window_taper(dat_tw, cfg%taper_percentage, cfg%itaper_type)
    call window_taper(syn_tw, cfg%taper_percentage, cfg%itaper_type)

    ! calculate misfit
    this%misfits(1) = 0.5_dp * simpson((syn_tw - dat_tw)**2, dt)
    this%residuals(1) = sum(syn_tw - dat_tw) / nlen_win_rf
    this%total_misfit = this%total_misfit + this%misfits(1)

    deallocate(adj_r, adj_z, dat_norm, syn_norm, adj_r_tw, adj_z_tw)

  end subroutine calc_adjoint_source

  subroutine initialize(this)
    class(RFMisfit), intent(inout) :: this

    ! Initialize the adjoint source and measurement arrays
    if (allocated(this%misfits)) deallocate(this%misfits)
    allocate(this%misfits(this%nwin))
    if (allocated(this%residuals)) deallocate(this%residuals)
    allocate(this%residuals(this%nwin))
    if (allocated(this%imeas)) deallocate(this%imeas)
    allocate(this%imeas(this%nwin))
    if (allocated(this%select_meas)) deallocate(this%select_meas)
    allocate(this%select_meas(this%nwin))
    if (allocated(this%errors)) deallocate(this%errors)
    allocate(this%errors(this%nwin))
    this%total_misfit = 0.0_dp
    this%imeas = IMEAS_RF
    this%misfits = 0.0_dp
    this%residuals = 0.0_dp
    this%errors = 1.0_dp
    this%select_meas = .true.
  end subroutine initialize
end module rf_misfit