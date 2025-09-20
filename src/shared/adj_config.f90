module adj_config
  use config
  use signal
  implicit none

  type :: AdjointMeasurement
    real(kind=dp), dimension(:), allocatable :: adj_src
    real(kind=dp), dimension(:,:), allocatable :: window_chi
    real(kind=dp), dimension(:), allocatable :: misfits, residuals, errors
    real(kind=dp) :: total_misfit
    integer :: nwin
    integer, dimension(:), allocatable :: imeas
    logical, dimension(:), allocatable :: select_meas

  end type AdjointMeasurement

  type :: ADJConfig
    real(kind=dp) :: min_period, max_period, taper_percentage=0.3_dp, &
                      dt_sigma_min=1.0_dp, dlna_sigma_min=0.5_dp, transfunc_waterlevel=1e-10_dp, &
                      water_threshold=0.02_dp, wtr_env=0.2_dp, mt_nw=4.0_dp, phase_step=1.5_dp, &
                      dt_fac=2.0_dp, err_fac=2.5_dp, dt_max_scale=3.5_dp, TSHIFT_MIN=-5.0_dp, & 
                      TSHIFT_MAX=5.0_dp, DLNA_MIN=-1.5_dp, DLNA_MAX=1.5_dp, CC_MIN=0.8_dp
    integer :: itaper_type=1, imeasure_type=1, min_cycle_in_window=3, &
                num_taper=5
    logical :: use_cc_error=.true., use_mt_error=.false.
    ! itaper_type: 1=Hanning, 2=Hamming, 3=cos, 4=cos^10
    ! imeasure_type:
    !   1: waveform difference
    !   2: reveiver function
    !   3: Cross-convoluton
    !   11: cross-correlation traveltime (CC-TT)
    !   12: cross-correlation amplitude (CC-DLNA)
    !   13: multitaper traveltime (MT-TT)
    !   14: multitaper amplitude (MT-DLNA)
    ! use_cc_error: use cross-correlation error estimates for weighting
    ! use_mt_error: use multitaper error estimates for weighting
  end type ADJConfig

  type(ADJConfig), save :: adj_config_global
end module adj_config