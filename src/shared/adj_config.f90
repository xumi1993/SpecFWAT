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
                      TSHIFT_MAX=5.0_dp, DLNA_MIN=-1.5_dp, DLNA_MAX=1.5_dp, CC_MIN=0.7_dp, &
                      target_dt=0.01_dp
    integer :: itaper_type=1, imeasure_type=1, min_cycle_in_window=3, &
                num_taper=5
    logical :: use_cc_error=.true., use_mt_error=.false., use_gpu=.false.
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

  
  type :: WINConfig
    real(kind=dp) :: jump_fac = 0.1, min_velocity = 2.4, sliding_win_len_fac = 3.0, &
                     threshold_shift_fac = 0.3, threshold_corr = 0.7, min_win_len_fac = 1.5, &
                     min_snr_window = 5.0, min_energy_ratio = 5.0, &
                     stalta_water_level = 0.08_dp, threshold_dlna = 0.8_dp, &
                     c_0 = 1.0_dp, c_2 = 0.0_dp, &
                     c_3a = 4.0_dp, c_3b = 2.5_dp, c_4a = 2.0_dp, c_4b = 6.0_dp, &
                     snr_integrate_base = 3.5_dp, &
                     max_time_before_first_arrival = 50.0_dp
    integer :: min_peaks_troughs = 3
    logical :: is_split_phases = .true.
  end type WINConfig

  type(ADJConfig), save :: adj_config_global
  type(WINConfig), save :: win_config_global
end module adj_config