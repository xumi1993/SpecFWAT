module adjoint_source
  use cc_tt_misfit, only: CCTTMisfit
  use mt_tt_misfit, only: MTTTMisfit  
  use exponentiated_phase_misfit, only: ExponentiatedPhaseMisfit
  use waveform_misfit, only: WaveformMisfit
  use waveform_conv_misfit, only: WaveformConvMisfit
  use rf_misfit, only: RFMisfit
  use adj_config, cfg => adj_config_global
  implicit none

contains

  subroutine calculate_adjoint_source(dat, syn, dt, windows, short_p, long_p, misfit_out)
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt, short_p, long_p
    real(kind=dp), dimension(:,:), intent(in) :: windows
    class(AdjointMeasurement), intent(out), allocatable :: misfit_out
    
    cfg%min_period = short_p
    cfg%max_period = long_p
    select case (cfg%imeasure_type)
      case (IMEAS_WAVEFORM) ! Waveform difference
        allocate(WaveformMisfit :: misfit_out)
        select type (misfit_out)
        type is (WaveformMisfit)
          call misfit_out%calc_adjoint_source(dat, syn, dt, windows)
        end select
        
      case (IMEAS_WAVEFORM_CONV) ! Waveform convolution
        ! Note: This requires radial and vertical components
        ! Use calc_adjoint_source_conv subroutine instead
        write(*,*) 'Warning: WAVEFORM_CONV misfit requires radial/vertical components'
        write(*,*) 'Use calc_adjoint_source_conv subroutine instead'
        error stop
        
      case (IMEAS_RF) ! Receiver function
        ! Note: This requires additional parameters for RF calculation
        write(*,*) 'Warning: RF misfit requires additional parameters'
        write(*,*) 'Use calc_adjoint_source_rf subroutine instead'
        error stop
        
      case (IMEAS_EXP_PHASE) ! Exponentiated phase
        allocate(ExponentiatedPhaseMisfit :: misfit_out)
        select type (misfit_out)
        type is (ExponentiatedPhaseMisfit)
          call misfit_out%calc_adjoint_source(dat, syn, dt, windows)
        end select
        
      case (IMEAS_CC_TT, IMEAS_CC_DLNA) ! Cross-correlation traveltime/amplitude
        allocate(CCTTMisfit :: misfit_out)
        select type (misfit_out)
        type is (CCTTMisfit)
          call misfit_out%calc_adjoint_source(dat, syn, dt, windows)
        end select
        
      case (IMEAS_CC_TT_MT, IMEAS_CC_DLNA_MT) ! Multitaper traveltime/amplitude
        allocate(MTTTMisfit :: misfit_out)
        select type (misfit_out)
        type is (MTTTMisfit)
          call misfit_out%calc_adjoint_source(dat, syn, dt, windows)
        end select
        
      case default
        write(*,'("Error: Unknown measurement type: ",I0)') cfg%imeasure_type
        write(*,*) 'Supported measurement types:'
        write(*,*) '  1 (IMEAS_WAVEFORM): Waveform difference'
        write(*,*) '  2 (IMEAS_WAVEFORM_CONV): Waveform convolution'
        write(*,*) '  3 (IMEAS_RF): Receiver function'
        write(*,*) '  4 (IMEAS_EXP_PHASE): Exponentiated phase'
        write(*,*) ' 11 (IMEAS_CC_TT): Cross-correlation traveltime'
        write(*,*) ' 12 (IMEAS_CC_DLNA): Cross-correlation amplitude'
        write(*,*) ' 13 (IMEAS_CC_TT_MT): Multitaper traveltime'
        write(*,*) ' 14 (IMEAS_CC_DLNA_MT): Multitaper amplitude'
        error stop
    end select

    call filter_adj(misfit_out%adj_src, dt, windows)


  end subroutine calculate_adjoint_source


  subroutine filter_adj(adj_src, dt, windows)
    use signal, only: bandpass_dp
    real(kind=dp), dimension(:), intent(inout) :: adj_src
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(in) :: windows
    real(kind=dp), dimension(2) :: window
    real(kind=dp), dimension(:), allocatable :: adj_tw
    integer :: nlen_win, nb, ne

    call bandpass_dp(adj_src, size(adj_src), dt, real(1/cfg%max_period), real(1/cfg%min_period), IORD)
    window(1) = minval(windows)
    window(2) = maxval(windows)

    call get_window_info(window, dt, nb, ne, nlen_win)
    allocate(adj_tw(nlen_win))
    adj_tw(:) = adj_src(nb:ne)

    ! taper the windows
    adj_src = 0.0_dp
    call window_taper(adj_tw, cfg%taper_percentage, cfg%itaper_type)
    adj_src(nb:ne) = adj_tw

  end subroutine filter_adj

end module adjoint_source