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

  subroutine calculate_adjoint_source(dat, syn, dt, windows, misfit_out)
    real(kind=dp), dimension(:), intent(in) :: dat, syn
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(in) :: windows
    class(AdjointMeasurement), intent(out), allocatable :: misfit_out
    
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

  end subroutine calculate_adjoint_source

end module adjoint_source