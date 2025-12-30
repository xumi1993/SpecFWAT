module leq_data
  use config
  use misfit_mod
  use adjoint_source, only: calculate_adjoint_source
  use adj_config, only: AdjointMeasurement
  use noise_data, only: NoiseData
  use input_params, fpar => fwat_par_global
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, DT, NSTEP,OUTPUT_FILES
  use fwat_mpi
  use utils, only: zeros_dp, zeros
  use distaz_lib
  use common_lib, only: get_band_name, get_icomp_syn, rotate_ZRT_to_ZNE, mkdir
  use signal, only: interpolate_func_dp, detrend, demean, bandpass_dp
  use win_sel, only: win_sel_type
  use flexwin, only: flexwin_type
  use travel_times_mod
  use sacio
  use logger, only: log
  implicit none

  character(len=MAX_STRING_LEN), private :: msg
 
  type, extends(NoiseData) :: LEQData
  contains
    procedure :: preprocess
    procedure, private :: measure_adj
  end type LEQData

contains

  subroutine preprocess(this, ievt)
    class(LEQData), intent(inout) :: this
    integer, intent(in) :: ievt

    call this%init(ievt)

    call this%od%read_stations(ievt, .true.)

    call this%od%read_obs_data()

    call this%read(this%od%baz)

    call this%measure_adj()

  end subroutine preprocess

  subroutine measure_adj(this)
    class(LEQData), intent(inout) :: this
    integer :: irec_local, irec, iflt, icomp, icomp_syn, iwin, nwin
    real(kind=dp), dimension(:,:,:,:), allocatable :: adj_src
    real(kind=dp), dimension(:), allocatable :: seismo_dat,seismo_syn
    real(kind=dp), dimension(:,:), allocatable :: windows
    real(kind=dp) :: tp
    type(RECMisfit), dimension(:), allocatable :: recm
    type(EVTMisfit) :: evtm
    class(AdjointMeasurement), allocatable :: misfit_out
    type(win_sel_type) :: win
    type(flexwin_type) :: fw

    if (this%nrec_loc > 0) then
      adj_src = zeros_dp(NSTEP, fpar%sim%NRCOMP, this%nrec_loc, fpar%sim%NUM_FILTER)
      seismo_dat = zeros_dp(NSTEP)
      seismo_syn = zeros_dp(NSTEP)
    endif
    call synchronize_all()
    do iflt = 1, fpar%sim%NUM_FILTER
      call get_band_name(fpar%sim%SHORT_P(iflt), fpar%sim%LONG_P(iflt), this%band_name)
      if (this%nrec_loc > 0) then
        ! allocate receiver misfit
        allocate(recm(this%nrec_loc))
        do irec_local = 1, this%nrec_loc
          irec = select_global_id_for_rec(irec_local)
          ! initialize recm
          recm(irec_local) = RECMisfit(fpar%sim%NRCOMP)
          recm(irec_local)%sta = trim(this%od%stnm(irec))
          recm(irec_local)%net = trim(this%od%netwk(irec))
          ! calculate theoretical P arrival time for LEQ
          if (fpar%sim%WIN_TYPE /= WIN_GROUPVEL_TYPE) then
            block
              integer :: nphases
              character(len=8), dimension(:), allocatable :: names
              real(kind=dp), dimension(:), allocatable :: times
              call ttimes(dble(this%od%gcarc(irec)),dble(fpar%acqui%evdp(irec)),fpar%sim%PHASE,nphases,names,times)
              tp = times(1)
              deallocate(names, times)
            end block
          endif
          do icomp = 1, fpar%sim%NRCOMP
            seismo_dat = 0.0_dp
            seismo_syn = 0.0_dp
            ! read observed data
            seismo_dat = interpolate_func_dp(this%od%data(:, icomp, irec), dble(this%od%tbeg(irec)), &
                                             dble(this%od%dt), this%od%npts, -dble(t0), dble(DT), NSTEP)
            call detrend(seismo_dat)
            call demean(seismo_dat)
            call window_taper(seismo_dat, taper_per, 1)
            call bandpass_dp(seismo_dat, fpar%sim%nstep, dble(fpar%sim%dt),&
                             1/fpar%sim%LONG_P(iflt), 1/fpar%sim%SHORT_P(iflt), IORD)
            ! read synthetic data
            icomp_syn = get_icomp_syn(fpar%sim%RCOMPS(icomp))
            seismo_syn = this%data(:, icomp_syn, irec)
            call window_taper(seismo_syn, taper_per, 1)
            call bandpass_dp(seismo_syn, fpar%sim%nstep, dble(fpar%sim%dt),&
                             1/fpar%sim%LONG_P(iflt), 1/fpar%sim%SHORT_P(iflt), IORD)

            ! Select windows
            select case (fpar%sim%WIN_TYPE)
              case (WIN_SELECTOR_TYPE)
                win = win_sel_type(seismo_dat, seismo_syn, dble(fpar%sim%dt), dble(T0), &
                              tp, dble(this%od%dist(irec)), dble(fpar%sim%SHORT_P(iflt)))
                call win%select_windows()
                nwin = win%n_win
                if (nwin > 0) windows = win%twin(:, :)
              case (WIN_FLEX_TYPE)
                ! dat, syn, dt, t0, tp, dis, min_period, max_period
                fw = flexwin_type(seismo_dat, seismo_syn, dble(fpar%sim%dt), dble(T0), &
                                  tp, dble(this%od%dist(irec)), dble(fpar%sim%SHORT_P(iflt)), dble(fpar%sim%LONG_P(iflt)))
                call fw%select_windows()
                nwin = fw%n_win
                if (nwin > 0) windows = fw%twin(:, :)
              case (WIN_GROUPVEL_TYPE)
                allocate(windows(1,2))
                windows(1,1) = this%od%dist(irec)/fpar%sim%GROUPVEL_MAX(iflt)-fpar%sim%LONG_P(iflt)/2.
                windows(1,2) = this%od%dist(irec)/fpar%sim%GROUPVEL_MIN(iflt)+fpar%sim%LONG_P(iflt)/2.
                windows(1,1) = max(windows(1,1), -t0)
                windows(1,2) = min(windows(1,2), (NSTEP-2)*dble(DT)-t0)
                nwin = 1
              case (WIN_ARRIVAL_TYPE)
                allocate(windows(1,2))
                windows(1,1) = tp + fpar%sim%TIME_WIN(1)
                windows(1,2) = tp + fpar%sim%TIME_WIN(2)
                windows(1,1) = max(windows(1,1), -t0)
                windows(1,2) = min(windows(1,2), (NSTEP-2)*dble(DT)-t0)
                nwin = 1
            end select

            ! calculate adjoint source for this component
            ! Initialize TraceMisfit for each component with n windows
            recm(irec_local)%trm(icomp) = TraceMisfit(nwin)
            if (nwin > 0) then
              if (IS_OUTPUT_PREPROC) then
                ! Write preprocessed data
                call this%write_in_preocess(irec, icomp, windows(1, 1), windows(nwin, 2), 'syn', seismo_syn)
                call this%write_in_preocess(irec, icomp, windows(1, 1), windows(nwin, 2), 'obs', seismo_dat)
              end if
              ! calculate adjoint source
              call calculate_adjoint_source(seismo_dat, seismo_syn, dble(fpar%sim%dt), windows+T0, &
                                            dble(fpar%sim%SHORT_P(iflt)), dble(fpar%sim%LONG_P(iflt)), misfit_out)
              adj_src(:, icomp, irec_local, iflt) = misfit_out%adj_src(:) * fpar%acqui%src_weight(this%ievt)
              ! store misfit results
              do iwin = 1, nwin
                recm(irec_local)%trm(icomp)%misfits(iwin) = misfit_out%misfits(iwin) * fpar%acqui%src_weight(this%ievt)
                recm(irec_local)%trm(icomp)%residuals(iwin) = misfit_out%residuals(iwin)
                recm(irec_local)%trm(icomp)%imeas(iwin) = misfit_out%imeas(iwin)
                recm(irec_local)%trm(icomp)%tstart(iwin) = windows(iwin, 1)
                recm(irec_local)%trm(icomp)%tend(iwin) = windows(iwin, 2)
                recm(irec_local)%total_misfit = recm(irec_local)%total_misfit + misfit_out%total_misfit &
                                                * fpar%acqui%src_weight(this%ievt)
              end do
            end if ! if (nwin > 0)
            recm(irec_local)%chan(icomp) = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))
            if (IS_OUTPUT_ADJ_SRC) call this%write_adj_src_sac(adj_src(:, icomp, irec_local, iflt), irec, icomp)
          enddo ! icomp
        enddo ! irec_local
      end if ! if (this%nrec_loc > 0)
      call synchronize_all()
      ! assemble and write event misfit for this frequency band
      evtm = EVTMisfit(fpar%acqui%evtid_names(this%ievt), this%nrec)
      evtm%band_name = trim(this%band_name)
      call evtm%assemble(recm)
      write(msg, '(a,f12.6)') 'Total misfit: of '//trim(this%band_name)//': ', evtm%total_misfit
      call log%write(msg, .true.)
      call evtm%write()
      if (allocated(recm)) deallocate(recm)
      this%total_misfit(iflt) = evtm%total_misfit
      call synchronize_all()
    enddo ! iflt

    call synchronize_all()
    if (allocated(seismo_dat)) deallocate(seismo_dat)
    if (allocated(seismo_syn)) deallocate(seismo_syn)

    call this%rotate_adj_src(adj_src)
    if(allocated(adj_src)) deallocate(adj_src)

    call synchronize_all()
  end subroutine measure_adj

end module leq_data