module noise_data
  use config
  use misfit_mod
  use adjoint_source, only: calculate_adjoint_source
  use adj_config, only: AdjointMeasurement
  use syn_data, only: SynData
  use input_params, fpar => fwat_par_global
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, DT, NSTEP,OUTPUT_FILES
  use fwat_mpi
  use utils, only: zeros_dp, zeros
  use distaz_lib
  use common_lib, only: get_band_name, get_icomp_syn, rotate_ZRT_to_ZNE, mkdir
  use signal, only: interpolate_func_dp, dif1, detrend, demean, bandpass_dp
  use sacio
  use logger, only: log
  implicit none

  character(len=MAX_STRING_LEN), private :: msg

  integer, private :: ier
  type, extends(SynData) :: NoiseData
  contains
    procedure :: semd2sac, preprocess, finalize, write_in_preocess, write_adj_src_sac
    procedure, private :: calc_distaz, measure_adj
  end type NoiseData

contains
  subroutine semd2sac(this, ievt)
    class(NoiseData), intent(inout) :: this
    integer, intent(in) :: ievt
    integer :: irec_local, irec, icomp, icomp_syn
    character(len=MAX_STRING_LEN) :: datafile
    type(sachead) :: header

    call this%init(ievt)

    call mkdir(fpar%acqui%in_dat_path(this%ievt))

    call this%od%read_stations(ievt, .true.)

    call this%calc_distaz()

    call this%read(this%od%baz)

    if (this%nrec_loc > 0) then
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
        do icomp = 1, fpar%sim%NRCOMP
          icomp_syn = get_icomp_syn(fpar%sim%RCOMPS(icomp))
          datafile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%od%netwk(irec))//'.'&
                      //trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)&
                      //trim(fpar%sim%RCOMPS(icomp))//'.sac'
          call sacio_newhead(header, real(fpar%sim%dt), fpar%sim%nstep, -real(T0))
          header%baz = this%od%baz(irec)
          header%dist = this%od%dist(irec)
          header%evla = fpar%acqui%evla(this%ievt)
          header%evlo = fpar%acqui%evlo(this%ievt)
          header%evdp = fpar%acqui%evdp(this%ievt)
          header%stla = this%od%stla(irec)
          header%stlo = this%od%stlo(irec)
          header%stel = this%od%stel(irec)
          header%knetwk = trim(this%od%netwk(irec))
          header%kstnm = trim(this%od%stnm(irec))
          header%kcmpnm = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))
          call sacio_writesac(datafile, header, this%data(:, icomp_syn, irec), ier)
          if (ier /= 0) call exit_MPI(0, 'Error writing SAC file '//trim(datafile))
        end do 
      end do
    end if
    call synchronize_all()

  end subroutine semd2sac

  subroutine preprocess(this, ievt)
    class(NoiseData), intent(inout) :: this
    integer, intent(in) :: ievt

    call this%init(ievt)

    call this%od%read_stations(ievt, .true.)

    call this%od%read_obs_data()

    call this%read(this%od%baz)

    call this%measure_adj()

  end subroutine preprocess

  subroutine measure_adj(this)
    class(NoiseData), intent(inout) :: this
    integer :: irec_local, irec, iflt, icomp, icomp_syn, nb, ne
    real(kind=dp), dimension(:,:), allocatable :: adj_loc
    real(kind=dp), dimension(:,:,:,:), allocatable :: adj_src
    real(kind=dp), dimension(:), allocatable :: seismo_dat,seismo_syn
    real(kind=dp), dimension(:), allocatable :: adj_2, adj_3, adj_1
    real(kind=dp), dimension(1, 2) :: windows 
    real(kind=dp) :: dist_min, max_amp, tstart, tend
    logical :: is_reject
    type(RECMisfit), dimension(:), allocatable :: recm
    type(EVTMisfit) :: evtm
    class(AdjointMeasurement), allocatable :: misfit_out

    ! allocate for adjoint src
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
          ! initialize recm
          recm(irec_local) = RECMisfit(fpar%sim%NRCOMP)
          ! Initialize TraceMisfit for each component with 1 window
          irec = select_global_id_for_rec(irec_local)
          ! time window
          tstart = this%od%dist(irec)/fpar%sim%GROUPVEL_MAX(iflt)-fpar%sim%LONG_P(iflt)/2.
          tend = this%od%dist(irec)/fpar%sim%GROUPVEL_MIN(iflt)+fpar%sim%LONG_P(iflt)/2.
          tstart = max(tstart, -t0)
          tend = min(tend, (NSTEP-2)*dble(DT)-t0)
          windows(1, :) = [tstart, tend] + t0
          recm(irec_local)%sta = trim(this%od%stnm(irec))
          recm(irec_local)%net = trim(this%od%netwk(irec))
          do icomp = 1, fpar%sim%NRCOMP
            ! get data component
            seismo_dat = 0.0_dp
            seismo_syn = 0.0_dp
            seismo_dat = interpolate_func_dp(this%od%data(:, icomp, irec), dble(this%od%tbeg(irec)), &
                                             dble(this%od%dt), this%od%npts, -dble(t0), dble(DT), NSTEP)
            if (.not. fpar%sim%SUPPRESS_EGF) call dif1(seismo_dat, dble(DT))
            call detrend(seismo_dat)
            call demean(seismo_dat)
            call bandpass_dp(seismo_dat, fpar%sim%nstep, dble(fpar%sim%dt),&
                             1/fpar%sim%LONG_P(iflt), 1/fpar%sim%SHORT_P(iflt), IORD)
            
            ! get synthetic component
            icomp_syn = get_icomp_syn(fpar%sim%RCOMPS(icomp))
            seismo_syn = this%data(:, icomp_syn, irec)
            call bandpass_dp(seismo_syn, fpar%sim%nstep, dble(fpar%sim%dt),&
                             1/fpar%sim%LONG_P(iflt), 1/fpar%sim%SHORT_P(iflt), IORD)

            ! reject near offsets
            is_reject = .false.
            if (.not. fpar%sim%USE_NEAR_OFFSET) then
              dist_min = (fpar%sim%GROUPVEL_MAX(iflt) + fpar%sim%GROUPVEL_MIN(iflt))*fpar%sim%LONG_P(iflt)/2.
              if (this%od%dist(irec) < dist_min) then
                seismo_syn = 1.0e-18_dp
                seismo_dat = 1.0e-18_dp
                is_reject = .true.
              end if
            end if

            ! normalization
            nb = floor((tstart+ t0)/dble(DT)+1)
            ne = floor((tend+ t0)/dble(DT)+1)
            seismo_dat = seismo_dat / maxval(abs(seismo_dat(nb:ne))) * maxval(abs(seismo_syn(nb:ne)))
            if (is_output_preproc) then
              call this%write_in_preocess(irec, icomp, tstart, tend, 'syn', seismo_syn)
              call this%write_in_preocess(irec, icomp, tstart, tend, 'obs', seismo_dat)
            end if

            ! calculate adjoint source
            recm(irec_local)%trm(icomp) = TraceMisfit(1)
            if (.not. is_reject) then
              call calculate_adjoint_source(seismo_dat, seismo_syn, dble(fpar%sim%dt), windows, &
                                            dble(fpar%sim%SHORT_P(iflt)), dble(fpar%sim%LONG_P(iflt)), misfit_out)
              adj_src(:, icomp, irec_local, iflt) = misfit_out%adj_src(:) * fpar%acqui%src_weight(this%ievt)
              recm(irec_local)%trm(icomp)%misfits(1) = misfit_out%misfits(1) * fpar%acqui%src_weight(this%ievt)
              recm(irec_local)%trm(icomp)%residuals(1) = misfit_out%residuals(1)
              recm(irec_local)%trm(icomp)%imeas(1) = misfit_out%imeas(1)
              recm(irec_local)%trm(icomp)%tstart(1) = tstart
              recm(irec_local)%trm(icomp)%tend(1) = tend
              recm(irec_local)%total_misfit = recm(irec_local)%total_misfit + misfit_out%total_misfit &
                                                  * fpar%acqui%src_weight(this%ievt)
            else
              recm(irec_local)%trm(icomp)%imeas(1) = 0
            endif
            recm(irec_local)%chan(icomp) = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))

            if (IS_OUTPUT_ADJ_SRC) call this%write_adj_src_sac(adj_src(:, icomp, irec_local, iflt), irec, icomp)
          end do ! icomp
        end do ! irec_local
      end if ! this%nrec_loc > 0
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
    end do ! iflt
    call synchronize_all()
    if (allocated(seismo_dat)) deallocate(seismo_dat)
    if (allocated(seismo_syn)) deallocate(seismo_syn)

    call this%get_comp_name_adj()
    if (this%nrec_loc > 0) then
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
        do icomp = 1, fpar%sim%NRCOMP
          adj_loc = zeros_dp(NSTEP, 3)
          if (fpar%sim%ADJ_SRC_NORM) then
            max_amp = maxval(abs(adj_src(:, icomp, irec_local, :)))
            do iflt = 1, fpar%sim%NUM_FILTER
              if (maxval(abs(adj_src(:, icomp, irec_local, iflt))) == 0.) cycle
              adj_loc(:, icomp) = adj_loc(:, icomp) + adj_src(:, icomp, irec_local, iflt)&
                                  / maxval(abs(adj_src(:, icomp, irec_local, iflt))) * max_amp
            end do
          else
            adj_loc(:, icomp) = sum(adj_src(:, icomp, irec_local, :), dim=2)
          end if
        end do
        adj_1 = zeros_dp(NSTEP)
        adj_2 = zeros_dp(NSTEP)
        adj_3 = zeros_dp(NSTEP)
        call rotate_ZRT_to_ZNE(adj_loc(:, 1), adj_loc(:, 2), adj_loc(:, 3), &
                               adj_1, adj_2, adj_3, NSTEP, dble(this%od%baz(irec)))
        call this%write_adj(adj_1, trim(this%comp_name(1)), irec)
        call this%write_adj(adj_2, trim(this%comp_name(2)), irec)
        call this%write_adj(adj_3, trim(this%comp_name(3)), irec)
      end do
    end if
    if(allocated(adj_loc)) deallocate(adj_loc)
    if(allocated(adj_src)) deallocate(adj_src)

    call synchronize_all()
  end subroutine measure_adj

  subroutine finalize(this)
    class(NoiseData), intent(inout) :: this

    call this%od%finalize()
    call free_shm_array(this%dat_win)
  end subroutine finalize

  subroutine write_in_preocess(this, irec, icomp, tb, te, label, trace)
    class(NoiseData), intent(inout) :: this
    integer, intent(in) :: irec, icomp
    real(kind=dp), intent(in) :: tb, te
    real(kind=dp), dimension(:), intent(in) :: trace
    character(len=*), intent(in) :: label
    type(sachead) :: header
    character(len=MAX_STRING_LEN) :: datafile

    datafile = trim(OUTPUT_FILES)//'/'//trim(this%od%netwk(irec))//'.'//&
               trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
               trim(fpar%sim%RCOMPS(icomp))//'.'//trim(label)//'.sac.'//trim(this%band_name)

    call sacio_newhead(header, real(DT), NSTEP, -real(T0))
    header%baz = this%od%baz(irec)
    header%dist = this%od%dist(irec)
    header%evla = fpar%acqui%evla(this%ievt)
    header%evlo = fpar%acqui%evlo(this%ievt)
    header%evdp = fpar%acqui%evdp(this%ievt)
    header%stla = this%od%stla(irec)
    header%stlo = this%od%stlo(irec)
    header%stel = this%od%stel(irec)
    header%knetwk = trim(this%od%netwk(irec))
    header%kstnm = trim(this%od%stnm(irec))
    header%kcmpnm = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))
    header%kevnm = trim(fpar%acqui%evtid_names(this%ievt))
    header%t1 = real(tb)
    header%t2 = real(te)
    call sacio_writesac(datafile, header, trace, ier)

  end subroutine write_in_preocess

  subroutine calc_distaz(this)
    class(NoiseData), intent(inout) :: this
    integer :: irec
    real(kind=dp) :: dist, bazi, azi, delta

    if (worldrank == 0) then
      do irec = 1, this%od%nsta
        if (SUPPRESS_UTM_PROJECTION) then
          call distaz_cart(dble(fpar%acqui%evlo(this%ievt)), dble(fpar%acqui%evla(this%ievt)), &
                           dble(this%od%stlo(irec)), dble(this%od%stla(irec)), azi, bazi, dist)
          dist = dist / 1000.0_dp
        else
          call distaz(dble(this%od%stla(irec)), dble(this%od%stlo(irec)), &
                      dble(fpar%acqui%evla(this%ievt)), dble(fpar%acqui%evlo(this%ievt)),&
                      azi, bazi, delta, dist)
        endif
        this%od%baz(irec) = real(bazi)
        this%od%dist(irec) = real(dist)
      end do
    end if
    call synchronize_all()
    call sync_from_main_rank_cr_1d(this%od%baz, this%od%nsta)
    call sync_from_main_rank_cr_1d(this%od%dist, this%od%nsta)

  end subroutine calc_distaz

  subroutine write_adj_src_sac(this, adj_src, irec, icomp)
    class(NoiseData), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: adj_src ! nstep, ncomp, nrec_loc
    integer, intent(in) :: irec, icomp
    type(sachead) :: header
    character(len=MAX_STRING_LEN) :: sacfile

    call sacio_newhead(header, real(DT), NSTEP, -real(T0))
    header%knetwk = trim(this%od%netwk(irec))
    header%kstnm = trim(this%od%stnm(irec))
    header%kcmpnm = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))
    sacfile = trim(fpar%acqui%out_fwd_path(this%ievt))//'/'//trim(ADJOINT_PATH)//&
          '/'//trim(this%od%netwk(irec))//'.'//trim(this%od%stnm(irec))//&
          '.'//trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))//&
          '.adj.sac.'//trim(this%band_name)
    call sacio_writesac(sacfile, header, adj_src, ier)
  end subroutine write_adj_src_sac

end module noise_data