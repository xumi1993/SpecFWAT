module tele_data
  use config
  use ma_constants
  use signal, only: bandpass_dp, interpolate_syn_dp, detrend, demean, myconvolution_dp
  use syn_data, only: SynData
  use obs_data, only: ObsData
  use input_params, fpar => fwat_par_global
  use fk_coupling
  use fwat_mpi
  use sacio
  use logger, only: log
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, nrec, nrec_local, &
                        number_receiver_global, ispec_selected_rec, islice_selected_rec

  implicit none

  integer, private :: ier

  type, extends(SynData) :: TeleData
    real(kind=cr), dimension(:), pointer :: ttp ! travel time of P
    real(kind=cr) :: baz, az
    integer :: ttp_win
    contains
    procedure :: preprocess, calc_fktimes, semd2sac!, calc_stf
    procedure, private :: interpolate
    procedure :: finalize
  end type TeleData

contains

  subroutine semd2sac(this, ievt, is_conv_stf)
    class(TeleData), intent(inout) :: this
    integer, intent(in) :: ievt
    logical, intent(in), optional :: is_conv_stf
    integer :: irec, irec_local, icomp
    logical :: is_conv_stf_local
    real(kind=dp), dimension(:), allocatable :: stf_array, seismo_syn
    character(len=MAX_STRING_LEN) :: datafile
    type(sachead) :: header
    real(kind=cr) :: az, baz, dist, gcarc

    if (present(is_conv_stf)) then
      is_conv_stf_local = is_conv_stf
    else
      is_conv_stf_local = .false.
    endif

    call this%read(ievt)

    call this%od%read_stations(this%ievt)
    
    call this%calc_fktimes()

    if (worldrank == 0) call system('mkdir -p '//trim(fpar%acqui%in_dat_path(this%ievt)))
    call synchronize_all()

    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        if (is_conv_stf_local) then
          ! read stf
          call read_stf(stf_array)
        endif
        do icomp = 1, fpar%sim%NRCOMP
          if (is_conv_stf_local) then
            ! convolve stf with seismogram
            call myconvolution_dp(this%data(:, icomp, irec), stf_array, seismo_syn, 0)
            seismo_syn = seismo_syn * fpar%sim%dt
          else
            seismo_syn = this%data(:, icomp, irec)
          endif
          
          ! write SAC files
          datafile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%od%netwk(irec))//'.'&
                     //trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)&
                     //trim(fpar%sim%RCOMPS(icomp))//'.sac'
          call sacio_newhead(header, real(fpar%sim%dt), fpar%sim%nstep, -real(T0))
          ! header%dist = dist
          header%az = this%az
          header%baz = this%baz
          ! header%gcarc = gcarc
          ! header%evla = fpar%acqui%evla(this%ievt)
          ! header%evlo = fpar%acqui%evlo(this%ievt)
          ! header%evdp = fpar%acqui%evdp(this%ievt)
          header%stla = this%od%stla(irec)
          header%stlo = this%od%stlo(irec)
          header%stel = this%od%stel(irec)
          header%knetwk = this%od%netwk(irec)
          header%kstnm = this%od%stnm(irec)
          header%kcmpnm = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))
          header%t0 = this%ttp(irec)
          call sacio_writesac(datafile, header, seismo_syn, ier)
          if (ier /= 0) call exit_MPI(0, 'Error writing SAC file '//trim(datafile))
        enddo
      enddo
    endif
    call synchronize_all()

    contains
      subroutine read_stf(datarray)
        character(len=MAX_STRING_LEN) :: datafile
        type(sachead) :: header
        real(kind=dp), dimension(:), allocatable, intent(out) :: datarray
        real(kind=dp) :: thalf

        datafile = trim(SRC_REC_DIR)//'/STF_'//trim(fpar%acqui%evtid_names(this%ievt))//'.sac'
        call sacio_readsac(datafile, header, datarray, ier)
        if (ier /= 0) call exit_MPI(0, 'Error reading STF file '//trim(datafile))
        ! interpolate to the same time step
        if (abs(header%delta - fpar%sim%dt) > 1.0e-3 .or. header%npts /= fpar%sim%nstep) then
          thalf = (header%npts - 1) * header%delta / 2.0_dp
          call interpolate_syn_dp(datarray, -thalf, dble(header%delta), header%npts, &
                                  -dble(T0), dble(fpar%sim%dt), fpar%sim%nstep)
        endif
      end subroutine read_stf
  end subroutine semd2sac

  subroutine preprocess(this, ievt)
    class(TeleData), intent(inout) :: this
    integer, intent(in) :: ievt

    call this%read(ievt)
    call this%od%read_stations(this%ievt)
    call this%od%read_obs_data()

    ! calculate fk times
    call this%calc_fktimes()

    ! interpolate observed data to synthetic data
  
  end subroutine preprocess

  subroutine interpolate(this)
    class(TeleData), intent(inout) :: this
    real(kind=cr) :: tstart
    integer :: irec, icomp

    if (noderank == 0) then
      do irec = 1, this%nrec
        tstart = this%ttp(irec) - (this%od%tarr(irec) - this%od%tbeg(irec))
        do icomp = 1, fpar%sim%NRCOMP
          call interpolate_syn_dp(this%od%data(:, icomp, irec), dble(tstart), dble(this%od%dt),&
                                  this%od%npts, -dble(T0), dble(fpar%sim%dt), fpar%sim%nstep)
          call detrend(this%od%data(:, icomp, irec))
          call demean(this%od%data(:, icomp, irec))
          call bandpass_dp(this%od%data(:, icomp, irec), fpar%sim%nstep, dble(fpar%sim%dt),&
                           1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
        enddo
      enddo
    endif
    call synchronize_all()
    call this%filter(1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)

  end subroutine interpolate

  subroutine calc_fktimes(this)
    class(TeleData), intent(inout) :: this

    ! calculate fault rupture times
    real(kind=dp) :: stla, stlo
    real(kind=dp) :: xx, yy, zz
    integer :: irec

    call prepare_shm_array_cr_1d(this%ttp, this%nrec, this%ttp_win)

    call read_fk_model(fpar%acqui%evtid_names(this%ievt))

    this%baz = -phi_FK - 90.d0
    this%az = 90.d0 - phi_FK
    
    if (worldrank == 0) then
      do irec = 1, this%nrec
        stlo = dble(this%od%stlo(irec))
        stla = dble(this%od%stla(irec))
        if (.not. SUPPRESS_UTM_PROJECTION) then
          call utm_geo(stlo, stla, xx,yy, ILONGLAT2UTM)
        else
          xx = this%od%stlo(irec)
          yy = this%od%stla(irec)
        endif
        zz = this%od%stel(irec)
        call fktime(real(xx), real(yy), real(zz), this%ttp(irec))
      end do
      this%ttp = this%ttp - T0
    endif
    call sync_from_main_rank_cr_1d(this%ttp, this%nrec)
    call free_fk_arrays()
  end subroutine calc_fktimes

  subroutine finalize(this)
    class(TeleData), intent(inout) :: this

    call this%od%finalize()
    call free_shm_array(this%ttp_win)
    call free_shm_array(this%dat_win)
  end subroutine finalize

end module tele_data