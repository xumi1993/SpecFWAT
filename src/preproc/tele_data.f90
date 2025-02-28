module tele_data
  use config
  use ma_constants
  use singal, only: bandpass_dp, interpolate_syn_dp
  use syn_data, only: SynData
  use obs_data, only: ObsData
  use input_params, fpar => fwat_par_global
  use fk_coupling
  use fwat_mpi
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION, ILONGLAT2UTM
  use specfem_par, only: T0

  implicit none

  type, extends(SynData) :: TeleData
    type(ObsData) :: od
    real(kind=cr), dimension(:), pointer :: ttp ! travel time of P
    integer :: ttp_win
    contains
    procedure :: preprocess, calc_stf, calc_fktimes
    procedure, private :: interpolate
    procedure :: finalize
  end type TeleData

contains

  subroutine preprocess(this)
    class(TeleData), intent(inout) :: this

    call this%read(this%ievt)
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
        tstart = this%ttp(irec) - (this%od%arr(irec) - this%od%tbeg(irec))
        do icomp = 1, fpar%acqui%NRCOMP
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
    
    if (worldrank == 0) then
      do irec = 1, this%nrec
        stlo = dble(this%od%stlo(irec))
        stla = dble(this%od%stla(irec))
        if (.not. SUPPRESS_UTM_PROJECTION) then
          call utm_geo(stlo(irec),stla(irec),xx,yy,ILONGLAT2UTM)
        else
          xx = stlo(irec)
          yy = stla(irec)
        endif
        zz = this%od%stel(irec)
        call fktime(real(xx), real(yy), real(zz), this%ttp(irec))
      end do
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