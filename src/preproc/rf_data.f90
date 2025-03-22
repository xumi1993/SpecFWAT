module rf_data
  use config
  use ma_constants
  use common_lib, only: get_band_name, rotate_R_to_NE_dp, dwascii
  use signal, only: bandpass_dp, interpolate_syn_dp, detrend, demean, &
                    myconvolution_dp, time_deconv
  use syn_data, only: SynData, average_amp_scale
  use obs_data, only: ObsData
  use input_params, fpar => fwat_par_global
  use decon_mod, only: deconit
  use fk_coupling
  use fwat_mpi
  use utils, only: zeros_dp, zeros, interp1
  use sacio
  use logger, only: log
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, nrec, nrec_local,OUTPUT_FILES, &
                        number_receiver_global, ispec_selected_rec, islice_selected_rec

  implicit none

  integer :: ier
  type, extends(SynData) :: RFData
    real(kind=cr), dimension(:), pointer :: ttp
    real(kind=dp), dimension(:), allocatable :: tstart, tend
    real(kind=cr) :: baz, az
    integer :: ttp_win
    contains
    procedure :: semd2sac
    procedure, private :: calc_times
  end type RFData

  contains

  subroutine semd2sac(this, ievt)
    class(RFData), intent(inout) :: this
    integer, intent(in) :: ievt
    real(kind=cr), dimension(:), allocatable :: bazi
    real(kind=dp), dimension(:), allocatable :: uin, win, rfi
    type(sachead) :: header
    integer :: irec_local, irec, igauss
    character(len=MAX_STRING_LEN) :: datafile, bandname

    call this%init(ievt)
    call this%od%read_stations(ievt, .true.)
    call this%calc_times()

    bazi = zeros(nrec)+this%baz
    call this%read(bazi)

    if (worldrank == 0) call system('mkdir -p '//trim(fpar%acqui%in_dat_path(this%ievt)))
    call synchronize_all()

    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        ! calculate rf
        uin = this%data(:, 2, irec)
        win = this%data(:, 3, irec)
        call bandpass_dp(uin ,NSTEP, dble(DT),1/fpar%sim%LONG_P(1),&
                         1/fpar%sim%SHORT_P(1), IORD)
        call bandpass_dp(win ,NSTEP, dble(DT),1/fpar%sim%LONG_P(1),&
                         1/fpar%sim%SHORT_P(1), IORD)
        do igauss = 1, fpar%sim%rf%NGAUSS
          write(bandname,'(a1,f3.1)') 'F',fpar%sim%rf%f0(igauss)
          datafile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%od%netwk(irec))//'.'&
                     //trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//'R.'&
                     //trim(bandname)//'.rf.sac'
          call deconit(uin, win, real(DT), fpar%sim%rf%tshift, fpar%sim%rf%f0(igauss),&
                       fpar%sim%rf%maxit, fpar%sim%rf%minderr, 0, rfi)
          call sacio_newhead(header, real(fpar%sim%dt), fpar%sim%nstep, -real(T0))
          header%az = this%az
          header%baz = this%baz
          header%stla = this%od%stla(irec)
          header%stlo = this%od%stlo(irec)
          header%stel = this%od%stel(irec)
          header%knetwk = this%od%netwk(irec)
          header%kstnm = this%od%stnm(irec)
          header%kcmpnm = trim(fpar%sim%CH_CODE)//'R'
          header%user1 = fpar%sim%rf%f0(igauss)
          header%kuser1 = 'gauss'
          call sacio_writesac(datafile, header, rfi, ier)
        enddo
      enddo
    endif
    call synchronize_all()

  end subroutine semd2sac

  subroutine calc_times(this)
    class(RFData), intent(inout) :: this
    real(kind=cr), dimension(:), allocatable :: ttp_local
    integer :: irec_local, irec

    ttp_local = zeros(this%nrec)
    call prepare_shm_array_cr_1d(this%ttp, this%nrec, this%ttp_win)
    call read_fk_model(fpar%acqui%evtid_names(this%ievt))

    this%baz = -phi_FK - 90.d0
    this%az = 90.d0 - phi_FK
    
    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        ttp_local(irec) = maxloc(seismograms_d(3, irec_local, 1:NSTEP), dim=1) * DT - T0
      end do
    endif
    call sum_all_1Darray_cr(ttp_local, this%ttp, this%nrec)
    call sync_from_main_rank_cr_1d(this%ttp, this%nrec)
    call free_fk_arrays()

  end subroutine calc_times

end module rf_data