module obs_data
  use config
  use fwat_mpi
  use common_lib, only: find_string
  use sacio
  use input_params, fpar => fwat_par_global

  implicit none

  integer, private :: ier

  type ObsData
    character(len=MAX_STRING_LEN), dimension(:), pointer :: netwk, stnm
    character(len=MAX_STRING_LEN) :: evt_id
    real(kind=cr), dimension(:), pointer :: stla, stlo, stel, baz, tarr, tbeg
    real(kind=dp) :: dt
    integer :: nsta, npts, ievt
    real(kind=dp), dimension(:, :, :), pointer :: data ! npts, ncomp, nsta
    integer :: net_win, sta_win, stla_win, stlo_win, stel_win, dat_win, &
                baz_win, t0_win, tb_win
    contains
    procedure :: read_stations, read_obs_data
    procedure, private :: alloc_sta_info
    procedure :: finalize
  end type ObsData

  type(ObsData) :: fwat_evt_data_global

contains
  subroutine read_stations(this, ievt, is_filted)
    class(ObsData), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: evt_id
    integer, intent(in) :: ievt
    logical, intent(in), optional :: is_filted
    character(len=MAX_STRING_LEN) :: line, fname
    integer :: idx, ista
    integer, parameter :: FID = 868

    this%ievt = ievt
    this%evt_id = fpar%acqui%evtid_names(ievt)

    if (present(is_filted)) then
      fname = fpar%acqui%station_file(ievt)//'_FILTERED'
    else
      fname = fpar%acqui%station_file(ievt)
    endif

    if (worldrank == 0) then
      open(unit=FID, file=fname, status='old', iostat=ier)
      if (ier /= 0) call exit_MPI(0, 'Error opening station file'//trim(fname))
      ista = 0
      do
        read(FID, *, iostat=ier) line
        if (ier /= 0) exit
        ista = ista + 1
      enddo
      this%nsta = ista
      rewind(FID)
    endif
    call bcast_all_singlei(this%nsta)
    call this%alloc_sta_info()

    if (worldrank == 0) then
      do ista = 1, this%nsta
        read(FID, *, iostat=ier) this%stnm(ista), this%netwk(ista), this%stla(ista), this%stlo(ista), this%stel(ista), line
      enddo
      close(FID)
    endif
    call sync_from_main_rank_ch(this%netwk, this%nsta, MAX_STRING_LEN)
    call sync_from_main_rank_ch(this%stnm, this%nsta, MAX_STRING_LEN)
    call sync_from_main_rank_cr_1d(this%stla, this%nsta)
    call sync_from_main_rank_cr_1d(this%stlo, this%nsta)
    call sync_from_main_rank_cr_1d(this%stel, this%nsta)

  end subroutine read_stations

  subroutine alloc_sta_info(this)
    class(ObsData), intent(inout) :: this

    call prepare_shm_array_ch_1d(this%netwk, this%nsta, MAX_STRING_LEN, this%net_win)
    call prepare_shm_array_ch_1d(this%stnm, this%nsta, MAX_STRING_LEN, this%sta_win)
    call prepare_shm_array_cr_1d(this%stla, this%nsta, this%stla_win)
    call prepare_shm_array_cr_1d(this%stlo, this%nsta, this%stlo_win)
    call prepare_shm_array_cr_1d(this%stel, this%nsta, this%stel_win)
    call prepare_shm_array_cr_1d(this%baz, this%nsta, this%baz_win)
    call prepare_shm_array_cr_1d(this%tbeg, this%nsta, this%tb_win)
    if (index(dat_type, 'tele') /= 0) &
      call prepare_shm_array_cr_1d(this%tarr, this%nsta, this%t0_win)
  end subroutine alloc_sta_info

  subroutine finalize(this)
    class(ObsData), intent(inout) :: this

    call free_shm_array(this%net_win)
    call free_shm_array(this%sta_win)
    call free_shm_array(this%stla_win)
    call free_shm_array(this%stlo_win)
    call free_shm_array(this%stel_win)
    call free_shm_array(this%baz_win)
    call free_shm_array(this%t0_win)
    call free_shm_array(this%dat_win)
    call free_shm_array(this%tb_win)

  end subroutine finalize

  subroutine read_obs_data(this)
    class(ObsData), intent(inout) :: this
    integer :: ista, icomp, ncomp
    type(sachead) :: header
    character(len=MAX_STRING_LEN) :: sacfile, gaus_str
    real(kind=dp), dimension(:), allocatable :: datarray

    if (worldrank == 0) then
      ! read sample to get header info
      if (dat_type /= 'rf') then
        sacfile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%netwk(ista))//'.'//trim(this%stnm(ista))//&
                  '.'//trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))//'.sac'
        ncomp = fpar%sim%NRCOMP
      else
        write(gaus_str, '("F",F3.1)') fpar%sim%rf%F0(1)
        sacfile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%netwk(ista))//'.'//trim(this%stnm(ista))//&
                  '.'//trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))//'.'//trim(gaus_str)//'.rf.sac'
        ncomp = fpar%sim%rf%NGAUSS
      endif
      call sacio_readhead(sacfile, header, ier)
      if (ier /= 0) call exit_MPI(0, 'Error reading SAC header '//trim(sacfile))
      this%npts = header%npts
      this%dt = dble(header%delta)
    endif
    call bcast_all_singlei(this%npts)
    call bcast_all_singledp(this%dt)
    call bcast_all_singlei(ncomp)
    call prepare_shm_array_dp_3d(this%data, this%npts, ncomp, this%nsta, this%dat_win)

    if (worldrank == 0) then
      do ista = 1, this%nsta
        do icomp = 1, ncomp
          ! get sac file name
          if (dat_type /= 'rf') then
            sacfile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%netwk(ista))//'.'//trim(this%stnm(ista))//&
                      '.'//trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))//'.sac'
          else
            write(gaus_str, '("F",F3.1)') fpar%sim%rf%F0(icomp)
            sacfile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%netwk(ista))//'.'//trim(this%stnm(ista))//&
                      '.'//trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))//'.'//trim(gaus_str)//'.rf.sac'
          endif

          ! read sac file
          call sacio_readsac(sacfile, header, datarray, ier)
          if (ier /= 0) call exit_MPI(0, 'Error reading SAC file '//trim(sacfile))

          ! assign data to this%data
          this%data(:, icomp, ista) = datarray

          ! assign header info
          if (icomp == 1) then
            this%tbeg(ista) = header%b
            if (header%baz == SAC_rnull) then
              call exit_MPI(0, 'Back azimuth not found in SAC header')
            else
              this%baz(ista) = header%baz
            endif
            if (index(dat_type, 'tele') /= 0 ) then
              if (header%t0 /= SAC_rnull) then
                this%tarr(ista) = header%t0
              else
                call exit_MPI(0, 't0 not found in SAC header')
              endif
            endif
          endif
        enddo
      enddo
    endif

    call sync_from_main_rank_cr_1d(this%baz, this%nsta)
    call sync_from_main_rank_cr_1d(this%tarr, this%nsta)
    call sync_from_main_rank_cr_1d(this%tbeg, this%nsta)
    call sync_from_main_rank_dp_3d(this%data, this%npts, ncomp, this%nsta)
    call synchronize_all()

  end subroutine read_obs_data
  
end module obs_data