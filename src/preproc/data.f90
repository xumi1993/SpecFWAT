module obs_data
  use config
  use fwat_mpi
  use common_lib, only: find_string
  use imput_params, fpar => fwat_par_global

  implicit none

  integer, private :: ier

  type evt_data
    character(len=MAX_STRING_LEN), dimension(:), pointer :: netwk, stnm
    character(len=MAX_STRING_LEN) :: evt_id
    real(kind=cr), dimension(:), pointer :: stla, stlo, stel
    integer :: nsta
    real(kind=dp), dimension(:, :, :), pointer :: data
    integer :: net_win, sta_win, stla_win, stlo_win, stel_win, dat_win
    contains
    procedure :: read_stations
    procedure, private :: alloc_sta_info
  end type evt_data

  type(evt_data) :: fwat_evt_data_global

contains
  subroutine read_stations(this, evt_id)
    class(evt_data), intent(inout) :: this
    character(len=MAX_STRING_LEN), intent(in) :: evt_id
    character(len=MAX_STRING_LEN) :: line
    integer :: idx, ista
    integer, parameter :: FID = 868

    this%evt_id = evt_id
    idx = find_string(fpar%acqui%evtid_names, evt_id)
    if (idx == -1) then
      if (worldrank == 0) call exit_MPI(0, 'Event ID not found')
    endif
    if (worldrank == 0) then
      open(unit=FID, file=fpar%acqui%station_file(idx), status='old', iostat=ier)
      if (ier /= 0) call exit_MPI(0, 'Error opening station file'//trim(fpar%acqui%station_file(idx)))
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
    class(evt_data), intent(inout) :: this

    call prepare_shm_array_ch_1d(this%netwk, this%nsta, MAX_STRING_LEN, this%net_win)
    call prepare_shm_array_ch_1d(this%stnm, this%nsta, MAX_STRING_LEN, this%sta_win)
    call prepare_shm_array_cr_1d(this%stla, this%nsta, this%stla_win)
    call prepare_shm_array_cr_1d(this%stlo, this%nsta, this%stlo_win)
    call prepare_shm_array_cr_1d(this%stel, this%nsta, this%stel_win)
    call prepare_shm_array_dp_3d(this%data, fpar%sim%NSTEP, fpar%sim%NRCOMP, this%nsta, this%dat_win)
  end subroutine alloc_sta_info
end module obs_data