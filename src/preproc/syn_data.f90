module syn_data
  use config
  use fwat_constants, only: targ
  use specfem_par, only: NSTEP, DT, ELASTIC_SIMULATION, seismograms_d, nrec, nrec_local, &
                        number_receiver_global, ispec_selected_rec, islice_selected_rec
  use specfem_par_elastic, only: ispec_is_elastic
  use fwat_mpi
  use input_params, fpar => fwat_par_global
  use common_lib, only: rotate_NE_to_RT

  implicit none

  integer, private :: ier
  integer, parameter, private :: NCOMP=3

  type :: SynData
    integer :: ievt, nrec
    real(kind=dp), dimension(:, :, :), pointer :: data ! npts, ncomp(zrt), nsta
    character(len=MAX_STRING_LEN) :: band_name
    integer :: dat_win
    contains
    procedure :: read=>read_syn_data, rotate_to_rt, filter
    procedure :: finalize

  end type SynData

contains

  subroutine read_syn_data(this, ievt)
    class(SynData), intent(inout) :: this
    integer, intent(in) :: ievt
    integer :: irec, irec_local, ispec, iproc, nsta_irank
    real(kind=dp), dimension(:, :, :), allocatable :: data_local, dat_sum

    ! read source and reciever files
    this%ievt = ievt
    this%nrec = nrec
    call prepare_shm_array_dp_3d(this%data, NSTEP, NCOMP, nrec, this%dat_win)
    allocate(data_local(NSTEP, NCOMP, nrec_local))
    data_local = 0.0_dp

    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        ispec = ispec_selected_rec(irec)
        if (ELASTIC_SIMULATION) then
          if (ispec_is_elastic(ispec)) then
            data_local(:, 1, irec_local) = dble(seismograms_d(3, irec_local, 1:NSTEP))
            data_local(:, 2, irec_local) = dble(seismograms_d(2, irec_local, 1:NSTEP))
            data_local(:, 3, irec_local) = dble(seismograms_d(1, irec_local, 1:NSTEP))
          endif
        endif
      enddo
    endif

    ! collect data to rank 0
    if (worldrank == 0) then
      if (nrec_local > 0) then
        do irec_local = 1, nrec_local
          irec = number_receiver_global(irec_local)
          this%data(:, :, irec) = data_local(:, :, irec_local)
        enddo
      endif
      do iproc = 1, worldsize-1
        nsta_irank=0
        do irec = 1, nrec
          if (islice_selected_rec(irec) == iproc) nsta_irank = nsta_irank + 1
        enddo
        if (nsta_irank > 0) then
          do irec_local = 1, nsta_irank
            call recv_singlei(irec, 0, targ)
            call recv_dp(this%data(:, :, irec), NSTEP*NCOMP, iproc, targ)
          enddo
        endif
      enddo
    else
      if (nrec_local > 0) then
        do irec_local = 1, nrec_local
          call send_singlei(number_receiver_global(irec_local), 0, targ)
          call send_dp(data_local(:, :, irec_local), NSTEP*NCOMP, 0, targ)
        enddo
      endif
    endif
    call sync_from_main_rank_dp_3d(this%data, NSTEP, NCOMP, nrec)

  end subroutine read_syn_data

  subroutine rotate_to_rt(this, baz)
    class(SynData), intent(inout) :: this
    real(kind=cr), dimension(:), intent(in) :: baz
    real(kind=cr), dimension(:), allocatable :: data_r, data_t
    integer :: irec, it

    allocate(data_r(NSTEP), data_t(NSTEP))
    if (noderank == 0) then
      do irec = 1, nrec
        call rotate_NE_to_RT(real(this%data(:, 2, irec)), real(this%data(:, 3, irec)), data_r, data_t, baz(irec))
        this%data(:, 2, irec) = dble(data_r)
        this%data(:, 3, irec) = dble(data_t)
      enddo
    endif
    call synchronize_all()

  end subroutine rotate_to_rt

  subroutine filter(this, freqmin, freqmax, order)
    use signal, only: bandpass_dp

    class(SynData), intent(inout) :: this
    real(kind=cr), intent(in) :: freqmin, freqmax
    integer, intent(in), optional :: order
    integer :: irec, icomp, order_local

    if (present(order)) then
      order_local = order
    else
      order_local = 2
    endif

    if (noderank == 0) then
      do irec = 1, nrec
        do icomp = 1, NCOMP
          call bandpass_dp(this%data(:, icomp, irec), NSTEP, DT, freqmin, freqmax, order_local)
        enddo
      enddo
    endif
    call synchronize_all()

  end subroutine filter

  subroutine finalize(this)
    class(SynData), intent(inout) :: this

    call free_shm_array(this%dat_win)

  end subroutine finalize

end module syn_data