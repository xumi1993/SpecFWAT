module misfit_mod
  use config

  implicit none

  type :: RECMisfit
    character(len=MAX_STRING_LEN) :: net, sta
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: chan, band_name
    real(kind=dp), dimension(:,:), allocatable :: tstart, tend
    real(kind=dp), dimension(:,:,:), allocatable :: misfits, residuals, imeas
    real(kind=dp), dimension(:), allocatable :: total_misfit
    integer :: nwin
  end type RECMisfit

  interface RECMisfit
    module procedure initialize_recm
  end interface RECMisfit

  type :: EVTMisfit
    integer :: nflt, nrec
    type(RECMisfit), dimension(:), allocatable :: rec_misfits
    real(kind=dp), dimension(:), allocatable :: total_misfit
    character(len=MAX_STRING_LEN) :: evtname
    contains
      procedure :: assemble => assemble_misfit_from_all_rec, write
  end type EVTMisfit

  interface EVTMisfit
    module procedure initialize_evtm
  end interface EVTMisfit

contains
  function initialize_recm(ncomp, nwin, nflt) result(recm)
    ! Initialize the RECMisfit object
    ! ncomp: number of components (e.g., 2 for radial and vertical)
    ! nwin: number of measurement windows
    ! nflt: number of frequency bands (Gaussian filters)
    integer, intent(in) :: nwin, ncomp, nflt
    type(RECMisfit) :: recm
    
    recm%nwin = nwin
    recm%net = ''
    recm%sta = ''
    allocate(recm%tstart(recm%nwin, nflt))
    allocate(recm%tend(recm%nwin, nflt))
    allocate(recm%misfits(ncomp, recm%nwin, nflt))
    allocate(recm%residuals(ncomp, recm%nwin, nflt))
    allocate(recm%imeas(ncomp, recm%nwin, nflt))
    allocate(recm%total_misfit(nflt))
    allocate(recm%chan(ncomp))
    allocate(recm%band_name(nflt))
    recm%imeas = 0
    recm%tstart = 0.0_dp
    recm%tend = 0.0_dp
    recm%misfits = 0.0_dp
    recm%residuals = 0.0_dp
    recm%total_misfit = 0.0_dp

  end function initialize_recm

  function initialize_evtm(evtname, nrec, nflt) result(evtm)
    ! Initialize the EVTMisfit object
    ! nrec: number of records
    ! nflt: number of frequency bands (Gaussian filters)
    integer, intent(in) :: nrec, nflt
    character(len=*), intent(in) :: evtname
    type(EVTMisfit) :: evtm

    evtm%evtname = trim(evtname)
    evtm%nrec = nrec
    evtm%nflt = nflt
    allocate(evtm%rec_misfits(nrec))
    allocate(evtm%total_misfit(nflt))
    evtm%total_misfit = 0.0_dp

  end function initialize_evtm

  subroutine assemble_misfit_from_all_rec(this, misfit_loc)
    ! Assemble misfit from all local RECMisfit objects to global arrays
    use fwat_mpi
    class(EVTMisfit), intent(inout) :: this
    type(RECMisfit), dimension(:), allocatable, intent(in) :: misfit_loc
    integer :: nrec_loc, ncomp, nwin, nflt, nsta_irank
    integer :: irec_loc, irec_glob, iproc, i
    integer, dimension(:), allocatable :: recv_indices, send_indices, recv_nwin
    real(kind=dp), dimension(:,:), allocatable :: recv_tstart, recv_tend, recv_total_misfit
    real(kind=dp), dimension(:,:,:), allocatable :: recv_misfits, recv_residuals, recv_imeas
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: recv_net, recv_sta
    character(len=MAX_STRING_LEN), dimension(:,:), allocatable :: recv_chan, recv_band_name
    
    ncomp = size(misfit_loc(1)%misfits, 1)
    nflt = this%nflt
    if (allocated(misfit_loc)) then
      nrec_loc = size(misfit_loc)
    else
      nrec_loc = 0
    endif

    if (worldrank == 0) then
      ! First copy local data on rank 0
      if (nrec_loc > 0) then
        do irec_loc = 1, nrec_loc
          irec_glob = select_global_id_for_rec(irec_loc)
          nwin = misfit_loc(irec_loc)%nwin
          this%rec_misfits(irec_glob) = RECMisfit(ncomp, nwin, nflt)
          this%rec_misfits(irec_glob)%net = misfit_loc(irec_loc)%net
          this%rec_misfits(irec_glob)%sta = misfit_loc(irec_loc)%sta
          this%rec_misfits(irec_glob)%chan = misfit_loc(irec_loc)%chan
          this%rec_misfits(irec_glob)%band_name = misfit_loc(irec_loc)%band_name
          this%rec_misfits(irec_glob)%tstart = misfit_loc(irec_loc)%tstart
          this%rec_misfits(irec_glob)%tend = misfit_loc(irec_loc)%tend
          this%rec_misfits(irec_glob)%misfits = misfit_loc(irec_loc)%misfits
          this%rec_misfits(irec_glob)%residuals = misfit_loc(irec_loc)%residuals
          this%rec_misfits(irec_glob)%imeas = misfit_loc(irec_loc)%imeas
          this%rec_misfits(irec_glob)%total_misfit = misfit_loc(irec_loc)%total_misfit
        enddo
      endif
      
      ! Then receive data from other ranks
      do iproc = 1, worldsize-1
        nsta_irank = get_num_recs_per_proc(this%nrec, iproc)
        if (nsta_irank > 0) then
          ! Receive indices and nwin values
          allocate(recv_indices(nsta_irank))
          allocate(recv_nwin(nsta_irank))
          allocate(recv_net(nsta_irank))
          allocate(recv_sta(nsta_irank))
          allocate(recv_chan(nsta_irank, ncomp))
          allocate(recv_band_name(nsta_irank, nflt))
          allocate(recv_total_misfit(nsta_irank, nflt))
          
          call recv_i(recv_indices, nsta_irank, iproc, targ)
          call recv_i(recv_nwin, nsta_irank, iproc, targ)
          call recv_ch_array(recv_net, nsta_irank, MAX_STRING_LEN, iproc, targ)
          call recv_ch_array(recv_sta, nsta_irank, MAX_STRING_LEN, iproc, targ)
          call recv_ch_array(recv_chan, nsta_irank*ncomp, MAX_STRING_LEN, iproc, targ)
          call recv_ch_array(recv_band_name, nsta_irank*nflt, MAX_STRING_LEN, iproc, targ)
          call recv_dp(recv_total_misfit, nsta_irank*nflt, iproc, targ)
          
          ! Receive data for each record
          do i = 1, nsta_irank
            irec_glob = recv_indices(i)
            nwin = recv_nwin(i)
            
            this%rec_misfits(irec_glob) = RECMisfit(ncomp, nwin, nflt)
            this%rec_misfits(irec_glob)%net = recv_net(i)
            this%rec_misfits(irec_glob)%sta = recv_sta(i)
            this%rec_misfits(irec_glob)%chan = recv_chan(i, :)
            this%rec_misfits(irec_glob)%band_name = recv_band_name(i, :)
            this%rec_misfits(irec_glob)%total_misfit = recv_total_misfit(i, :)
            
            ! Allocate and receive arrays for this record
            allocate(recv_tstart(nwin, nflt))
            allocate(recv_tend(nwin, nflt))
            allocate(recv_misfits(ncomp, nwin, nflt))
            allocate(recv_residuals(ncomp, nwin, nflt))
            allocate(recv_imeas(ncomp, nwin, nflt))
            
            call recv_dp(recv_tstart, nwin*nflt, iproc, targ)
            call recv_dp(recv_tend, nwin*nflt, iproc, targ)
            call recv_dp(recv_misfits, ncomp*nwin*nflt, iproc, targ)
            call recv_dp(recv_residuals, ncomp*nwin*nflt, iproc, targ)
            call recv_dp(recv_imeas, ncomp*nwin*nflt, iproc, targ)
            
            ! Direct assignment - no reshape needed
            this%rec_misfits(irec_glob)%tstart = recv_tstart
            this%rec_misfits(irec_glob)%tend = recv_tend
            this%rec_misfits(irec_glob)%misfits = recv_misfits
            this%rec_misfits(irec_glob)%residuals = recv_residuals
            this%rec_misfits(irec_glob)%imeas = recv_imeas
            
            deallocate(recv_tstart, recv_tend, recv_misfits, recv_residuals, recv_imeas)
          enddo
          
          deallocate(recv_indices, recv_nwin, recv_net, recv_sta, recv_chan, recv_band_name, recv_total_misfit)
        end if ! nsta_irank > 0
      enddo ! iproc
      
    else ! worldrank != 0
      if (nrec_loc > 0) then
        allocate(send_indices(nrec_loc))
        allocate(recv_nwin(nrec_loc))
        allocate(recv_net(nrec_loc))
        allocate(recv_sta(nrec_loc))
        allocate(recv_chan(nrec_loc, ncomp))
        allocate(recv_band_name(nrec_loc, nflt))
        allocate(recv_total_misfit(nrec_loc, nflt))
        
        ! Prepare indices and basic data for each record
        do irec_loc = 1, nrec_loc
          send_indices(irec_loc) = select_global_id_for_rec(irec_loc)
          recv_nwin(irec_loc) = misfit_loc(irec_loc)%nwin
          recv_net(irec_loc) = misfit_loc(irec_loc)%net
          recv_sta(irec_loc) = misfit_loc(irec_loc)%sta
          recv_chan(irec_loc, :) = misfit_loc(irec_loc)%chan
          recv_band_name(irec_loc, :) = misfit_loc(irec_loc)%band_name
          recv_total_misfit(irec_loc, :) = misfit_loc(irec_loc)%total_misfit
        enddo
        
        ! Send basic information first
        call send_i(send_indices, nrec_loc, 0, targ)
        call send_i(recv_nwin, nrec_loc, 0, targ)
        call send_ch_array(recv_net, nrec_loc, MAX_STRING_LEN, 0, targ)
        call send_ch_array(recv_sta, nrec_loc, MAX_STRING_LEN, 0, targ)
        call send_ch_array(recv_chan, nrec_loc*ncomp, MAX_STRING_LEN, 0, targ)
        call send_ch_array(recv_band_name, nrec_loc*nflt, MAX_STRING_LEN, 0, targ)
        call send_dp(recv_total_misfit, nrec_loc*nflt, 0, targ)
        
        ! Send array data for each record
        do irec_loc = 1, nrec_loc
          nwin = misfit_loc(irec_loc)%nwin
          
          ! Send arrays directly - Fortran will handle memory layout
          call send_dp(misfit_loc(irec_loc)%tstart, nwin*nflt, 0, targ)
          call send_dp(misfit_loc(irec_loc)%tend, nwin*nflt, 0, targ)
          call send_dp(misfit_loc(irec_loc)%misfits, ncomp*nwin*nflt, 0, targ)
          call send_dp(misfit_loc(irec_loc)%residuals, ncomp*nwin*nflt, 0, targ)
          call send_dp(misfit_loc(irec_loc)%imeas, ncomp*nwin*nflt, 0, targ)
        enddo
        
        deallocate(send_indices, recv_nwin, recv_net, recv_sta, recv_chan, recv_band_name, recv_total_misfit)
      endif
    endif

    ! Calculate total misfit across all ranks
    this%total_misfit = 0.0_dp
    if (worldrank == 0) then
      do irec_glob = 1, this%nrec
        this%total_misfit = this%total_misfit + this%rec_misfits(irec_glob)%total_misfit
      enddo
    endif

  end subroutine assemble_misfit_from_all_rec

  subroutine write(this)
    class(EVTMisfit), intent(inout) :: this
    integer :: unit, ier
    integer :: irec, icomp, iwin, iflt, ncomp
    character(len=MAX_STRING_LEN) :: chifile

    if (worldrank /= 0) return

    ncomp = size(this%rec_misfits(1)%misfits, 1)
    do iflt = 1, this%nflt
      chifile = trim(MISFITS_DIR)//'/'//trim(model_name)//'.'//&
                trim(this%evtname)//'_'//trim(this%rec_misfits(1)%band_name(iflt))//&
                '_window_chi'
      open(newunit=unit, file=chifile, status='unknown', iostat=ier)
      if (ier /= 0) then
        write(*,*) 'Error opening file: ', trim(chifile)
        error stop
      end if

      do irec = 1, this%nrec
        do icomp = 1, ncomp
          do iwin = 1, this%rec_misfits(irec)%nwin
            write(unit, '(a,1x,a,1x,a,1x,i4,1x,f8.3,1x,f8.3,1x,f12.6,1x,f12.6)') &
                         trim(this%rec_misfits(irec)%net), &
                         trim(this%rec_misfits(irec)%sta), &
                         trim(this%rec_misfits(irec)%chan(icomp)), &
                         this%rec_misfits(irec)%imeas(icomp, iwin, iflt), &
                         this%rec_misfits(irec)%tstart(iwin, iflt), &
                         this%rec_misfits(irec)%tend(iwin, iflt), &
                         this%rec_misfits(irec)%residuals(icomp, iwin, iflt), &
                         this%rec_misfits(irec)%misfits(icomp, iwin, iflt)
          enddo ! iwin
        enddo ! icomp
      enddo ! irec
      close(unit)
    enddo ! iflt

  end subroutine write

  function read_evt_misfit(model_name_in, ievt) result(total_misfit)
    use input_params, only: fpar => fwat_par_global
    use common_lib, only: get_band_name
    use utils, only: split_by_spaces

    character(len=*), intent(in) :: model_name_in
    integer, intent(in) :: ievt
    character(len=MAX_STRING_LEN) :: chifile, band_name, line
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: line_sp
    real(kind=dp) :: total_misfit
    real(kind=dp) :: misfit
    integer :: unit, ier, iflt, nflt

    if (dat_type == 'rf') then
      nflt = fpar%sim%rf%NGAUSS
      write(band_name, '("F",F3.1)') fpar%sim%rf%f0(iflt)
    else
      nflt = fpar%sim%NUM_FILTER
      call get_band_name(fpar%sim%SHORT_P(iflt), fpar%sim%LONG_P(iflt), band_name)
    end if

    total_misfit = 0.0_dp
    if (worldrank == 0) then
      do iflt = 1, nflt
        chifile = trim(MISFITS_DIR)//'/'//trim(model_name_in)//'.'//&
                  trim(fpar%acqui%evtid_names(ievt))//'_'//trim(band_name)//'_window_chi'
        open(newunit=unit, file=chifile, status='old', iostat=ier)
        if (ier /= 0) then
          write(*,*) 'Error opening file: ', trim(chifile)
          error stop
        end if
        do
          read(unit, '(a)', iostat=ier) line
          line_sp = split_by_spaces(trim(line))
          read(line_sp(8), *) misfit
          if (ier /= 0) exit
          total_misfit = total_misfit + misfit
        end do
      end do
    endif
    call synchronize_all()
    call bcast_all_singledp(total_misfit)

  end function read_evt_misfit

end module misfit_mod