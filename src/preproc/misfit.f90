module misfit_mod
  use config

  implicit none

  type TraceMisfit
    real(kind=dp), dimension(:), allocatable :: misfits, residuals, tstart, tend
    integer, dimension(:), allocatable :: imeas
    integer :: nwin
  end type TraceMisfit

  interface TraceMisfit
    module procedure initialize_trm
  end interface TraceMisfit

  type :: RECMisfit
    character(len=MAX_STRING_LEN) :: net, sta
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: chan
    type(TraceMisfit), dimension(:), allocatable :: trm
    ! real(kind=dp), dimension(:,:), allocatable :: tstart, tend
    ! real(kind=dp), dimension(:,:,:), allocatable :: misfits, residuals
    ! integer, dimension(:,:,:), allocatable :: imeas
    real(kind=dp) :: total_misfit
    integer :: ncomp
  end type RECMisfit

  interface RECMisfit
    module procedure initialize_recm
  end interface RECMisfit

  type :: EVTMisfit
    integer :: nrec
    type(RECMisfit), dimension(:), allocatable :: rec_misfits
    real(kind=dp) :: total_misfit
    character(len=MAX_STRING_LEN) :: evtname, band_name
    contains
      procedure :: assemble => assemble_misfit_from_all_rec, write
  end type EVTMisfit

  interface EVTMisfit
    module procedure initialize_evtm
  end interface EVTMisfit

contains
  function initialize_trm(nwin) result(trm)
    ! Initialize the TraceMisfit object
    integer, intent(in) :: nwin
    type(TraceMisfit) :: trm

    trm%nwin = nwin
    allocate(trm%misfits(nwin))
    allocate(trm%residuals(nwin))
    allocate(trm%tstart(nwin))
    allocate(trm%tend(nwin))
    allocate(trm%imeas(nwin))

    trm%misfits = 0.0_dp
    trm%residuals = 0.0_dp
    trm%tstart = 0.0_dp
    trm%tend = 0.0_dp
    trm%imeas = 0

  end function initialize_trm

  function initialize_recm(ncomp) result(recm)
    ! Initialize the RECMisfit object
    ! ncomp: number of components (e.g., 2 for radial and vertical)
    integer, intent(in) :: ncomp
    type(RECMisfit) :: recm
    
    recm%ncomp = ncomp
    recm%net = ''
    recm%sta = ''
    ! allocate(recm%tstart(recm%nwin, nflt))
    ! allocate(recm%tend(recm%nwin, nflt))
    ! allocate(recm%misfits(ncomp, recm%nwin, nflt))
    ! allocate(recm%residuals(ncomp, recm%nwin, nflt))
    ! allocate(recm%imeas(ncomp, recm%nwin, nflt))
    ! allocate(recm%total_misfit(nflt))
    allocate(recm%chan(ncomp))
    allocate(recm%trm(ncomp))
    ! allocate(recm%band_name(nflt))
    ! recm%imeas = 0
    ! recm%tstart = 0.0_dp
    ! recm%tend = 0.0_dp
    ! recm%misfits = 0.0_dp
    ! recm%residuals = 0.0_dp
    recm%total_misfit = 0.0_dp

  end function initialize_recm

  function initialize_evtm(evtname, nrec) result(evtm)
    ! Initialize the EVTMisfit object
    ! nrec: number of records
    integer, intent(in) :: nrec
    character(len=*), intent(in) :: evtname
    type(EVTMisfit) :: evtm

    evtm%evtname = trim(evtname)
    evtm%nrec = nrec
    if (worldrank == 0) allocate(evtm%rec_misfits(nrec))
    evtm%total_misfit = 0.0_dp

  end function initialize_evtm

  subroutine assemble_misfit_from_all_rec(this, misfit_loc)
    ! Assemble misfit from all local RECMisfit objects to global arrays
    use fwat_mpi
    class(EVTMisfit), intent(inout) :: this
    type(RECMisfit), dimension(:), allocatable, intent(in) :: misfit_loc
    integer :: nrec_loc, ncomp, nwin, nsta_irank
    integer :: irec_loc, irec_glob, iproc, i, icomp
    integer, dimension(:), allocatable :: recv_indices, send_indices, recv_ncomp, recv_nwin_arr
    real(kind=dp), dimension(:), allocatable :: recv_total_misfit
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: recv_net, recv_sta
    character(len=MAX_STRING_LEN), dimension(:,:), allocatable :: recv_chan
    ! Temporary arrays for each TraceMisfit
    real(kind=dp), dimension(:), allocatable :: temp_misfits, temp_residuals, temp_tstart, temp_tend
    integer, dimension(:), allocatable :: temp_imeas
    
    nrec_loc = 0
    ncomp = 0
    if (allocated(misfit_loc)) then
      nrec_loc = size(misfit_loc)
      if (nrec_loc > 0) ncomp = misfit_loc(1)%ncomp
    endif

    if (worldrank == 0) then
      ! First copy local data on rank 0
      if (nrec_loc > 0) then
        do irec_loc = 1, nrec_loc
          irec_glob = select_global_id_for_rec(irec_loc)
          this%rec_misfits(irec_glob) = RECMisfit(misfit_loc(irec_loc)%ncomp)
          this%rec_misfits(irec_glob)%net = misfit_loc(irec_loc)%net
          this%rec_misfits(irec_glob)%sta = misfit_loc(irec_loc)%sta
          this%rec_misfits(irec_glob)%chan = misfit_loc(irec_loc)%chan
          this%rec_misfits(irec_glob)%total_misfit = misfit_loc(irec_loc)%total_misfit
          
          ! Copy TraceMisfit objects for each component
          do icomp = 1, misfit_loc(irec_loc)%ncomp
            nwin = misfit_loc(irec_loc)%trm(icomp)%nwin
            this%rec_misfits(irec_glob)%trm(icomp) = TraceMisfit(nwin)
            this%rec_misfits(irec_glob)%trm(icomp)%misfits = misfit_loc(irec_loc)%trm(icomp)%misfits
            this%rec_misfits(irec_glob)%trm(icomp)%residuals = misfit_loc(irec_loc)%trm(icomp)%residuals
            this%rec_misfits(irec_glob)%trm(icomp)%tstart = misfit_loc(irec_loc)%trm(icomp)%tstart
            this%rec_misfits(irec_glob)%trm(icomp)%tend = misfit_loc(irec_loc)%trm(icomp)%tend
            this%rec_misfits(irec_glob)%trm(icomp)%imeas = misfit_loc(irec_loc)%trm(icomp)%imeas
          enddo
        enddo
      endif
      
      ! Then receive data from other ranks
      do iproc = 1, worldsize-1
        nsta_irank = get_num_recs_per_proc(this%nrec, iproc)
        if (nsta_irank > 0) then
          ! Allocate arrays for basic info
          allocate(recv_indices(nsta_irank))
          allocate(recv_ncomp(nsta_irank))
          allocate(recv_net(nsta_irank))
          allocate(recv_sta(nsta_irank))
          allocate(recv_total_misfit(nsta_irank))
          
          ! Receive basic information
          call recv_i(recv_indices, nsta_irank, iproc, targ)
          call recv_i(recv_ncomp, nsta_irank, iproc, targ)
          call recv_ch_array(recv_net, nsta_irank, MAX_STRING_LEN, iproc, targ)
          call recv_ch_array(recv_sta, nsta_irank, MAX_STRING_LEN, iproc, targ)
          call recv_dp(recv_total_misfit, nsta_irank, iproc, targ)
          
          ! For each record, receive channel names and TraceMisfit data
          do i = 1, nsta_irank
            irec_glob = recv_indices(i)
            ncomp = recv_ncomp(i)
            
            ! Initialize the RECMisfit object
            this%rec_misfits(irec_glob) = RECMisfit(ncomp)
            this%rec_misfits(irec_glob)%net = recv_net(i)
            this%rec_misfits(irec_glob)%sta = recv_sta(i)
            this%rec_misfits(irec_glob)%total_misfit = recv_total_misfit(i)
            
            ! Receive channel names
            allocate(recv_chan(1, ncomp))
            call recv_ch_array(recv_chan, ncomp, MAX_STRING_LEN, iproc, targ)
            this%rec_misfits(irec_glob)%chan = recv_chan(1, :)
            deallocate(recv_chan)
            
            ! Receive nwin for each component
            allocate(recv_nwin_arr(ncomp))
            call recv_i(recv_nwin_arr, ncomp, iproc, targ)
            
            ! Receive TraceMisfit data for each component
            do icomp = 1, ncomp
              nwin = recv_nwin_arr(icomp)
              this%rec_misfits(irec_glob)%trm(icomp) = TraceMisfit(nwin)
              
              ! Receive arrays for this component
              allocate(temp_misfits(nwin))
              allocate(temp_residuals(nwin))
              allocate(temp_tstart(nwin))
              allocate(temp_tend(nwin))
              allocate(temp_imeas(nwin))
              
              call recv_dp(temp_misfits, nwin, iproc, targ)
              call recv_dp(temp_residuals, nwin, iproc, targ)
              call recv_dp(temp_tstart, nwin, iproc, targ)
              call recv_dp(temp_tend, nwin, iproc, targ)
              call recv_i(temp_imeas, nwin, iproc, targ)
              
              ! Copy to the TraceMisfit object
              this%rec_misfits(irec_glob)%trm(icomp)%misfits = temp_misfits
              this%rec_misfits(irec_glob)%trm(icomp)%residuals = temp_residuals
              this%rec_misfits(irec_glob)%trm(icomp)%tstart = temp_tstart
              this%rec_misfits(irec_glob)%trm(icomp)%tend = temp_tend
              this%rec_misfits(irec_glob)%trm(icomp)%imeas = temp_imeas
              
              deallocate(temp_misfits, temp_residuals, temp_tstart, temp_tend, temp_imeas)
            enddo
            deallocate(recv_nwin_arr)
          enddo

          deallocate(recv_indices, recv_ncomp, recv_net, recv_sta, recv_total_misfit)
        end if ! nsta_irank > 0
      enddo ! iproc
      
    else ! worldrank != 0
      if (nrec_loc > 0) then
        allocate(send_indices(nrec_loc))
        allocate(recv_ncomp(nrec_loc))
        allocate(recv_net(nrec_loc))
        allocate(recv_sta(nrec_loc))
        allocate(recv_total_misfit(nrec_loc))
        
        ! Prepare basic data for each record
        do irec_loc = 1, nrec_loc
          send_indices(irec_loc) = select_global_id_for_rec(irec_loc)
          recv_ncomp(irec_loc) = misfit_loc(irec_loc)%ncomp
          recv_net(irec_loc) = misfit_loc(irec_loc)%net
          recv_sta(irec_loc) = misfit_loc(irec_loc)%sta
          recv_total_misfit(irec_loc) = misfit_loc(irec_loc)%total_misfit
        enddo
        
        ! Send basic information
        call send_i(send_indices, nrec_loc, 0, targ)
        call send_i(recv_ncomp, nrec_loc, 0, targ)
        call send_ch_array(recv_net, nrec_loc, MAX_STRING_LEN, 0, targ)
        call send_ch_array(recv_sta, nrec_loc, MAX_STRING_LEN, 0, targ)
        call send_dp(recv_total_misfit, nrec_loc, 0, targ)
        
        ! Send detailed data for each record
        do irec_loc = 1, nrec_loc
          ncomp = misfit_loc(irec_loc)%ncomp
          
          ! Send channel names
          allocate(recv_chan(1, ncomp))
          recv_chan(1, :) = misfit_loc(irec_loc)%chan
          call send_ch_array(recv_chan, ncomp, MAX_STRING_LEN, 0, targ)
          deallocate(recv_chan)
          
          ! Send nwin for each component
          allocate(recv_nwin_arr(ncomp))
          do icomp = 1, ncomp
            recv_nwin_arr(icomp) = misfit_loc(irec_loc)%trm(icomp)%nwin
          enddo
          call send_i(recv_nwin_arr, ncomp, 0, targ)
          
          ! Send TraceMisfit data for each component
          do icomp = 1, ncomp
            nwin = recv_nwin_arr(icomp)
            call send_dp(misfit_loc(irec_loc)%trm(icomp)%misfits, nwin, 0, targ)
            call send_dp(misfit_loc(irec_loc)%trm(icomp)%residuals, nwin, 0, targ)
            call send_dp(misfit_loc(irec_loc)%trm(icomp)%tstart, nwin, 0, targ)
            call send_dp(misfit_loc(irec_loc)%trm(icomp)%tend, nwin, 0, targ)
            call send_i(misfit_loc(irec_loc)%trm(icomp)%imeas, nwin, 0, targ)
          enddo
          deallocate(recv_nwin_arr)
        enddo
        
        deallocate(send_indices, recv_ncomp, recv_net, recv_sta, recv_total_misfit)
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
    integer :: irec, icomp, iwin
    character(len=MAX_STRING_LEN) :: chifile

    if (worldrank /= 0) return

    ! Write window chi file
    chifile = trim(MISFITS_DIR)//'/'//trim(model_name)//'.'//&
              trim(this%evtname)//'_'//trim(this%band_name)//&
              '_window_chi'
    open(newunit=unit, file=chifile, status='unknown', iostat=ier)
    if (ier /= 0) then
      write(*,*) 'Error opening file: ', trim(chifile)
      error stop
    end if

    do irec = 1, this%nrec
      do icomp = 1, this%rec_misfits(irec)%ncomp
        do iwin = 1, this%rec_misfits(irec)%trm(icomp)%nwin
          write(unit, '(a,1x,a,1x,a,1x,i4,1x,f8.3,1x,f8.3,1x,f12.6,1x,f12.6)') &
                       trim(this%rec_misfits(irec)%net), &
                       trim(this%rec_misfits(irec)%sta), &
                       trim(this%rec_misfits(irec)%chan(icomp)), &
                       this%rec_misfits(irec)%trm(icomp)%imeas(iwin), &
                       this%rec_misfits(irec)%trm(icomp)%tstart(iwin), &
                       this%rec_misfits(irec)%trm(icomp)%tend(iwin), &
                       this%rec_misfits(irec)%trm(icomp)%residuals(iwin), &
                       this%rec_misfits(irec)%trm(icomp)%misfits(iwin)
        enddo ! iwin
      enddo ! icomp
    enddo ! irec
    close(unit)

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
    else
      nflt = fpar%sim%NUM_FILTER
    end if

    total_misfit = 0.0_dp
    if (worldrank == 0) then
      do iflt = 1, nflt
        if (dat_type == 'rf') then
          write(band_name, '("F",F3.1)') fpar%sim%rf%f0(iflt)
        else
          call get_band_name(fpar%sim%SHORT_P(iflt), fpar%sim%LONG_P(iflt), band_name)
        end if
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