module syn_data
  use config
  use fwat_constants, only: targ
  use specfem_par, only: NSTEP, DT, T0, ELASTIC_SIMULATION, seismograms_d, nrec, nrec_local, &
                        number_receiver_global, ispec_selected_rec, islice_selected_rec
  use specfem_par_elastic, only: ispec_is_elastic
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use fwat_mpi
  use input_params, fpar => fwat_par_global
  use common_lib, only: rotate_NE_to_RT, rotate_NE_to_RT_dp, dwascii
  use obs_data, only: ObsData
  use window_chi, only: WindowChi
  use utils, only: zeros_dp

  implicit none

  integer, private :: ier
  integer, parameter, private :: NCOMP=3

  type :: SynData
    character(len=MAX_STRING_LEN), dimension(3) :: comp_name
    integer :: ievt, nrec
    type(ObsData) :: od
    type(WindowChi), dimension(:), allocatable :: wchi
    real(kind=dp), dimension(:, :, :), pointer :: data ! npts, ncomp(zrt), nsta
    character(len=MAX_STRING_LEN) :: band_name
    real(kind=dp), dimension(:), allocatable :: total_misfit
    integer :: dat_win
    contains
    procedure :: read=>read_syn_data, filter, assemble_2d, assemble_3d, init
    procedure :: get_comp_name_adj, finalize, write_adj

  end type SynData

contains

  subroutine init(this, ievt)
    class(SynData), intent(inout) :: this
    integer, intent(in) :: ievt

    this%ievt = ievt
    this%nrec = nrec

  end subroutine init

  subroutine read_syn_data(this, bazi)
    class(SynData), intent(inout) :: this
    real(kind=cr), dimension(:), intent(in) :: bazi
    integer :: irec, irec_local, ispec, iproc, nsta_irank, i
    real(kind=dp), dimension(:, :, :), allocatable :: data_local, dat_sum
    real(kind=dp), dimension(:,:,:), allocatable :: recv_buffer
    integer, dimension(:), allocatable :: recv_indices, send_indices

    ! read source and reciever files
    call prepare_shm_array_dp_3d(this%data, NSTEP, NCOMP, nrec, this%dat_win)
    ! allocate(data_local(NSTEP, NCOMP, nrec_local))
    ! data_local = 0.0_dp

    if (nrec_local > 0) then
      data_local = zeros_dp(NSTEP, NCOMP, nrec_local)
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        ispec = ispec_selected_rec(irec)
        if (ELASTIC_SIMULATION) then
          if (ispec_is_elastic(ispec)) then
            call rotate_NE_to_RT_dp(dble(seismograms_d(2, irec_local, 1:NSTEP)), dble(seismograms_d(1, irec_local, 1:NSTEP)),&
                                    data_local(:, 2, irec_local), data_local(:, 3, irec_local), bazi(irec))
            data_local(:, 1, irec_local) = dble(seismograms_d(3, irec_local, 1:NSTEP))
            ! data_local(:, 2, irec_local) = dble(seismograms_d(2, irec_local, 1:NSTEP))
            ! data_local(:, 3, irec_local) = dble(seismograms_d(1, irec_local, 1:NSTEP))

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
          allocate(recv_buffer(NSTEP, NCOMP, nsta_irank))  ! Allocate a buffer to receive data
          allocate(recv_indices(nsta_irank)) 
          ! Receive the indices first
          call recv_i(recv_indices, nsta_irank, iproc, targ)
          ! Receive the data
          call recv_dp(recv_buffer, NSTEP*NCOMP*nsta_irank, iproc, targ)
          ! Copy the received data to the correct location
          do i = 1, nsta_irank
            irec = recv_indices(i)
            this%data(:, :, irec) = recv_buffer(:, :, i)
          enddo
          deallocate(recv_buffer)
          deallocate(recv_indices)
        endif
      enddo
    else
      if (nrec_local > 0) then
        allocate(send_indices(nrec_local))
        do irec_local = 1, nrec_local
          send_indices(irec_local) = number_receiver_global(irec_local)
        enddo
        call send_i(send_indices, nrec_local, 0, targ)
        call send_dp(data_local(:, :, :), NSTEP*NCOMP*nrec_local, 0, targ)
        deallocate(send_indices)
      endif
    endif
    call sync_from_main_rank_dp_3d(this%data, NSTEP, NCOMP, nrec)

  end subroutine read_syn_data

  subroutine assemble_2d(this, array_local, array_global)
    class(SynData), intent(inout) :: this
    real(kind=dp), dimension(:,:), intent(in) :: array_local
    real(kind=dp), dimension(:,:), allocatable, intent(out) :: array_global
    integer :: irec, irec_local, nsta_irank, iproc, i
    real(kind=dp), dimension(:,:), allocatable :: recv_buffer
    integer, dimension(:), allocatable :: recv_indices, send_indices

    if (worldrank == 0) then
      array_global = zeros_dp(NSTEP, nrec)
      if (nrec_local > 0) then
        do irec_local = 1, nrec_local
          irec = number_receiver_global(irec_local)
          array_global(:, irec) = array_local(:, irec_local)
        enddo
      endif
      do iproc = 1, worldsize-1
        nsta_irank=0
        do irec = 1, nrec
          if (islice_selected_rec(irec) == iproc) nsta_irank = nsta_irank + 1
        enddo
      
        if (nsta_irank > 0) then
          allocate(recv_buffer(NSTEP, nsta_irank))  ! Allocate a buffer to receive data
          allocate(recv_indices(nsta_irank)) 
          ! Receive the indices first
          call recv_i(recv_indices, nsta_irank, iproc, targ)
          ! Receive the data
          call recv_dp(recv_buffer, NSTEP*nsta_irank, iproc, targ)
          ! Copy the received data to the correct location
          do i = 1, nsta_irank
            irec = recv_indices(i)
            array_global(:, irec) = recv_buffer(:, i)
          enddo
          deallocate(recv_buffer)
          deallocate(recv_indices)
        endif
      enddo
    else
      if (nrec_local > 0) then
        allocate(send_indices(nrec_local))
        do irec_local = 1, nrec_local
          send_indices(irec_local) = number_receiver_global(irec_local)
        enddo
        call send_i(send_indices, nrec_local, 0, targ)
        call send_dp(array_local(:, :), NSTEP*nrec_local, 0, targ)
        deallocate(send_indices)
      endif
    endif

  end subroutine assemble_2d

  subroutine assemble_3d(this, array_local, array_global, nc)
    ! Assemble 3D array from local 3D array
    ! Input:
    !   array_local: local 3D array
    !   nc: number of components
    ! Output:
    !   array_global: global 3D array
    class(SynData), intent(inout) :: this
    real(kind=dp), dimension(:,:,:), intent(in) :: array_local
    real(kind=dp), dimension(:,:,:), allocatable, intent(out) :: array_global
    integer, intent(in) :: nc
    integer :: irec, irec_local, nsta_irank, iproc, i
    real(kind=dp), dimension(:,:,:), allocatable :: recv_buffer
    integer, dimension(:), allocatable :: recv_indices, send_indices

    ! if (.not. allocated(array_global)) allocate(array_global(NSTEP, nc, nrec))
    ! array_global = 0.0_dp

    if (worldrank == 0) then
      array_global = zeros_dp(NSTEP, nc, nrec)
      if (nrec_local > 0) then
        do irec_local = 1, nrec_local
          irec = number_receiver_global(irec_local)
          array_global(:, :, irec) = array_local(:, :, irec_local)
        enddo
      endif
      do iproc = 1, worldsize-1
        nsta_irank=0
        do irec = 1, nrec
          if (islice_selected_rec(irec) == iproc) nsta_irank = nsta_irank + 1
        enddo
        if (nsta_irank > 0) then
          allocate(recv_buffer(NSTEP, nc, nsta_irank))  ! Allocate a buffer to receive data
          allocate(recv_indices(nsta_irank)) 
          ! Receive the indices first
          call recv_i(recv_indices, nsta_irank, iproc, targ)
          ! Receive the data
          call recv_dp(recv_buffer, NSTEP*nc*nsta_irank, iproc, targ)
          ! Copy the received data to the correct location
          do i = 1, nsta_irank
            irec = recv_indices(i)
            array_global(:, :, irec) = recv_buffer(:, :, i)
          enddo
          deallocate(recv_buffer)
          deallocate(recv_indices)
        endif
      enddo
    else
      if (nrec_local > 0) then
        allocate(send_indices(nrec_local))
        do irec_local = 1, nrec_local
          send_indices(irec_local) = number_receiver_global(irec_local)
        enddo
        call send_i(send_indices, nrec_local, 0, targ)
        call send_dp(array_local(:, :, :), NSTEP*nc*nrec_local, 0, targ)
        deallocate(send_indices)
      endif
    endif
  end subroutine assemble_3d

  subroutine write_adj(this, adj_data, kcmp, irec)
    class(SynData), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: adj_data
    character(len=*), intent(in) :: kcmp
    integer, intent(in) :: irec
    character(len=MAX_STRING_LEN) :: adj_file

    adj_file = trim(fpar%acqui%out_fwd_path(this%ievt))//'/'//trim(ADJOINT_PATH)//&
               '/'//trim(this%od%netwk(irec))//'.'//trim(this%od%stnm(irec))//&
               '.'//trim(fpar%sim%CH_CODE)//trim(kcmp)//'.adj'
    call dwascii(adj_file, adj_data, NSTEP, -dble(T0), dble(DT))
  end subroutine write_adj

  ! subroutine rotate_to_rt(this, baz)
  !   class(SynData), intent(inout) :: this
  !   real(kind=cr), dimension(:), intent(in) :: baz
  !   real(kind=cr), dimension(:), allocatable :: data_r, data_t
  !   integer :: irec, it

  !   allocate(data_r(NSTEP), data_t(NSTEP))
  !   if (noderank == 0) then
  !     do irec = 1, nrec
  !       call rotate_NE_to_RT(real(this%data(:, 3, irec)), real(this%data(:, 2, irec)), data_r, data_t, baz(irec))
  !       this%data(:, 2, irec) = dble(data_r)
  !       this%data(:, 3, irec) = dble(data_t)
  !     enddo
  !   endif
  !   call synchronize_all()

  ! end subroutine rotate_to_rt

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

  function average_amp_scale(glob_dat_tw, icomp) result(avgamp)
    real(kind=cr) :: avgamp
    real(kind=cr) :: avgamp0
    integer :: igood, icomp, irec
    real(kind=dp), dimension(:,:,:)   :: glob_dat_tw

    ! use only Z component for amplitude scale
    avgamp0=0.
    do irec =1 ,nrec
      avgamp0=avgamp0+maxval(abs(glob_dat_tw(:,icomp,irec))) 
    enddo
    avgamp0=avgamp0/nrec
    avgamp=0
    igood=0
    do irec =1, nrec
      if ((maxval(abs(glob_dat_tw(:,icomp,irec)))-avgamp0)<0.2*avgamp0) then
        avgamp=avgamp+maxval(abs(glob_dat_tw(:,icomp,irec)))
        igood=igood+1
      endif
    enddo
    avgamp=avgamp/igood
  end function average_amp_scale

  subroutine get_comp_name_adj(this)
    class(SynData), intent(inout) :: this
    
    if (SUPPRESS_UTM_PROJECTION) then
      this%comp_name = ['Z', 'Y', 'X']
    else
      this%comp_name = ['Z', 'N', 'E']
    endif
    call synchronize_all()
  end subroutine get_comp_name_adj

  subroutine finalize(this)
    class(SynData), intent(inout) :: this
    integer :: i

    call this%od%finalize()
    do i = 1, size(this%wchi)
      call this%wchi(i)%finalize()
    enddo
    call free_shm_array(this%dat_win)

  end subroutine finalize

end module syn_data