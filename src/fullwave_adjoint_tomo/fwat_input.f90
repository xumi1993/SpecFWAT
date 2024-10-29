!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!                  
!=====================================================================
!
module fwat_input
  use constants, only :MAX_STRING_LEN 
  use fullwave_adjoint_tomo_par
  use ma_variables
  use meshgrid, only: ReglGrid

  type(fwat_noise_parameters), save, public :: noise_par
  type(fwat_tele_parameters), save, public :: tele_par
  type(fwat_rf_parameters), save, public :: rf_par
  type(fwat_leq_parameters), save, public :: leq_par
  type(fwat_tomo_parameters), save, public :: tomo_par
  type(fwat_acqui), save, public :: acqui_par
  type(ReglGrid), save, public :: rg

contains

  subroutine read_fwat_par_file()
    implicit none
    
    ! read common parameters
    call fwat_common_read_par_file()

    ! read noise para
    call noise_par%read_noise_par_file()

    ! read tele para
    call tele_par%read_tele_par_file()
    
    ! read rf para
    call rf_par%read_rf_par_file()

    call leq_par%read_leq_par_file()

    call tomo_par%read_tomo_par_file()

    call fwat_read_meas_adj_par_file()

  end subroutine read_fwat_par_file

  subroutine fwat_parameters_free()
    implicit none
    
    ! free common variables
    if (allocated(RCOMPS)) deallocate(RCOMPS)
    if (allocated(RCOMPS)) deallocate(SHORT_P)
    if (allocated(RCOMPS)) deallocate(LONG_P)
    call acqui_par%free()
  
    ! ! free noise variables
    ! if (allocated(noise_par%SCOMPS)) deallocate(noise_par%SCOMPS)
    ! if (allocated(noise_par%RCOMPS)) deallocate(noise_par%RCOMPS)
    ! if (allocated(noise_par%SHORT_P)) deallocate(noise_par%SHORT_P)
    ! if (allocated(noise_par%LONG_P)) deallocate(noise_par%LONG_P)
    ! if (allocated(noise_par%GROUPVEL_MIN)) deallocate(noise_par%GROUPVEL_MIN)
    ! if (allocated(noise_par%GROUPVEL_MAX)) deallocate(noise_par%GROUPVEL_MAX)

    ! ! free tele variables
    ! if (allocated(tele_par%RCOMPS)) deallocate(tele_par%RCOMPS)
    ! if (allocated(tele_par%SHORT_P)) deallocate(tele_par%SHORT_P)
    ! if (allocated(tele_par%LONG_P)) deallocate(tele_par%LONG_P)

    ! ! free rf variables
    ! if (allocated(rf_par%F0)) deallocate(rf_par%F0)

    ! ! free leq variables
    ! if (allocated(leq_par%SHORT_P)) deallocate(leq_par%SHORT_P)
    ! if (allocated(leq_par%LONG_P)) deallocate(leq_par%LONG_P)

    ! ! free tomo variables
    ! if (allocated(tomo_par%STEP_LENS)) deallocate(tomo_par%STEP_LENS)
    ! if (allocated(tomo_par%NOISE_SET_NAMES)) deallocate(tomo_par%NOISE_SET_NAMES)
    ! if (allocated(tomo_par%TELE_SET_NAMES)) deallocate(tomo_par%TELE_SET_NAMES)

  end subroutine fwat_parameters_free

  subroutine fwat_read_meas_adj_par_file()
    use ma_constants

    implicit none
    integer            :: myrank,ipos0,ipos1
    character(len=MAX_STRING_LEN)  :: line,keyw  
    ! get current proc
    call world_rank(myrank)

    DISPLAY_DETAILS = .false.
    OUTPUT_MEASUREMENT_FILES = .false.
    RUN_BANDPASS = .false.
    COMPUTE_ADJOINT_SOURCE = .true.
    USE_PHYSICAL_DISPERSION = .false.
    ! begin reading params 
    if(myrank == 0) then 
      open(666,file=FWAT_PAR_FILE) 
      do 
        read(666,'(a)',end=99) line
        if (is_blank_line(line) .or. line(1:1) == '#') cycle
        !! INDICES TO READ line -----------------------------------------------
        ipos0=index(line,':')+1
        ipos1=index(line,'#')-1
        if (ipos1 < 0 ) ipos1=len_trim(line)

        !! STORE KEYWORD ITEM -------------------------------------------------
        keyw=trim(adjustl(line(1:ipos0-2)))
        !! DIFFERENT ITEM TO READ ---------------------------------------------
        select case (trim(keyw))
          case('TSHIFT_MIN')
            read(line(ipos0:ipos1),*) TSHIFT_MIN
          case('TSHIFT_MAX')
            read(line(ipos0:ipos1),*) TSHIFT_MAX
          case('DLNA_MIN')
            read(line(ipos0:ipos1),*) DLNA_MIN
          case('DLNA_MAX')
            read(line(ipos0:ipos1),*) DLNA_MAX
          case('CC_MIN')
            read(line(ipos0:ipos1),*) CC_MIN
          case('ERROR_TYPE')
            read(line(ipos0:ipos1),*) ERROR_TYPE
          case('DT_SIGMA_MIN')
            read(line(ipos0:ipos1),*) DT_SIGMA_MIN
          case('DLNA_SIGMA_MIN')
            read(line(ipos0:ipos1),*) DLNA_SIGMA_MIN
          case('WTR')
            read(line(ipos0:ipos1),*) WTR
          case('NPI')
            read(line(ipos0:ipos1),*) NPI
          case('DT_FAC')
            read(line(ipos0:ipos1),*) DT_FAC
          case('ERR_FAC')
            read(line(ipos0:ipos1),*) ERR_FAC
          case('DT_MAX_SCALE')
            read(line(ipos0:ipos1),*) DT_MAX_SCALE
          case('NCYCLE_IN_WINDOW')
            read(line(ipos0:ipos1),*) NCYCLE_IN_WINDOW
          case('USE_PHYSICAL_DISPERSION')
            read(line(ipos0:ipos1),*) USE_PHYSICAL_DISPERSION
        end select
      enddo
      if( DO_RAY_DENSITY_SOURCE ) ERROR_TYPE = 0
    endif
99 close(666) ! close par file
    call bcast_all_singledp(TSHIFT_MIN)
    call bcast_all_singledp(TSHIFT_MAX)
    call bcast_all_singledp(DLNA_MIN)
    call bcast_all_singledp(DLNA_MAX)
    call bcast_all_singledp(CC_MIN)
    call bcast_all_singledp(DT_SIGMA_MIN)
    call bcast_all_singledp(DLNA_SIGMA_MIN)
    call bcast_all_singlei(ERROR_TYPE)
    call bcast_all_singledp(WTR)
    call bcast_all_singledp(NPI)
    call bcast_all_singledp(DT_FAC)
    call bcast_all_singledp(ERR_FAC)
    call bcast_all_singledp(DT_MAX_SCALE)
    call bcast_all_singledp(NCYCLE_IN_WINDOW)
    call bcast_all_singlel(USE_PHYSICAL_DISPERSION)
    
  end subroutine fwat_read_meas_adj_par_file

  subroutine setup_common_variables(simu_type)
    use specfem_par, only: NSTEP, DT
    character(len=MAX_STRING_LEN), intent(in) :: simu_type

    ! if(.not. allocated(RCOMPS)) allocate(RCOMPS(NRCOMP))
    ! if(.not. allocated(SHORT_P)) allocate(SHORT_P(NUM_FILTER))
    ! if(.not. allocated(LONG_P)) allocate(LONG_P(NUM_FILTER))
    ! SHORT_P(:) = 0.
    ! LONG_P(:) = 0.
    if (allocated(RCOMPS)) deallocate(RCOMPS)
    if (allocated(SHORT_P)) deallocate(SHORT_P)
    if (allocated(LONG_P)) deallocate(LONG_P)

    if (simu_type == 'rf' .or. index(simu_type, 'tele') /= 0) then
      NRCOMP = tele_par%NRCOMP
      allocate(RCOMPS(NRCOMP))
      RCOMPS = tele_par%RCOMPS
      CH_CODE = tele_par%CH_CODE
      dat_coord = tele_par%dat_coord
      NUM_FILTER = tele_par%NUM_FILTER
      allocate(SHORT_P(NUM_FILTER))
      allocate(LONG_P(NUM_FILTER))
      SHORT_P = tele_par%SHORT_P
      LONG_P = tele_par%LONG_P
      imeas0 = tele_par%IMEAS
      imeas = imeas0
      is_mtm0 = 0
      is_mtm = is_mtm0
      DT = tele_par%DT
      NSTEP = tele_par%NSTEP
      PRECOND_TYPE = tele_par%PRECOND_TYPE
    elseif (simu_type == 'noise') then
      NRCOMP = noise_par%NRCOMP
      allocate(RCOMPS(NRCOMP))
      RCOMPS = noise_par%RCOMPS
      CH_CODE = noise_par%CH_CODE
      dat_coord = noise_par%dat_coord
      NUM_FILTER = noise_par%NUM_FILTER
      allocate(SHORT_P(NUM_FILTER))
      allocate(LONG_P(NUM_FILTER))
      SHORT_P = noise_par%SHORT_P
      LONG_P = noise_par%LONG_P
      imeas0 = noise_par%IMEAS
      imeas = imeas0
      is_mtm0 = 1
      is_mtm = is_mtm0
      DT = noise_par%DT
      NSTEP = noise_par%NSTEP
      PRECOND_TYPE = noise_par%PRECOND_TYPE
    elseif (simu_type == 'leq') then
      NRCOMP = leq_par%NRCOMP
      allocate(RCOMPS(NRCOMP))
      RCOMPS = leq_par%RCOMPS
      CH_CODE = leq_par%CH_CODE
      dat_coord = leq_par%dat_coord
      NUM_FILTER = leq_par%NUM_FILTER
      allocate(SHORT_P(NUM_FILTER))
      allocate(LONG_P(NUM_FILTER))
      SHORT_P = leq_par%SHORT_P
      LONG_P = leq_par%LONG_P
      imeas0 = leq_par%IMEAS
      imeas = imeas0
      if ( imeas == 1 .or. imeas == 2 ) then ! waveforms
        is_mtm0 = 0
      else if ( imeas >= 3 .and. imeas <= 6 ) then ! CC Dt/DlnA
        ! for CC kernels, ITAPER must be a single taper (2 or 3)
        if ( leq_par%ITAPER == 1 ) stop  'Error: Change ITAPER to 2/3 for CC measurements'
        is_mtm0 = leq_par%ITAPER     ! 2 or 3 for CC adjoint sources
      else if ( imeas==7 .or. imeas==8 ) then
        is_mtm0 = 1          ! multitaper required for MT adjoint source
      else
        stop 'Error: imeas must by 1-8'
      endif
      is_mtm = is_mtm0
      DT = leq_par%DT
      NSTEP = leq_par%NSTEP
      PRECOND_TYPE = leq_par%PRECOND_TYPE
    endif
    call synchronize_all()
  end subroutine setup_common_variables

  subroutine fwat1_out_log(model, evtset, simu_type)
    implicit none

    character(len=MAX_STRING_LEN) :: model,evtset,simu_type
    integer :: myrank
    logical :: io_stat

    call world_rank(myrank)

    if (myrank == 0) then
      inquire(unit=OUT_FWAT_LOG, opened=io_stat)
      if (.not. io_stat) then
        open(unit=OUT_FWAT_LOG,file='output_fwat1_log_'//trim(model)//'.'//trim(evtset)//'.txt')
      endif
      ! write(OUT_FWAT_LOG,*) 'This is run_fwat1_fwd_measure_adj !!!'
      write(OUT_FWAT_LOG,*) 'model,evtset,nevents: ',trim(model),' ', trim(evtset),' ',acqui_par%nevents
      if (simu_type == 'noise') then
        write(OUT_FWAT_LOG,*) '============= Noise parameters =============='
        write(OUT_FWAT_LOG,*) 'NOISE_SHORT_P=',noise_par%SHORT_P(:)
        write(OUT_FWAT_LOG,*) 'NOISE_LONG_P=',noise_par%LONG_P(:)
        write(OUT_FWAT_LOG,*) 'NOISE_GROUPVEL_MIN=',noise_par%GROUPVEL_MIN(:)
        write(OUT_FWAT_LOG,*) 'NOISE_GROUPVEL_MAX=',noise_par%GROUPVEL_MAX(:)
        write(OUT_FWAT_LOG,*) 'NOISE_ADJ_SRC_NORM: ',noise_par%ADJ_SRC_NORM
        write(OUT_FWAT_LOG,*) 'NOISE_USE_NEAR_OFFSET: ',noise_par%USE_NEAR_OFFSET
        write(OUT_FWAT_LOG,*) 'NOISE_SUPPRESS_EGF: ',noise_par%SUPPRESS_EGF
        write(OUT_FWAT_LOG,*) 'NOISE_IMEAS=',noise_par%IMEAS
        write(OUT_FWAT_LOG,*) 'NOISE_ITAPER=',noise_par%ITAPER
      endif
      if (simu_type == 'rf' .or. index(simu_type, 'tele') /= 0) then
        write(OUT_FWAT_LOG,*) '============= Tele parameters ==============='
        write(OUT_FWAT_LOG,*) 'TELE_SHORT_P=',tele_par%SHORT_P(:)
        write(OUT_FWAT_LOG,*) 'TELE_LONG_P=',tele_par%LONG_P(:)
        write(OUT_FWAT_LOG,*) 'TELE_TW_BEFORE,TELE_TW_AFTER=',tele_par%TW_BEFORE,tele_par%TW_AFTER
        write(OUT_FWAT_LOG,*) 'TELE_IMEAS=',tele_par%IMEAS
        write(OUT_FWAT_LOG,*) 'TELE_ITAPER=',tele_par%ITAPER
      endif
      if (simu_type == 'rf') then
        write(OUT_FWAT_LOG,*) '============== RF parameters ================'
        write(OUT_FWAT_LOG,*) 'RF_F0=',rf_par%F0(:)
        write(OUT_FWAT_LOG,*) 'RF_TSHIFT=',rf_par%RF_TSHIFT
        write(OUT_FWAT_LOG,*) 'RF_MAXIT=',rf_par%MAXIT
        write(OUT_FWAT_LOG,*) 'RF_MINDERR=',rf_par%MINDERR
      endif
      if (simu_type == 'leq') then
        write(OUT_FWAT_LOG,*) '=========== Local EQ parameters =============='
        write(OUT_FWAT_LOG,*) 'LEQ_SHORT_P=',leq_par%SHORT_P(:)
        write(OUT_FWAT_LOG,*) 'LEQ_LONG_P=',leq_par%LONG_P(:)
      endif
      write(OUT_FWAT_LOG,*) '*******************************************************'
      flush(OUT_FWAT_LOG)
    endif

  end subroutine fwat1_out_log

end module fwat_input
