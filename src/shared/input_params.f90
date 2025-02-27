module input_params
  use config
  use fwat_mpi

  implicit none

  integer, private :: ier

  type acqui_params
    integer               :: nevents ! # of events
    character(len= MAX_STRING_LEN) :: evtset_file
    character(len= MAX_STRING_LEN), dimension(:),       pointer  :: station_file
    character(len= MAX_STRING_LEN), dimension(:),       pointer  :: src_solution_file
    character(len= MAX_STRING_LEN), dimension(:),       pointer  :: evtid_names
    character(len= MAX_STRING_LEN), dimension(:),       pointer  :: out_fwd_path
    character(len= MAX_STRING_LEN), dimension(:),       pointer  :: in_dat_path
    character(len= MAX_STRING_LEN), dimension(:),       pointer  :: fkmodel_file
    real(kind=cr), dimension(:), pointer :: evla,evlo,evdp,src_weight
    integer :: sta_win, src_win, evid_win, ofp_win, idp_win, evla_win, evlo_win,&
               evdp_win, w_win, fk_win

    contains
    procedure, private :: alloc => acqui_alloc
    procedure :: read => acqui_read_source_set
  end type acqui_params

  type rf_params
    real(kind=cr)                                    :: MINDERR, TSHIFT
    integer                                          :: NGAUSS, MAXIT
    real(kind=cr), dimension(:), allocatable :: F0
  end type rf_params

  type sim_params
    integer :: NRCOMP, NSCOMP, NUM_FILTER, NSTEP, IMEAS, ITAPER, PRECOND_TYPE
    character(len= MAX_STRING_LEN), dimension(:), allocatable :: RCOMPS, SCOMPS
    character(len= MAX_STRING_LEN) :: CH_CODE
    real(kind=cr) :: DT
    real(kind=cr), dimension(:), allocatable :: SHORT_P, LONG_P, GROUPVEL_MIN, GROUPVEL_MAX, TIME_WIN
    logical :: USE_NEAR_OFFSET, ADJ_SRC_NORM, SUPPRESS_EGF, USE_LOCAL_STF
    type(rf_params) :: rf
  end type sim_params

  type fwat_params
    type(acqui_params) :: acqui
    type(sim_params), pointer :: sim
    contains
    procedure :: read => read_parameter_file, select_simu_type
  end type fwat_params

  type(sim_params), target :: tele_par, noise_par
  type(fwat_params) :: fwat_par_global

contains
  subroutine acqui_read_source_set(this)
    implicit none

    class(acqui_params), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: evtnm, line
    real(kind=cr) :: junk_cr
    integer, parameter :: FID = 878
    integer :: ievt

    this%evtset_file = trim(SRC_REC_DIR)//'/'//trim(SRC_PREFIX)//'_'//trim(dat_type)//'.dat'
    if (worldrank == 0) then
      open(unit=FID, file=this%evtset_file, status='old', iostat=ier)
      if (ier /= 0) call exit_mpi(worldrank, 'ERROR: cannot open file '//trim(this%evtset_file))
      ievt = 0
      do 
        read(FID, *, iostat=ier) evtnm
        if (ier /= 0) exit
        ievt = ievt + 1
      end do
      close(FID)
      this%nevents = ievt
    end if ! worldrank == 0
    call bcast_all_singlei(this%nevents)

    ! allocate arrays
    call this%alloc()

    if (worldrank == 0) then
      open(unit=FID, file=this%evtset_file, status='old', iostat=ier)
      if (ier /= 0) call exit_mpi(worldrank, 'ERROR: cannot open file '//trim(this%evtset_file))
      do ievt = 1, this%nevents
        read(FID, *, iostat=ier) line
        read(line,*,iostat=ier) evtnm,this%evla(ievt),this%evlo(ievt), &
                                this%evdp(ievt),junk_cr,&
                                this%src_weight(ievt)
        if (ier /= 0) this%src_weight(ievt) = 1.0_cr
        this%station_file(ievt) = trim(SRC_REC_DIR)//'/'//trim(STATIONS_PREFIX)//'_'//trim(evtnm)
        if (simu_type == SIMU_TYPE_NOISE) then
          this%src_solution_file(ievt) = trim(SRC_REC_DIR)//'/'//trim(FORCESOLUTION_PREFIX)//'_'//trim(evtnm)
        else if (simu_type == SIMU_TYPE_TELE) then
          ! this%src_solution_file(ievt) = trim(SRC_REC_DIR)//'/'//trim(CMTSOLUTION_PREFIX)//'_'//trim(evtnm)
          this%src_solution_file(ievt) = 'DATA/'//trim(CMTSOLUTION_PREFIX)
          this%fkmodel_file(ievt) = trim(SRC_REC_DIR)//'/'//trim(FKMODEL_PREFIX)//'_'//trim(evtnm)
        endif
        this%evtid_names(ievt) = evtnm
        this%out_fwd_path(ievt)=trim(SOLVER_DIR)//'/'//trim(model_name)//'.'//trim(dat_type)//'/'//trim(evtnm)
        this%in_dat_path(ievt) = trim(DATA_DIR)//'/'//trim(evtnm)//'/'
      end do
      close(FID)
    end if ! worldrank == 0
    call synchronize_all()

    ! broadcast arrays
    call sync_from_main_rank_ch(this%station_file, this%nevents, MAX_STRING_LEN)
    call sync_from_main_rank_ch(this%src_solution_file, this%nevents, MAX_STRING_LEN)
    call sync_from_main_rank_ch(this%evtid_names, this%nevents, MAX_STRING_LEN)
    call sync_from_main_rank_ch(this%out_fwd_path, this%nevents, MAX_STRING_LEN)
    call sync_from_main_rank_ch(this%in_dat_path, this%nevents, MAX_STRING_LEN)
    if (simu_type == SIMU_TYPE_TELE) then
      call sync_from_main_rank_ch(this%fkmodel_file, this%nevents, MAX_STRING_LEN)
    endif
    call sync_from_main_rank_cr_1d(this%evla, this%nevents)
    call sync_from_main_rank_cr_1d(this%evlo, this%nevents)
    call sync_from_main_rank_cr_1d(this%evdp, this%nevents)
    call sync_from_main_rank_cr_1d(this%src_weight, this%nevents)

  end subroutine acqui_read_source_set

  subroutine acqui_alloc(this)
    class(acqui_params), intent(inout) :: this

    ! allocate arrays on shared memory
    call prepare_shm_array_ch_1d(this%station_file, this%nevents, MAX_STRING_LEN, this%sta_win)
    call prepare_shm_array_ch_1d(this%src_solution_file, this%nevents, MAX_STRING_LEN, this%src_win)
    call prepare_shm_array_ch_1d(this%evtid_names, this%nevents, MAX_STRING_LEN, this%evid_win)
    call prepare_shm_array_ch_1d(this%out_fwd_path, this%nevents, MAX_STRING_LEN, this%ofp_win)
    call prepare_shm_array_ch_1d(this%in_dat_path, this%nevents, MAX_STRING_LEN, this%idp_win)
    if (simu_type == SIMU_TYPE_TELE) then
      call prepare_shm_array_ch_1d(this%fkmodel_file, this%nevents, MAX_STRING_LEN, this%fk_win)
    endif
    call prepare_shm_array_cr_1d(this%evla, this%nevents, this%evla_win)
    call prepare_shm_array_cr_1d(this%evlo, this%nevents, this%evlo_win)
    call prepare_shm_array_cr_1d(this%evdp, this%nevents, this%evdp_win)
    call prepare_shm_array_cr_1d(this%src_weight, this%nevents, this%w_win)

  end subroutine acqui_alloc

  subroutine read_parameter_file(this, fname)
    use yaml, only: parse, error_length
    use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                        type_list, type_list_item, type_scalar
    class(fwat_params), intent(inout) :: this
    character(len=*), intent(in) :: fname
    class(type_node), pointer :: root
    class(type_dictionary), pointer :: noise, tele, tomo, rf
    class (type_list), pointer :: list
    character(len=error_length) :: error
    type (type_error), pointer :: io_err

    if (worldrank == 0) then
      root => parse(fname, error = error)
      if (error/='') call exit_mpi(worldrank, error)
      select type (root)
      class is (type_dictionary)
        ! read parameters for noise FWI
        this%sim => noise_par
        this%sim%IMEAS = 7
        this%sim%ITAPER = 1
        noise => root%get_dictionary('NOISE', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        this%sim%NSTEP = noise%get_integer('NSTEP', error=io_err)
        this%sim%PRECOND_TYPE = noise%get_integer('PRECOND_TYPE', error=io_err)
        list => noise%get_list('RCOMPS', required=.true., error=io_err)
        if(associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        call read_string_list(list, this%sim%RCOMPS)
        this%sim%NRCOMP = size(this%sim%RCOMPS)
        list => noise%get_list('SCOMPS', required=.true., error=io_err)
        if(associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        call read_string_list(list, this%sim%SCOMPS)
        this%sim%NSCOMP = size(this%sim%SCOMPS)
        this%sim%CH_CODE = noise%get_string('CH_CODE', error=io_err)
        this%sim%DT = noise%get_real('DT', error=io_err)
        list => noise%get_list('SHORT_P', required=.true., error=io_err)
        if(associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        call read_real_list(list, this%sim%SHORT_P)
        list => noise%get_list('LONG_P', required=.true., error=io_err)
        if(associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        call read_real_list(list, this%sim%LONG_P)
        list => noise%get_list('GROUPVEL_MIN', required=.true., error=io_err)
        if(associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        call read_real_list(list, this%sim%GROUPVEL_MIN)
        list => noise%get_list('GROUPVEL_MAX', required=.true., error=io_err)
        if(associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        call read_real_list(list, this%sim%GROUPVEL_MAX)
        if (size(this%sim%SHORT_P)      == size(this%sim%LONG_P) .and. &
            size(this%sim%LONG_P)       == size(this%sim%GROUPVEL_MIN) .and. &
            size(this%sim%GROUPVEL_MIN) == size(this%sim%GROUPVEL_MAX)) then
          this%sim%NUM_FILTER = size(this%sim%SHORT_P)
        else
          call exit_mpi(worldrank, 'ERROR: the number of filters for noise FWI is not consistent')
        endif
        this%sim%USE_NEAR_OFFSET = noise%get_logical('USE_NEAR_OFFSET', error=io_err)
        this%sim%ADJ_SRC_NORM = noise%get_logical('ADJ_SRC_NORM', error=io_err)
        this%sim%SUPPRESS_EGF = noise%get_logical('SUPPRESS_EGF', error=io_err)

        ! read parameters for teleseismic FWI
        this%sim => tele_par
        this%sim%IMEAS = 2
        this%sim%ITAPER = 2
        tele => root%get_dictionary('TELE', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        list => tele%get_list('RCOMPS', required=.true., error=io_err)
        if(associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        call read_string_list(list, this%sim%RCOMPS)
        this%sim%NRCOMP = size(this%sim%RCOMPS)
        this%sim%NSTEP = tele%get_integer('NSTEP', error=io_err)
        this%sim%PRECOND_TYPE = tele%get_integer('PRECOND_TYPE', error=io_err)
        this%sim%CH_CODE = tele%get_string('CH_CODE', error=io_err)
        this%sim%DT = tele%get_real('DT', error=io_err)
        list => tele%get_list('TIME_WIN', required=.true., error=io_err)
        if(associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        call read_real_list(list, this%sim%TIME_WIN)
        list => tele%get_list('SHORT_P', required=.true., error=io_err)
        if(associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        call read_real_list(list, this%sim%SHORT_P)
        list => tele%get_list('LONG_P', required=.true., error=io_err)
        if(associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        call read_real_list(list, this%sim%LONG_P)
        this%sim%NUM_FILTER = 1
        this%sim%USE_LOCAL_STF = tele%get_logical('USE_LOCAL_STF', error=io_err)
        ! read parameters for RF proc
        rf => tele%get_dictionary('RF', required=.true., error=io_err)
        this%sim%rf%MINDERR = rf%get_real('MINDERR', error=io_err)
        this%sim%rf%TSHIFT = rf%get_real('TSHIFT', error=io_err)
        this%sim%rf%MAXIT = rf%get_integer('MAXIT', error=io_err)
        list => rf%get_list('F0', required=.true., error=io_err)
        if (associated(io_err)) then
          this%sim%rf%NGAUSS = 0
        else
          call read_real_list(list, this%sim%rf%F0)
          this%sim%rf%NGAUSS = size(this%sim%rf%F0)
        endif
      end select
      call root%finalize()
      deallocate(root)
    endif
    call synchronize_all()

    ! broadcast noise parameters
    call bcast_all_singlei(noise_par%NSTEP)
    call bcast_all_singlei(noise_par%IMEAS)
    call bcast_all_singlei(noise_par%ITAPER)
    call bcast_all_singlei(noise_par%PRECOND_TYPE)
    call bcast_all_singlei(noise_par%NRCOMP)
    call bcast_all_singlei(noise_par%NSCOMP)
    call bcast_all_singlei(noise_par%NUM_FILTER)
    call bcast_all_singlel(noise_par%USE_NEAR_OFFSET)
    call bcast_all_singlel(noise_par%ADJ_SRC_NORM)
    call bcast_all_singlel(noise_par%SUPPRESS_EGF)
    if (worldrank > 0) then
      allocate(noise_par%RCOMPS(noise_par%NRCOMP))
      allocate(noise_par%SCOMPS(noise_par%NSCOMP))
      allocate(noise_par%SHORT_P(noise_par%NUM_FILTER))
      allocate(noise_par%LONG_P(noise_par%NUM_FILTER))
      allocate(noise_par%GROUPVEL_MIN(noise_par%NUM_FILTER))
      allocate(noise_par%GROUPVEL_MAX(noise_par%NUM_FILTER))
    endif
    call bcast_all_ch_array(noise_par%RCOMPS, noise_par%NRCOMP, MAX_STRING_LEN)
    call bcast_all_ch_array(noise_par%SCOMPS, noise_par%NSCOMP, MAX_STRING_LEN)
    call bcast_all_ch_array(noise_par%CH_CODE, 1, MAX_STRING_LEN)
    call bcast_all_singlecr(noise_par%DT)
    call bcast_all_r(noise_par%SHORT_P, noise_par%NUM_FILTER)
    call bcast_all_r(noise_par%LONG_P, noise_par%NUM_FILTER)
    call bcast_all_r(noise_par%GROUPVEL_MIN, noise_par%NUM_FILTER)
    call bcast_all_r(noise_par%GROUPVEL_MAX, noise_par%NUM_FILTER)

    ! broadcast tele parameters
    call bcast_all_singlei(tele_par%NSTEP)
    call bcast_all_singlei(tele_par%IMEAS)
    call bcast_all_singlei(tele_par%ITAPER)
    call bcast_all_singlei(tele_par%PRECOND_TYPE)
    call bcast_all_singlei(tele_par%NRCOMP)
    call bcast_all_singlei(tele_par%NUM_FILTER)
    call bcast_all_singlel(tele_par%USE_LOCAL_STF)
    call bcast_all_singlei(tele_par%rf%MAXIT)
    call bcast_all_singlecr(tele_par%rf%MINDERR)
    call bcast_all_singlecr(tele_par%rf%TSHIFT)
    call bcast_all_singlei(tele_par%rf%NGAUSS)
    if (worldrank > 0) then
      allocate(tele_par%RCOMPS(tele_par%NRCOMP))
      allocate(tele_par%TIME_WIN(2))
      allocate(tele_par%SHORT_P(1))
      allocate(tele_par%LONG_P(1))
      if (tele_par%rf%NGAUSS > 0) then
        allocate(tele_par%rf%F0(tele_par%rf%NGAUSS))
      endif
    endif
    call bcast_all_ch_array(tele_par%RCOMPS, tele_par%NRCOMP, MAX_STRING_LEN)
    call bcast_all_ch_array(tele_par%CH_CODE, 1, MAX_STRING_LEN)
    call bcast_all_singlecr(tele_par%DT)
    call bcast_all_r(tele_par%TIME_WIN, 2)
    call bcast_all_r(tele_par%SHORT_P, 1)
    call bcast_all_r(tele_par%LONG_P, 1)
    if (tele_par%rf%NGAUSS > 0) then
      call bcast_all_r(tele_par%rf%F0, tele_par%rf%NGAUSS)
    endif

    call synchronize_all()

  end subroutine read_parameter_file

  subroutine select_simu_type(this)
    use specfem_par, only: NSTEP, DT
    use ma_variables
    class(fwat_params), intent(inout) :: this
    integer :: icomp

    select case (simu_type)
    case (SIMU_TYPE_TELE)
      this%sim => tele_par
      is_mtm0 = 0
    case (SIMU_TYPE_NOISE)
      this%sim => noise_par
      is_mtm0 = 1
    end select
    imeas0 = this%sim%IMEAS
    imeas = imeas0
    itaper = this%sim%ITAPER
    is_mtm = is_mtm0
    DT = this%sim%DT
    NSTEP = this%sim%NSTEP
    do icomp = 1, this%sim%NRCOMP
      if (trim(this%sim%RCOMPS(icomp)) == 'R' .or. trim(this%sim%RCOMPS(icomp)) == 'T') then
        dat_coord = 'ZRT'
      else
        dat_coord = 'ZNE'
      end if
    end do
  end subroutine select_simu_type

  subroutine read_string_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    character(len=MAX_STRING_LEN), dimension(:), allocatable, intent(out) :: list_out
    integer :: i, count
    
    ! count the number of items in the list
    item => list%first
    i = 0
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        i = i + 1
      end select
      item => item%next
    end do
    count = i
    allocate(list_out(count))

    item => list%first
    i = 1
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        list_out(i) = trim(element%to_string())
        item => item%next
        i = i + 1
      end select
    enddo
  end subroutine read_string_list

  subroutine read_dp_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    real(kind=dp), dimension(:), allocatable, intent(out) :: list_out
    integer :: i, count
    
    ! count the number of items in the list
    item => list%first
    i = 0
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        i = i + 1
      end select
      item => item%next
    end do
    count = i
    allocate(list_out(count))

    item => list%first
    i = 1
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        list_out(i) = element%to_real(0.)
        item => item%next
        i = i + 1
      end select
    enddo
  end subroutine read_dp_list

  subroutine read_real_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    real(kind=cr), dimension(:), allocatable, intent(out) :: list_out
    integer :: i, count

    ! count the number of items in the list
    item => list%first
    i = 0
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        i = i + 1
      end select
      item => item%next
    end do
    count = i
    allocate(list_out(count))
    
    item => list%first
    i = 1
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        list_out(i) = element%to_real(0.)
        item => item%next
        i = i + 1
      end select
    enddo
  end subroutine read_real_list

  subroutine read_logi_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    logical, dimension(:), allocatable, intent(out) :: list_out
    integer :: i, count

    ! count the number of items in the list
    item => list%first
    i = 0
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        i = i + 1
      end select
      item => item%next
    end do
    count = i
    allocate(list_out(count))
    
    item => list%first
    i = 1
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        list_out(i) = element%to_logical(.true.)
        item => item%next
        i = i + 1
      end select
    enddo
  end subroutine read_logi_list

  subroutine read_i_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    integer, dimension(:), intent(out) :: list_out
    integer :: i
    
    item => list%first
    i = 1
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        list_out(i) = element%to_integer(0)
        item => item%next
        i = i + 1
      end select
    enddo
  end subroutine read_i_list

end module input_params