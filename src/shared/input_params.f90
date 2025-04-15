module input_params
  use config
  use fwat_mpi
  use ma_variables

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
    procedure :: read => acqui_read_source_set, finalize => acqui_finalize
  end type acqui_params

  type rf_params
    real(kind=cr)                                    :: MINDERR, TSHIFT
    integer                                          :: NGAUSS, MAXIT
    real(kind=cr), dimension(:), allocatable :: F0
  end type rf_params

  type sim_params
    character(len=MAX_STRING_LEN) :: mesh_par_file
    integer :: NRCOMP, NSCOMP, NUM_FILTER, NSTEP, IMEAS, ITAPER, PRECOND_TYPE, TELE_TYPE
    character(len= MAX_STRING_LEN), dimension(:), allocatable :: RCOMPS, SCOMPS
    character(len= MAX_STRING_LEN) :: CH_CODE
    real(kind=cr) :: DT, SIGMA_H, SIGMA_V
    real(kind=cr), dimension(:), allocatable :: SHORT_P, LONG_P, GROUPVEL_MIN, GROUPVEL_MAX, TIME_WIN
    logical :: USE_NEAR_OFFSET, ADJ_SRC_NORM, SUPPRESS_EGF, USE_LOCAL_STF, USE_RHO_SCALING, SAVE_FK
    type(rf_params) :: rf
  end type sim_params

  type model_grid
    integer :: regular_grid_size(3)
    real(kind=cr) :: regular_grid_min_coord(3), regular_grid_interval(3)
  end type model_grid

  type postproc_params
    logical, dimension(2) :: INV_TYPE
    integer :: SMOOTH_TYPE
    real(kind=cr), dimension(2) :: JOINT_WEIGHT
    real(kind=cr) :: TAPER_H_SUPPRESS, TAPER_V_SUPPRESS, TAPER_H_BUFFER, TAPER_V_BUFFER
    integer, dimension(3) :: ninv
    integer :: n_inversion_grid
    logical :: IS_PRECOND
  end type postproc_params

  type update_params
    integer :: MODEL_TYPE, ITER_START, LBFGS_M_STORE, OPT_METHOD, MAX_SUB_ITER
    real(kind=cr) :: MAX_SLEN, MAX_SHRINK
    logical :: DO_LS
    real(kind=cr), dimension(2) :: VPVS_RATIO_RANGE
  end type update_params

  type fwat_params
    type(acqui_params) :: acqui
    type(sim_params), pointer :: sim
    type(postproc_params) :: postproc
    type(update_params) :: update
    type(model_grid) :: grid
    contains
    procedure :: read => read_fwat_parameter_file, select_simu_type
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
        read(FID,*,iostat=ier) evtnm,this%evla(ievt),this%evlo(ievt), &
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

  subroutine acqui_finalize(this)
    class(acqui_params), intent(inout) :: this

    ! deallocate arrays on shared memory
    call free_shm_array(this%sta_win)
    call free_shm_array(this%src_win)
    call free_shm_array(this%evid_win)
    call free_shm_array(this%ofp_win)
    call free_shm_array(this%idp_win)
    if (simu_type == SIMU_TYPE_TELE) then
      call free_shm_array(this%fk_win)
    endif
    call free_shm_array(this%evla_win)
    call free_shm_array(this%evlo_win)
    call free_shm_array(this%evdp_win)
    call free_shm_array(this%w_win)

  end subroutine acqui_finalize

  subroutine read_fwat_parameter_file(this, fname)
    use yaml, only: parse, error_length
    use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                        type_list, type_list_item, type_scalar
    class(fwat_params), intent(inout) :: this
    character(len=*), intent(in) :: fname
    class(type_node), pointer :: root
    class(type_dictionary), pointer :: noise, tele, tomo, rf, output, post, update, ma, grid
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
        this%sim%mesh_par_file = noise%get_string('MESH_PAR_FILE', error=io_err)
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
        this%sim%SIGMA_H = noise%get_real('SIGMA_H', error=io_err)
        this%sim%SIGMA_V = noise%get_real('SIGMA_V', error=io_err)
        this%sim%USE_RHO_SCALING = noise%get_logical('USE_RHO_SCALING', error=io_err, default=.true.)

        ! read parameters for teleseismic FWI
        this%sim => tele_par
        this%sim%IMEAS = 2
        this%sim%ITAPER = 2
        this%sim%USE_RHO_SCALING = .false.
        tele => root%get_dictionary('TELE', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        this%sim%mesh_par_file = tele%get_string('MESH_PAR_FILE', error=io_err)
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
        this%sim%SIGMA_H = tele%get_real('SIGMA_H', error=io_err)
        this%sim%SIGMA_V = tele%get_real('SIGMA_V', error=io_err)
        this%sim%TELE_TYPE = tele%get_integer('TELE_TYPE', error=io_err)
        this%sim%SAVE_FK = tele%get_logical('SAVE_FK', error=io_err, default=.true.)
        compress_level = tele%get_integer('COMPRESS_LEVEL', error=io_err, default=0)
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

        ! measure adjoint source
        ma => root%get_dictionary('MEASURE_ADJ', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        TSHIFT_MIN = ma%get_real('TSHIFT_MIN', error=io_err, default=4.5)
        TSHIFT_MAX = ma%get_real('TSHIFT_MAX', error=io_err, default=4.5)
        DLNA_MIN = ma%get_real('DLNA_MIN', error=io_err)
        DLNA_MAX = ma%get_real('DLNA_MAX', error=io_err)
        CC_MIN = ma%get_real('CC_MIN', error=io_err)
        ERROR_TYPE = ma%get_integer('ERROR_TYPE', error=io_err)
        DT_SIGMA_MIN = ma%get_real('DT_SIGMA_MIN', error=io_err)
        DLNA_SIGMA_MIN = ma%get_real('DLNA_SIGMA_MIN', error=io_err)
        WTR = ma%get_real('WTR', error=io_err)
        NPI = ma%get_real('NPI', error=io_err)
        DT_FAC = ma%get_real('DT_FAC', error=io_err)
        ERR_FAC = ma%get_real('ERR_FAC', error=io_err)
        DT_MAX_SCALE = ma%get_real('DT_MAX_SCALE', error=io_err)
        NCYCLE_IN_WINDOW = ma%get_real('NCYCLE_IN_WINDOW', error=io_err)
        USE_PHYSICAL_DISPERSION = ma%get_logical('USE_PHYSICAL_DISPERSION', error=io_err, default=.false.)

        ! output
        output => root%get_dictionary('OUTPUT', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        is_output_preproc = output%get_logical('IS_OUTPUT_PREPROC', error=io_err, default=.false.)
        is_output_adj_src = output%get_logical('IS_OUTPUT_ADJ_SRC', error=io_err, default=.false.)
        is_output_event_kernel = output%get_logical('IS_OUTPUT_EVENT_KERNEL', error=io_err, default=.false.)
        is_output_sum_kernel = output%get_logical('IS_OUTPUT_SUM_KERNEL', error=io_err, default=.false.)
        is_output_hess_inv = output%get_logical('IS_OUTPUT_HESS_INV', error=io_err, default=.false.)
        is_output_inv_grid = output%get_logical('IS_OUTPUT_INV_GRID', error=io_err, default=.false.)

        ! model grid
        grid => root%get_dictionary('MODEL_GRID', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        list => grid%get_list('REGULAR_GRID_SIZE', required=.true., error=io_err)
        call read_static_int_list(list, this%grid%regular_grid_size)
        list => grid%get_list('REGULAR_GRID_MIN_COORD', required=.true., error=io_err)
        call read_static_real_list(list, this%grid%regular_grid_min_coord)
        list => grid%get_list('REGULAR_GRID_INTERVAL', required=.true., error=io_err)
        call read_static_real_list(list, this%grid%regular_grid_interval)

        ! POSTPROC
        post => root%get_dictionary('POSTPROC', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        list => post%get_list('INV_TYPE', required=.true., error=io_err)
        call read_static_logi_list(list, this%postproc%INV_TYPE)
        list => post%get_list('JOINT_WEIGHT', required=.true., error=io_err)
        call read_static_real_list(list, this%postproc%JOINT_WEIGHT)
        list => post%get_list('NINV', required=.true., error=io_err)
        call read_static_int_list(list, this%postproc%ninv)
        this%postproc%n_inversion_grid = post%get_integer('N_INVERSION_GRID', error=io_err, default=5)
        this%postproc%TAPER_H_SUPPRESS = post%get_real('TAPER_H_SUPPRESS', error=io_err)
        this%postproc%TAPER_V_SUPPRESS = post%get_real('TAPER_V_SUPPRESS', error=io_err)
        this%postproc%TAPER_H_BUFFER = post%get_real('TAPER_H_BUFFER', error=io_err)
        this%postproc%TAPER_V_BUFFER = post%get_real('TAPER_V_BUFFER', error=io_err)
        this%postproc%IS_PRECOND = post%get_logical('IS_PRECOND', error=io_err)
        this%postproc%SMOOTH_TYPE = post%get_integer('SMOOTH_TYPE', error=io_err, default=1)

        ! Model UPDATE
        update => root%get_dictionary('MODEL_UPDATE', required=.true., error=io_err)
        if (associated(io_err)) call exit_mpi(worldrank, trim(io_err%message))
        this%update%MODEL_TYPE = update%get_integer('MODEL_TYPE', error=io_err)
        this%update%ITER_START = update%get_integer('ITER_START', error=io_err)
        this%update%LBFGS_M_STORE = update%get_integer('LBFGS_M_STORE', error=io_err)
        this%update%OPT_METHOD = update%get_integer('OPT_METHOD', error=io_err)
        this%update%MAX_SLEN = update%get_real('MAX_SLEN', error=io_err)
        this%update%MAX_SHRINK = update%get_real('MAX_SHRINK', error=io_err)
        this%update%MAX_SUB_ITER = update%get_integer('MAX_SUB_ITER', error=io_err)
        this%update%DO_LS = update%get_logical('DO_LS', error=io_err)
        list => update%get_list('VPVS_RATIO_RANGE', required=.true., error=io_err)
        call read_static_real_list(list, this%update%VPVS_RATIO_RANGE)

      end select
      call root%finalize()
      deallocate(root)
    endif
    call synchronize_all()

    ! broadcast noise parameters
    call bcast_all_ch_array(noise_par%mesh_par_file, 1, MAX_STRING_LEN)
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
    call bcast_all_singlel(noise_par%USE_RHO_SCALING)
    call bcast_all_singlecr(noise_par%SIGMA_H)
    call bcast_all_singlecr(noise_par%SIGMA_V)
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
    call bcast_all_ch_array(tele_par%mesh_par_file, 1, MAX_STRING_LEN)
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
    call bcast_all_singlecr(tele_par%SIGMA_H)
    call bcast_all_singlecr(tele_par%SIGMA_V)
    call bcast_all_singlel(tele_par%USE_RHO_SCALING)
    call bcast_all_singlel(tele_par%SAVE_FK)
    call bcast_all_singlei(compress_level)
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
    call bcast_all_singlei(tele_par%TELE_TYPE)
    if (tele_par%rf%NGAUSS > 0) then
      call bcast_all_r(tele_par%rf%F0, tele_par%rf%NGAUSS)
    endif

    ! measure adjoint source
    call bcast_all_singledp(TSHIFT_MIN)
    call bcast_all_singledp(TSHIFT_MAX)
    call bcast_all_singledp(DLNA_MIN)
    call bcast_all_singledp(DLNA_MAX)
    call bcast_all_singledp(CC_MIN)
    call bcast_all_singlei(ERROR_TYPE)
    call bcast_all_singledp(DT_SIGMA_MIN)
    call bcast_all_singledp(DLNA_SIGMA_MIN)
    call bcast_all_singledp(WTR)
    call bcast_all_singledp(NPI)
    call bcast_all_singledp(DT_FAC)
    call bcast_all_singledp(ERR_FAC)
    call bcast_all_singledp(DT_MAX_SCALE)
    call bcast_all_singledp(NCYCLE_IN_WINDOW)
    call bcast_all_singlel(USE_PHYSICAL_DISPERSION)

    ! output
    call bcast_all_singlel(IS_OUTPUT_PREPROC)
    call bcast_all_singlel(IS_OUTPUT_ADJ_SRC)
    call bcast_all_singlel(IS_OUTPUT_EVENT_KERNEL)
    call bcast_all_singlel(IS_OUTPUT_SUM_KERNEL)
    call bcast_all_singlel(IS_OUTPUT_HESS_INV)
    call bcast_all_singlel(IS_OUTPUT_INV_GRID)

    ! model grid
    call bcast_all_i(this%grid%regular_grid_size, 3)
    call bcast_all_r(this%grid%regular_grid_min_coord, 3)
    call bcast_all_r(this%grid%regular_grid_interval, 3)

    ! POSTPROC
    call bcast_all_singlecr(this%postproc%TAPER_H_SUPPRESS)
    call bcast_all_singlecr(this%postproc%TAPER_V_SUPPRESS)
    call bcast_all_singlecr(this%postproc%TAPER_H_BUFFER)
    call bcast_all_singlecr(this%postproc%TAPER_V_BUFFER)
    call bcast_all_singlel(this%postproc%IS_PRECOND)
    call bcast_all_l_array(this%postproc%INV_TYPE, 2)
    call bcast_all_r(this%postproc%JOINT_WEIGHT, 2)
    call bcast_all_i(this%postproc%ninv, 3)
    call bcast_all_singlei(this%postproc%n_inversion_grid)
    call bcast_all_singlei(this%postproc%SMOOTH_TYPE)

    ! Model UPDATE
    call bcast_all_singlei(this%update%MODEL_TYPE)
    call bcast_all_singlei(this%update%ITER_START)
    call bcast_all_singlei(this%update%LBFGS_M_STORE)
    call bcast_all_singlei(this%update%OPT_METHOD)
    call bcast_all_singlecr(this%update%MAX_SLEN)
    call bcast_all_singlecr(this%update%MAX_SHRINK)
    call bcast_all_singlei(this%update%MAX_SUB_ITER)
    call bcast_all_singlel(this%update%DO_LS)
    call bcast_all_r(this%update%VPVS_RATIO_RANGE, 2)
    call synchronize_all()

  end subroutine read_fwat_parameter_file

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
    ! do icomp = 1, this%sim%NRCOMP
    !   if (trim(this%sim%RCOMPS(icomp)) == 'R' .or. trim(this%sim%RCOMPS(icomp)) == 'T') then
    !     dat_coord = 'ZRT'
    !   else
    !     dat_coord = 'ZNE'
    !   end if
    ! end do
  end subroutine select_simu_type

  subroutine read_static_int_list(list, list_out)
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
  end subroutine read_static_int_list

  subroutine read_static_real_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    real(kind=cr), dimension(:), intent(out) :: list_out
    integer :: i
    
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
  end subroutine read_static_real_list

  subroutine read_static_logi_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    logical, dimension(:), intent(out) :: list_out
    integer :: i
    
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
  end subroutine read_static_logi_list

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