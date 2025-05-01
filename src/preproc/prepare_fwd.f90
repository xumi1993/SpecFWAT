module preproc_fwd
  use config
  use fwat_mpi
  use specfem_par, only: ACOUSTIC_SIMULATION, APPROXIMATE_HESS_KL, ELASTIC_SIMULATION, MOVIE_VOLUME, NGLOB_AB,&
                         PRINT_SOURCE_TIME_FUNCTION, SAVE_FORWARD, SAVE_MESH_FILES, SIMULATION_TYPE, LOCAL_PATH,&
                         ANISOTROPIC_KL, NSPEC_ADJOINT, FKMODEL_FILE, prname, COUPLE_WITH_INJECTION_TECHNIQUE,&
                         USE_FORCE_POINT_SOURCE, OUTPUT_FILES, INJECTION_TECHNIQUE_TYPE,ADIOS_FOR_MESH, &
                         ATTENUATION
  use specfem_par_elastic, only: rmass, rmassx, rmassy, rmassz, COMPUTE_AND_STORE_STRAIN
  use specfem_par_acoustic, only: rmass_acoustic
  ! use specfem_par_poroelastic
  use input_params, fpar => fwat_par_global
  use fk_coupling, only: couple_with_injection_prepare_boundary_fwat, check_fk_files, read_fk_model
  use logger, only: log
  use tele_data, only: TeleData
  use rf_data, only: RFData
  use telecc_data, only: TeleCCData
  use noise_data, only: NoiseData
  use common_lib, only: get_dat_type, mkdir

  implicit none

  integer, private :: ier

  type :: PrepareFWD
    integer :: ievt=0, run_mode
    real(kind=dp) :: obj_func
    contains
    procedure :: init, prepare_for_event, destroy, simulation
    procedure, private :: initialize_kernel_matrice, semd2sac, measure_adj, run_simulation,&
                          postproc_adjoint
  end type PrepareFWD

contains

  subroutine init(this, is_init_log)
    use specfem_api, only: backup_rmass

    logical, optional, intent(in) :: is_init_log
    logical :: is_init_log_loc
    class(PrepareFWD), intent(inout) :: this

    this%run_mode = run_mode

    is_init_log_loc = .true.
    if (present(is_init_log)) then
      is_init_log_loc = is_init_log
    endif
    if (is_init_log_loc) then
      if (single_run) then
        ! if (worldrank == 0) call system('mkdir -p '//trim(fpar%acqui%out_fwd_path(this%ievt)))
        call mkdir(trim(fpar%acqui%out_fwd_path(this%ievt)))
        call log%init(trim(fpar%acqui%out_fwd_path(this%ievt))//'/output_fwd_measure_adj.log')
      else
        call log%init('output_fwd_measure_adj_'//trim(dat_type)//'_'//trim(model_name)//'.log')
      endif
    endif
    call log%write('*******************************************', .false.)

    call force_ftz()

    ! reads in parameters
    call initialize_simulation_fwat()
    
    call read_mesh_databases_fwat()

    ! copy mass matrixs for prepare_timerun()
    call backup_rmass()

    ! reads in moho mesh
    if (ADIOS_FOR_MESH) then
      call read_mesh_databases_moho_adios()
    else
      call read_mesh_databases_moho_fwat()
    endif

    call read_mesh_databases_adjoint_fwat()

    ! call couple_with_injection_setup()

    ! sets up reference element GLL points/weights/derivatives
    call setup_GLL_points()

    call detect_mesh_surfaces()

    if(this%run_mode == FORWARD_ADJOINT) call this%initialize_kernel_matrice()

  end subroutine init

  subroutine destroy(this)
    class(PrepareFWD), intent(inout) :: this

    if (allocated(rmass_acoustic_copy)) deallocate(rmass_acoustic_copy)
    if (allocated(rmass_copy)) deallocate(rmass_copy)
    if (allocated(rmassx_copy)) deallocate(rmassx_copy)
    if (allocated(rmassy_copy)) deallocate(rmassy_copy)
    if (allocated(rmassz_copy)) deallocate(rmassz_copy)
    if (allocated(sum_rho_kl)) deallocate(sum_rho_kl)
    if (allocated(sum_cijkl_kl)) deallocate(sum_cijkl_kl)
    if (allocated(sum_mu_kl)) deallocate(sum_mu_kl)
    if (allocated(sum_kappa_kl)) deallocate(sum_kappa_kl)
    if (allocated(sum_hess_kl)) deallocate(sum_hess_kl)
  end subroutine destroy

  subroutine initialize_kernel_matrice(this)
    class(PrepareFWD), intent(inout) :: this
    if (ELASTIC_SIMULATION) then
      allocate(sum_rho_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      sum_rho_kl = 0.
      if (ier /= 0) call exit_MPI(worldrank, 'Error allocating array sum_rho_kl')
      if (ANISOTROPIC_KL) then
        allocate(sum_cijkl_kl(21,NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
        sum_cijkl_kl = 0.
        if (ier /= 0) call exit_MPI(worldrank, 'Error allocating array sum_cijkl_kl')
      else
        allocate(sum_mu_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
        sum_mu_kl = 0.
        if (ier /= 0) call exit_MPI(worldrank, 'Error allocating array sum_mu_kl')
        allocate(sum_kappa_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
        sum_kappa_kl = 0.
        if (ier /= 0) call exit_MPI(worldrank, 'Error allocating array sum_kappa_kl')
      endif
      if (APPROXIMATE_HESS_KL) then
        allocate(sum_hess_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
        sum_hess_kl = 0.
        if (ier /= 0) call exit_MPI(worldrank, 'Error allocating array sum_hess_kl')
      endif
    endif
  end subroutine initialize_kernel_matrice

  subroutine prepare_for_event(this)
    class(PrepareFWD), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: evtid
    character(len=MAX_STRING_LEN) :: path, source_fname, station_fname

    evtid = fpar%acqui%evtid_names(this%ievt)
    path = trim(fpar%acqui%out_fwd_path(this%ievt))
    call log%write('Creating directory for event '//trim(evtid), .true.)

    ! create directory for forward simulation
    call mkdir(trim(path))
    call mkdir(trim(path)//'/'//trim(OUTPUT_PATH))
    if (run_mode > FORWARD_ONLY) call mkdir(trim(path)//'/'//trim(ADJOINT_PATH))
    if (run_mode > FORWARD_MEASADJ) call mkdir(trim(path)//'/'//trim(EKERNEL_PATH))
    if (run_mode >= FORWARD_MEASADJ) call mkdir('misfits')
    call synchronize_all()

    ! assign local path
    OUTPUT_FILES=trim(path)//'/OUTPUT_FILES'
    LOCAL_PATH=trim(path)//'/EKERNEL'
    call create_name_database(prname,worldrank,LOCAL_PATH)

    ! change simulation parameters for simu_type
    if (simu_type == SIMU_TYPE_TELE) then
      USE_FORCE_POINT_SOURCE = .false.
      COUPLE_WITH_INJECTION_TECHNIQUE=.true.
      INJECTION_TECHNIQUE_TYPE=3
      FKMODEL_FILE = fpar%acqui%fkmodel_file(this%ievt)
      if (dat_type == 'rf') then
        fpar%sim%NUM_FILTER = fpar%sim%rf%NGAUSS
      endif
    else if (simu_type == SIMU_TYPE_NOISE) then
      USE_FORCE_POINT_SOURCE = .true.
      COUPLE_WITH_INJECTION_TECHNIQUE=.false.
    endif
    source_fname = fpar%acqui%src_solution_file(this%ievt)
    station_fname = fpar%acqui%station_file(this%ievt)

    ! print info
    block 
      character(len=MAX_STRING_LEN) :: msg
      write(msg, '(a,I4,a)') 'Index:',this%ievt,', Event ID: '//trim(evtid)
      call log%write(msg)
      msg = 'Source file: '//trim(source_fname)//', Station file: '//trim(station_fname)
      call log%write(msg)
      write(msg, '(A, f6.4)') 'Weight: ', fpar%acqui%src_weight(this%ievt)
      call log%write(msg)
      call log%write('out_fwd_path = '//trim(path))
    end block

    call synchronize_all()

  end subroutine prepare_for_event

  subroutine semd2sac(this)
    ! This subroutine saves the simulation results to sac files
    class(PrepareFWD), intent(inout) :: this
    type(TeleData) :: td
    type(RFData) :: rd
    type(NoiseData) :: nd

    select case (dat_type)
      case ('tele') 
        call td%semd2sac(this%ievt)
      case ('telecc')
        call td%semd2sac(this%ievt)
      case ('rf')
        call rd%semd2sac(this%ievt)
      case ('noise')
        call nd%semd2sac(this%ievt)
    end select

  end subroutine semd2sac

  subroutine measure_adj(this)
    class(PrepareFWD), intent(inout) :: this
    type(TeleData) :: td
    type(RFData) :: rd
    type(TeleCCData) :: tc
    type(NoiseData) :: nd
    
    select case(dat_type)
    case ('tele')
      call td%preprocess(this%ievt)
      call td%od%copy_adjoint_stations()
      this%obj_func = sum(td%total_misfit)
      call td%finalize()
    case ('telecc')
      call tc%preprocess(this%ievt)
      call tc%od%copy_adjoint_stations()
      this%obj_func = sum(tc%total_misfit)
      call tc%finalize()
    case ('rf')
      call rd%preprocess(this%ievt)
      call rd%od%copy_adjoint_stations()
      this%obj_func = sum(rd%total_misfit)
      call rd%finalize()
    case ('noise')
      call nd%preprocess(this%ievt)
      call nd%od%copy_adjoint_stations()
      this%obj_func = sum(nd%total_misfit)
      call nd%finalize()
    end select
    call bcast_all_singledp(this%obj_func)
  end subroutine measure_adj

  subroutine simulation(this)
    class(PrepareFWD), intent(inout) :: this
    type(TeleData) :: td

    call log%write('This is forward simulations ...', .true.)
    
    ! run forward simulation
    call this%run_simulation(1)

    ! preprocess data
    if (this%run_mode == FORWARD_ONLY) then
      ! save simulation results
      call log%write('Writing synthetic data ...', .true.)
      call this%semd2sac()
    elseif (this%run_mode >= FORWARD_MEASADJ) then
      ! save simulation results
      call log%write('Measuring adjoint source ...', .true.)
      call this%measure_adj()
    endif

    ! run adjoint simulation
    if (this%run_mode == FORWARD_ADJOINT) then
      ! save adjoint source
      call log%write('This is adjoint simulations...', .true.)
      call this%run_simulation(3)
      call this%postproc_adjoint()
    endif
    OUTPUT_FILES = output_files_backup
    call log%write('-------------------------------------------', .false.)
  end subroutine simulation

  subroutine run_simulation(this, run_opt)
    use specfem_api, only: InitSpecfem, FinalizeSpecfem, restore_rmass
    class(PrepareFWD), intent(inout) :: this
    integer, intent(in) :: run_opt

    SIMULATION_TYPE = run_opt
    if (run_mode == FORWARD_ADJOINT) then
      SAVE_FORWARD = .true.
    else
      SAVE_FORWARD = .false.
    endif
    if (SIMULATION_TYPE == 1) then
      COMPUTE_AND_STORE_STRAIN = .false.
    else
      COMPUTE_AND_STORE_STRAIN = .true.
      ATTENUATION = .false.
    endif

    ! Initialize variables of Specfem
    call InitSpecfem()

    ! set up sources and receivers
    call setup_sources_receivers_fwat(this%ievt)

    ! restore mass matrixs from rmass_copy
    call restore_rmass()

    ! prepare for time run
    call prepare_timerun_fwat(this%ievt)

    ! run forward simulation
    call iterate_time()

    ! save absobing boundary wavefields
    call FinalizeSpecfem()
  end subroutine run_simulation

  subroutine postproc_adjoint(this)
    class(PrepareFWD), intent(inout) :: this

    open(unit=1234, iostat=ier, file=trim(prname)//'save_forward_arrays.bin', status='old')
    if (ier == 0) close(1234, status='delete') 
    open(unit=1234, iostat=ier, file=trim(prname)//'absorb_field.bin', status='old')
    if (ier == 0) close(1234, status='delete') 
    
    call synchronize_all()
    call log%write('This is saving kernels ...', .true.)
    call save_adjoint_kernels() 

  end subroutine postproc_adjoint

end module