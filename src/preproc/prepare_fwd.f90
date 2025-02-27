module preproc_fwd
  use config
  use fwat_mpi
  use specfem_par, only: ACOUSTIC_SIMULATION, APPROXIMATE_HESS_KL, ELASTIC_SIMULATION, MOVIE_VOLUME, NGLOB_AB,&
                         PRINT_SOURCE_TIME_FUNCTION, SAVE_FORWARD, SAVE_MESH_FILES, SIMULATION_TYPE, LOCAL_PATH,&
                         ANISOTROPIC_KL, NSPEC_ADJOINT, FKMODEL_FILE, prname, COUPLE_WITH_INJECTION_TECHNIQUE,&
                         USE_FORCE_POINT_SOURCE, OUTPUT_FILES, INJECTION_TECHNIQUE_TYPE
  use specfem_par_elastic, only: rmass, rmassx, rmassy, rmassz
  use specfem_par_acoustic, only: rmass_acoustic
  ! use specfem_par_poroelastic
  use imput_params, fpar => fwat_par_global
  use fk_coupling, only: couple_with_injection_prepare_boundary_fwat, check_fk_files
  use logger, only: logger_type

  implicit none

  integer :: ier
  type(logger_type) :: log

  type :: PrepareFWD
    integer :: ievt=0, simu_opt
    contains
    procedure :: init, calc_fk_wavefield, prepare_for_event, destroy
    procedure, private :: initialize_kernel_matrice
  end type PrepareFWD

contains

  subroutine init(this, simu_opt)
    class(PrepareFWD), intent(inout) :: this
    integer, intent(in) :: simu_opt

    this%simu_opt = simu_opt

    if (this%ievt == 0) then
      call log%write('ERROR: Event index not set', .true.)
      call exit_MPI(0, 'ERROR: Event index not set')
    else if (this%ievt > fpar%acqui%nevents) then
      call log%write('ERROR: Event index out of range', .true.)
      call exit_MPI(0, 'ERROR: Event index out of range')
    endif

    if (single_run) then
      if (worldrank == 0) call system('mkdir -p '//trim(fpar%acqui%out_fwd_path(this%ievt)))
      call log%init(trim(fpar%acqui%out_fwd_path(this%ievt))//'/output_fwd_measure_adj.log')
    else
      call log%init('output_fwd_measure_adj.log')
    endif

    call force_ftz()

    ! reads in parameters
    call initialize_simulation_fwat()
    
    SAVE_MESH_FILES=.false.
    PRINT_SOURCE_TIME_FUNCTION=.false.
    MOVIE_VOLUME=.false.
    local_path_backup = LOCAL_PATH
    if (simu_opt < 3) then
      SIMULATION_TYPE = 1
      if (simu_opt == 1) then
        SAVE_FORWARD = .true.
      else
        SAVE_FORWARD = .false.
      endif
      APPROXIMATE_HESS_KL=.false.
    else
      SIMULATION_TYPE = 3
      APPROXIMATE_HESS_KL=.true.
    endif

    ! copy mass matrixs for prepare_timerun()
    if (ACOUSTIC_SIMULATION) then
      allocate(rmass_acoustic_copy(NGLOB_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating array rmass_acoustic_copy'
      rmass_acoustic_copy=rmass_acoustic
    endif
    if (ELASTIC_SIMULATION) then
      allocate(rmass_copy(NGLOB_AB),stat=ier)
      allocate(rmassx_copy(NGLOB_AB),stat=ier)
      allocate(rmassy_copy(NGLOB_AB),stat=ier)
      allocate(rmassz_copy(NGLOB_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating array rmass_copy'
      rmass_copy=rmass
      rmassx_copy=rmassx
      rmassy_copy=rmassy
      rmassz_copy=rmassz
    endif

    call read_mesh_databases_fwat()

    ! reads adjoint parameters
    call read_mesh_databases_adjoint_fwat()

    ! sets up reference element GLL points/weights/derivatives
    call setup_GLL_points()

    call this%initialize_kernel_matrice()

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

  subroutine calc_fk_wavefield(this)
    class(PrepareFWD), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: evtid
    integer :: iev

    evtid = fpar%acqui%evtid_names(this%ievt)
    ! do iev = 1, fpar%acqui%nevents
    if (.not. check_fk_files(evtid)) then
      call log%write('Calculating FK wavefield for event '//trim(evtid), .true.)
      call couple_with_injection_prepare_boundary_fwat(evtid)
    else
      call log%write('FK wavefield already calculated for event '//trim(evtid), .true.)
    endif
    ! enddo

  end subroutine calc_fk_wavefield

  subroutine prepare_for_event(this)
    class(PrepareFWD), intent(inout) :: this
    character(len=MAX_STRING_LEN) :: evtid
    character(len=MAX_STRING_LEN) :: path, source_fname, station_fname

    evtid = fpar%acqui%evtid_names(this%ievt)
    path = trim(fpar%acqui%out_fwd_path(this%ievt))
    call log%write('Creating directory for event '//trim(evtid), .true.)

    ! create directory for forward simulation
    if (worldrank == 0) then
      call system('mkdir -p '//trim(path))
      call system('mkdir -p '//trim(path)//'/'//trim(OUTPUT_PATH))
      if (this%simu_opt > 1) call system('mkdir -p '//trim(path)//'/'//trim(ADJOINT_PATH))
      if (this%simu_opt > 2) call system('mkdir -p '//trim(path)//'/'//trim(EKERNEL_PATH))
    endif

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
    else if (simu_type == SIMU_TYPE_NOISE) then
      USE_FORCE_POINT_SOURCE = .true.
      COUPLE_WITH_INJECTION_TECHNIQUE=.false.
    ! else
    !   call exit_MPI(worldrank, 'Unknown simulation type')
    endif
    source_fname=fpar%acqui%src_solution_file(this%ievt)
    station_fname=fpar%acqui%station_file(this%ievt)

    call synchronize_all()

  end subroutine prepare_for_event


end module