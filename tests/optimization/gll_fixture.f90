! Link-time fixtures replace only configuration, mesh loading and wave solves.
! The command parser, GLL optimizer, MPI reductions and binary I/O remain real.
module gll_fixture
  use config
  implicit none
  character(len=32) :: scenario
contains
  subroutine start_mpi()
    use fwat_mpi, only: init_mpi_single_group
    call init_mpi_single_group()
  end subroutine

  subroutine read_config()
    use input_params, only: fpar => fwat_par_global, tele_par
    call get_environment_variable('GLL_TEST_CASE', scenario)
    if (.not. use_gll) then
      print *, 'regular-grid mode selected'
      call finalize_mpi()
      stop
    endif
    fpar%postproc%INV_TYPE = [.false., .true., .false.]
    fpar%postproc%IS_HESS_PRECOND = scenario /= 'precond'
    fpar%update%MODEL_TYPE = 1
    fpar%update%ITER_START = 0
    fpar%update%OPT_METHOD = 1
    fpar%update%LBFGS_M_STORE = 3
    fpar%update%MAX_SLEN = 0.1_cr
    fpar%update%DO_LS = scenario == 'linesearch'
    fpar%update%MAX_SUB_ITER = 3
    fpar%update%MAX_SHRINK = 0.5_cr
    fpar%update%C1 = 0.01_cr
    fpar%update%VPVS_RATIO_RANGE = [1.3_cr, 2.5_cr]
    parameter_type = 1
    is_output_direction = .true.
    tele_par%IMEAS = IMEAS_RF
    tele_par%TELE_TYPE = 2
    tele_par%DT = 0.01_cr
    tele_par%NSTEP = 10
    tele_par%rf%NGAUSS = 1
    tele_par%rf%F0 = [1.0_cr]
    if (scenario == 'lbfgs') fpar%update%OPT_METHOD = 2
    if (scenario == 'bad_method') fpar%update%OPT_METHOD = 3
    if (scenario == 'joint') fpar%postproc%INV_TYPE(1) = .true.
  end subroutine

  subroutine read_parameters()
    use specfem_par, only: LOCAL_PATH, MODEL, IMODEL, NPROC, NUMBER_OF_SIMULTANEOUS_RUNS, &
      HDF5_IO_NODES, ADIOS_ENABLED, HDF5_ENABLED, MESH_A_CHUNK_OF_THE_EARTH, &
      ANISOTROPY, ANISOTROPIC_KL, GPU_MODE, ATTENUATION
    LOCAL_PATH = 'DATABASES_MPI'
    MODEL = 'gll'
    IMODEL = IMODEL_GLL
    if (scenario == 'external') IMODEL = IMODEL_USER_EXTERNAL
    NPROC = worldsize
    if (scenario == 'bad_mpi') NPROC = worldsize + 1
    NUMBER_OF_SIMULTANEOUS_RUNS = 1
    HDF5_IO_NODES = 0
    ADIOS_ENABLED = .false.
    HDF5_ENABLED = .false.
    MESH_A_CHUNK_OF_THE_EARTH = .false.
    ANISOTROPY = .false.
    ANISOTROPIC_KL = .false.
    GPU_MODE = .false.
    ATTENUATION = .false.
  end subroutine

  subroutine read_mesh()
    use specfem_par, only: NSPEC_AB, irregular_element_number, jacobian_regular, &
      jacobianstore, kappastore, mustore, rhostore, ELASTIC_SIMULATION, &
      ACOUSTIC_SIMULATION, POROELASTIC_SIMULATION
    use optimize_gll, only: write_gll_vector, gll_model_path
    use common_lib, only: mkdir
    real(kind=cr) :: model(NGLLX,NGLLY,NGLLZ,1,3), grad(NGLLX,NGLLY,NGLLZ,1,3)
    real(kind=cr) :: hess(NGLLX,NGLLY,NGLLZ,1,1)
    integer :: ipar, unit
    character(len=MAX_STRING_LEN) :: names(3), path
    NSPEC_AB = 1
    allocate(irregular_element_number(1), jacobianstore(NGLLX,NGLLY,NGLLZ,1))
    irregular_element_number = 0
    jacobian_regular = 0.125_cr
    ELASTIC_SIMULATION = .true.
    ACOUSTIC_SIMULATION = .false.
    POROELASTIC_SIMULATION = .false.
    model(:,:,:,:,1) = 6000.0_cr
    model(:,:,:,:,2) = 3500.0_cr
    model(:,:,:,:,3) = 2700.0_cr
    allocate(kappastore(NGLLX,NGLLY,NGLLZ,1), mustore(NGLLX,NGLLY,NGLLZ,1), &
             rhostore(NGLLX,NGLLY,NGLLZ,1))
    rhostore = model(:,:,:,:,3)
    mustore = rhostore*model(:,:,:,:,2)**2
    kappastore = rhostore*model(:,:,:,:,1)**2 - (4.0_cr/3.0_cr)*mustore
    do ipar = 1, 3
      grad(:,:,:,:,ipar) = real(ipar*(worldrank+1),cr)
      names(ipar) = trim(kernel_names(ipar))//'_kernel_smooth'
    enddo
    if (scenario == 'zero') grad = 0.0_cr
    call write_gll_vector(local_path_fwat, parameter_names, model)
    path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)
    call write_gll_vector(path, names, grad)
    hess = 2.0_cr
    if (worldrank == 1) hess = 0.5_cr
    call write_gll_vector(path, [HESS_PREFIX], hess)
    if (scenario == 'lbfgs') then
      call write_gll_vector(gll_model_path(0), parameter_names, model*exp(-0.02_cr))
      call write_gll_vector(trim(OPT_DIR)//'/SUM_KERNELS_M00', names, grad-0.25_cr)
    endif
    if (scenario == 'linesearch') then
      do ipar = 1, 3
        names(ipar) = trim(kernel_names(ipar))//'_kernel'
      enddo
      call write_gll_vector(path, names, grad)
      call mkdir(MISFITS_DIR)
      if (worldrank == 0) then
        open(newunit=unit, file=trim(MISFITS_DIR)//'/M00.event_F1.0_window_chi', status='replace')
        write(unit, '(A)') 'a b c d e f g 1.25'
        write(unit, '(A)') ''
        write(unit, '(A)') 'a b c d e f g 2.75'
        close(unit)
      endif
      call synchronize_all()
    endif
  end subroutine

  subroutine read_events()
    use input_params, only: fpar => fwat_par_global
    fpar%acqui%nevents = 1
    allocate(fpar%acqui%evtid_names(1))
    fpar%acqui%evtid_names = 'event'
  end subroutine

  subroutine free_events()
    use input_params, only: fpar => fwat_par_global
    deallocate(fpar%acqui%evtid_names)
  end subroutine

  subroutine generate_databases()
  end subroutine

  subroutine check_resolution()
  end subroutine

  subroutine forward(total, previous, current)
    real(kind=dp), intent(out) :: total, previous
    real(kind=dp), intent(in), optional :: current
    integer, save :: attempts = 0
    if (.not. present(current)) error stop 'Missing baseline misfit'
    if (abs(current-4.0_dp) > 1.e-12_dp) error stop 'Incorrect baseline misfit'
    attempts = attempts + 1
    previous = current
    total = 5.0_dp
    if (attempts > 1) total = 3.0_dp
  end subroutine
end module gll_fixture
