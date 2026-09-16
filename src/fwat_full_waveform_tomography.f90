program fwat_full_waveform_tomography
  ! Single-data inversion on a fixed GLL mesh. Models and gradient histories
  ! stay in per-rank SPECFEM binary files throughout SD/L-BFGS optimization.
  use config
  use fwat_mpi
  use argparse, only: parse_invert_args
  use input_params, only: fpar => fwat_par_global
  use common_lib, only: get_dat_type, get_kernel_names, mkdir, cp
  use meshfem3D_subs, only: meshfem3D_fwat
  use generate_databases_subs, only: generate_databases_fwat
  use preproc_fwd, only: PrepareFWD
  use post_processing, only: PostFlow, remove_ekernel
  use optimize_gll, only: OptGLLFlow, write_gll_vector
  use shared_parameters, only: LOCAL_PATH, MODEL, IMODEL, ANISOTROPY, ANISOTROPIC_KL, &
    SIMULATION_TYPE, SAVE_FORWARD, SAVE_MESH_FILES, NPROC, ATTENUATION, ADIOS_ENABLED, HDF5_ENABLED, &
    COUPLE_WITH_INJECTION_TECHNIQUE, INJECTION_TECHNIQUE_TYPE, MESH_A_CHUNK_OF_THE_EARTH, &
    NUMBER_OF_SIMULTANEOUS_RUNS, HDF5_IO_NODES, TOMOGRAPHY_PATH
  use specfem_par, only: OUTPUT_FILES, ELASTIC_SIMULATION, ACOUSTIC_SIMULATION, POROELASTIC_SIMULATION
  use logger, only: log
  implicit none

  type(PrepareFWD) :: fwd
  type(PostFlow) :: post
  type(OptGLLFlow) :: opt
  integer :: first, iteration_count, iter, itype, i
  character(len=MAX_STRING_LEN) :: database_path, solver_output_path, msg
  real(kind=dp) :: current_misfit
  logical :: attenuation_model

  ! Initialize one MPI group before parsing, so --help needs no input files.
  call init_mpi_single_group()
  call init_mpi_fwat()
  call parse_invert_args(first, iteration_count)
  ! Select GLL mode before reading YAML; MODEL_GRID is then unnecessary.
  use_gll = .true.
  call fpar%read(FWAT_PAR_FILE)
  call read_parameter_file(.true.)
  attenuation_model = ATTENUATION

  ! Validate the selected data type, elastic parameterization, optimizer,
  ! and MPI/database layout before writing model or solver files.
  if (count(fpar%postproc%INV_TYPE) /= 1) &
    call exit_MPI(worldrank, 'xspecfwat requires exactly one POSTPROC.INV_TYPE')
  itype = 0
  do i = 1, NUM_INV_TYPE
    if (INV_TYPE_NAMES(i) == simu_type) itype = i
  enddo
  if (.not. fpar%postproc%INV_TYPE(itype)) &
    call exit_MPI(worldrank, '-s must match the enabled POSTPROC.INV_TYPE')
  if (fpar%update%MODEL_TYPE /= 1) call exit_MPI(worldrank, 'GLL inversion currently supports MODEL_TYPE: 1 (vp/vs/rho)')
  if (fpar%update%OPT_METHOD /= 1 .and. fpar%update%OPT_METHOD /= 2) &
    call exit_MPI(worldrank, 'GLL inversion supports OPT_METHOD: 1 (SD) or 2 (L-BFGS)')
  if (fpar%update%ITER_START < 0 .or. first < fpar%update%ITER_START) &
    call exit_MPI(worldrank, 'Starting model must be at or after non-negative ITER_START')
  if (fpar%update%LBFGS_M_STORE < 1 .or. fpar%update%MAX_SLEN <= 0.0_cr) &
    call exit_MPI(worldrank, 'LBFGS_M_STORE and MAX_SLEN must be positive')
  if (fpar%update%VPVS_RATIO_RANGE(1) <= sqrt(4.0_cr/3.0_cr) .or. &
      fpar%update%VPVS_RATIO_RANGE(2) < fpar%update%VPVS_RATIO_RANGE(1)) &
    call exit_MPI(worldrank, 'Invalid VPVS_RATIO_RANGE for an elastic GLL model')
  if (fpar%update%DO_LS) then
    if (fpar%update%MAX_SUB_ITER < 1 .or. fpar%update%MAX_SHRINK <= 0.0_cr .or. &
        fpar%update%MAX_SHRINK >= 1.0_cr .or. fpar%update%C1 <= 0.0_cr .or. fpar%update%C1 >= 1.0_cr) &
      call exit_MPI(worldrank, 'Invalid GLL line-search parameters')
  endif
  if (NPROC /= worldsize .or. NUMBER_OF_SIMULTANEOUS_RUNS /= 1 .or. HDF5_IO_NODES /= 0) &
    call exit_MPI(worldrank, 'GLL inversion requires one MPI group with exactly NPROC ranks and no separate I/O ranks')
  if (ADIOS_ENABLED .or. HDF5_ENABLED .or. MESH_A_CHUNK_OF_THE_EARTH) &
    call exit_MPI(worldrank, 'GLL inversion requires binary databases and the standard internal mesher')

  ! Event simulations change the global paths; retain the base directories
  ! so later stages and iterations can restore them.
  database_path = LOCAL_PATH
  solver_output_path = OUTPUT_FILES
  local_path_backup = database_path
  call fpar%select_simu_type()
  if (fpar%sim%SIGMA_H <= 0.0_cr .or. fpar%sim%SIGMA_V <= 0.0_cr) &
    call exit_MPI(worldrank, 'GLL PDE smoothing requires positive SIGMA_H and SIGMA_V')
  if (min(fpar%postproc%TAPER_H_SUPPRESS, fpar%postproc%TAPER_H_BUFFER, &
          fpar%postproc%TAPER_V_SUPPRESS, fpar%postproc%TAPER_V_BUFFER) < 0.0_cr) &
    call exit_MPI(worldrank, 'Taper distances must be non-negative')
  if (fpar%sim%PRECOND_TYPE < 1 .or. fpar%sim%PRECOND_TYPE > 3) &
    call exit_MPI(worldrank, 'GLL inversion requires PRECOND_TYPE: 1, 2 or 3')
  call get_kernel_names()

  call mkdir(database_path)
  call mkdir(solver_output_path)
  call synchronize_all()

  ! Generate the mesh once. Keeping its partition and element order fixed
  ! lets every L-BFGS history vector refer to the same GLL points.
  if (first == 0) then
    if (IMODEL /= IMODEL_USER_EXTERNAL .and. IMODEL /= IMODEL_GLL) &
      call exit_MPI(worldrank, 'Initial MODEL in Par_file must be external or gll')
    if (IMODEL == IMODEL_USER_EXTERNAL .and. len_trim(fpar%update%INIT_MODEL_PATH) == 0) &
      call exit_MPI(worldrank, 'INIT_MODEL_PATH is required for an external initial model')
    if (IMODEL == IMODEL_USER_EXTERNAL) &
      call cp(fpar%update%INIT_MODEL_PATH, trim(TOMOGRAPHY_PATH) // '/tomography_model.h5')
  else
    MODEL = 'gll'
    IMODEL = IMODEL_GLL
  endif
  ANISOTROPY = .false.
  ANISOTROPIC_KL = .false.
  if (simu_type == SIMU_TYPE_TELE) then
    COUPLE_WITH_INJECTION_TECHNIQUE = .true.
    INJECTION_TECHNIQUE_TYPE = 3
  endif
  SIMULATION_TYPE = 1
  SAVE_FORWARD = .true.
  SAVE_MESH_FILES = .false.
  call meshfem3D_fwat(fpar%sim%mesh_par_file)

  do iter = first, first+iteration_count-1
    write(model_name, '(A,I2.2)') 'M', iter
    ! 1. Build from the chosen external/GLL initial model, then from the
    ! updated GLL files in LOCAL_PATH on every subsequent iteration.
    LOCAL_PATH = database_path
    OUTPUT_FILES = solver_output_path
    local_path_backup = database_path
    local_path_fwat = database_path
    SIMULATION_TYPE = 1
    SAVE_FORWARD = .true.
    call generate_databases_fwat()

    ! 2. Run forward/measurement/adjoint simulations for every selected event.
    run_mode = FORWARD_ADJOINT
    call get_dat_type()
    call fpar%acqui%read()
    if (fpar%acqui%nevents < 1) call exit_MPI(worldrank, 'No events in the selected source list')
    call fwd%init()
    if (.not. ELASTIC_SIMULATION .or. ACOUSTIC_SIMULATION .or. POROELASTIC_SIMULATION) &
      call exit_MPI(worldrank, 'GLL inversion currently requires a purely elastic mesh')
    ! Keep the objective in double precision for line search; reading the
    ! text misfit files would introduce rounding into its acceptance test.
    current_misfit = 0.0_dp
    do i = 1, fpar%acqui%nevents
      fwd%ievt = i
      call fwd%prepare_for_event()
      call fwd%simulation()
      current_misfit = current_misfit + fwd%obj_func
    enddo
    call fwd%destroy()
    call fpar%acqui%finalize()
    call log%finalize()

    ! 3. Sum event kernels, smooth and taper on GLL points, then save the
    ! processed gradient used by the optimizer and its future history pairs.
    LOCAL_PATH = database_path
    OUTPUT_FILES = solver_output_path
    run_mode = 1  ! Sum event kernels rather than reload a previous sum.
    SIMULATION_TYPE = 1
    call post%init()
    call post%init_for_type(itype)
    call post%sum_kernel()
    ! Apply preconditioning to the kernels here, or save the inverse Hessian
    ! for application when the optimizer constructs its search direction.
    if (fpar%postproc%IS_HESS_PRECOND) then
      call post%apply_precond()
    else
      call post%sum_precond()
    endif
    call post%pde_smooth()
    call post%taper_kernel_gll()
    call post%write_gradient_gll()
    call remove_ekernel()
    call fpar%acqui%finalize()
    call post%finalize()
    call log%finalize()

    ! 4. Archive this iteration's GLL model when optimization starts,
    ! read the history, construct an SD/L-BFGS direction, and choose
    ! either a fixed step or an Armijo step using new forward evaluations.
    call log%init('output_optimize_gll_'//trim(model_name)//'.log')
    ! Adjoint simulations temporarily disable attenuation; retain the model
    ! setting when saving quality factors and building later GLL databases.
    ATTENUATION = attenuation_model
    call opt%init(iter)
    ! Line searches and later iterations must use the updated LOCAL_PATH model.
    MODEL = 'gll'
    IMODEL = IMODEL_GLL
    call opt%get_direction()
    if (opt%converged) then
      call log%write('GLL gradient is zero; inversion finished')
      if (worldrank == 0) print *, 'GLL gradient is zero; inversion finished at ', trim(model_name)
      exit
    endif
    if (is_output_direction) &
      call write_gll_vector(trim(OPT_DIR)//'/DIRECTION_'//trim(model_name), parameter_names, opt%direction)
    step_len = fpar%update%MAX_SLEN
    if (fpar%update%DO_LS) then
      call opt%run_linesearch(current_misfit)
    else
      call opt%update_trial()
    endif
    ! Install the update for the next database build. Its history copy is
    ! written only when the next iteration reaches optimizer initialization.
    call write_gll_vector(database_path, parameter_names, opt%trial)
    write(msg, '(A,ES12.4,A,I2.2,A)') &
      'Accepted step ', step_len, '; M', iter+1, ' saved to '//trim(database_path)
    call log%write(msg, .true.)
    if (worldrank == 0) print *, trim(msg)
    ! Release optimization arrays before allocating the next wave simulation.
    call opt%finalize()
    call log%finalize()
  enddo
  call log%finalize()
  call synchronize_all()
  call finalize_mpi()
end program fwat_full_waveform_tomography
