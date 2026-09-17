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
    SIMULATION_TYPE, SAVE_FORWARD, SAVE_MESH_FILES, ATTENUATION, TOMOGRAPHY_PATH, &
    COUPLE_WITH_INJECTION_TECHNIQUE, INJECTION_TECHNIQUE_TYPE
  use specfem_par, only: OUTPUT_FILES
  use logger, only: log
  use param_check, only: check_inversion_params, check_nevents, check_elastic_mesh
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

  ! Event simulations change the global paths; retain the base directories
  ! so later stages and iterations can restore them.
  database_path = LOCAL_PATH
  solver_output_path = OUTPUT_FILES
  local_path_backup = database_path
  call fpar%select_simu_type()
  ! Resuming rebuilds the model from the GLL binaries of the previous
  ! iteration, whatever initial model the Par_file asks for, and vp/vs/rho
  ! updates always run on isotropic databases.
  if (first > 0) then
    MODEL = 'gll'
    IMODEL = IMODEL_GLL
  endif
  ANISOTROPY = .false.
  ANISOTROPIC_KL = .false.
  ! One driver process runs every step, so validate all of them before
  ! writing any model, database or solver file.
  call check_inversion_params(first, iteration_count, itype)
  call get_kernel_names()

  call mkdir(database_path)
  call mkdir(solver_output_path)
  call synchronize_all()

  ! Stage the external starting model where generate_databases expects it.
  if (first == 0 .and. IMODEL == IMODEL_USER_EXTERNAL) &
    call cp(fpar%update%INIT_MODEL_PATH, trim(TOMOGRAPHY_PATH) // '/tomography_model.h5')
  if (simu_type == SIMU_TYPE_TELE) then
    COUPLE_WITH_INJECTION_TECHNIQUE = .true.
    INJECTION_TECHNIQUE_TYPE = 3
  endif
  SIMULATION_TYPE = 1
  SAVE_FORWARD = .true.
  SAVE_MESH_FILES = .false.
  ! Generate the mesh once. Keeping its partition and element order fixed
  ! lets every L-BFGS history vector refer to the same GLL points.
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
    call check_nevents()
    call fwd%init()
    call check_elastic_mesh()
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
    call log%write('*******************************************', .false.)
    call log%write('*********** PRE-PROCESSING DONE ***********', .false.)
    call log%write('*******************************************', .false.)
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

    call log%write('*******************************************', .false.)
    call log%write('********** POST-PROCESSING DONE ***********', .false.)
    call log%write('*******************************************', .false.)
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
    call log%write('*******************************************', .false.)
    call log%write('*********** OPTIMIZATION DONE *************', .false.)
    call log%write('*******************************************', .false.)
    ! Release optimization arrays before allocating the next wave simulation.
    call opt%finalize()
    call log%finalize()
  enddo
  call log%finalize()
  call synchronize_all()
  call finalize_mpi()
end program fwat_full_waveform_tomography
