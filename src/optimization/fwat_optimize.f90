program fwat_optimize
  use optimize_grid, only: OptGridFlow
  use optimize_gll, only: OptGLLFlow, write_gll_vector, gll_model_path
  use fwat_mpi
  use config
  use input_params, fpar => fwat_par_global
  use argparse, only: parse_args_optimize
  use model_grid_data, only: write_grid_model
  use logger, only: log

  implicit none

  call init_mpi()
  call init_mpi_fwat()

  call parse_args_optimize()

  call fpar%read(FWAT_PAR_FILE)
  call read_parameter_file(.true.)

  if (use_gll) then
    call update_gll_model()
  else
    call update_grid_model()
  endif

  call synchronize_all()
  call log%write('*******************************************', .false.)
  call log%write('************ OPTIMIZATION DONE ************', .false.)
  call log%write('*******************************************', .false.)
  call log%finalize()

  call synchronize_all()
  call finalize_MPI()

contains

  subroutine update_grid_model()
    type(OptGridFlow) :: fop

    call fop%init()
    if (fpar%update%OPT_METHOD == 1) then
      call fop%get_SD_direction()
    elseif (fpar%update%OPT_METHOD == 2) then
      call fop%get_lbfgs_direction()
      if (fop%angle > 90) then
        call log%write('Stop optimization here...')
        call finalize_MPI()
        stop
      endif
    elseif (fpar%update%OPT_METHOD == 3) then
      call fop%get_CG_direction()
    else
      call exit_MPI(0, 'Unknown optimization method')
    endif

    if (fpar%update%DO_LS) call fop%run_linesearch()
    call fop%model_update()
    call write_grid_model(fop%model_fname, fop%model)
    call write_grid_model(fop%output_model_path, fop%model)
  end subroutine update_grid_model

  subroutine update_gll_model()
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use common_lib, only: get_kernel_names, get_dat_type
    use kernel_io, only: read_mesh_databases_for_init
    use misfit_mod, only: read_evt_misfit
    use specfem_par, only: LOCAL_PATH, MODEL, IMODEL, SIMULATION_TYPE, SAVE_FORWARD, &
      ANISOTROPY, ANISOTROPIC_KL, NPROC, NUMBER_OF_SIMULTANEOUS_RUNS, HDF5_IO_NODES, &
      ADIOS_ENABLED, HDF5_ENABLED, MESH_A_CHUNK_OF_THE_EARTH, ELASTIC_SIMULATION, &
      ACOUSTIC_SIMULATION, POROELASTIC_SIMULATION
    type(OptGLLFlow) :: opt
    integer :: iter, ios, itype, ievt
    real(kind=dp) :: current_misfit
    character(len=MAX_STRING_LEN) :: database_path, msg

    if (len_trim(model_name) /= 3 .or. model_name(1:1) /= 'M' .or. &
        verify(model_name(2:3), '0123456789') /= 0) &
      call exit_MPI(worldrank, 'GLL model must have the form M00 through M98')
    read(model_name(2:3), '(I2)', iostat=ios) iter
    if (ios /= 0) call exit_MPI(worldrank, 'Invalid GLL model index')
    if (iter > 98 .or. fpar%update%ITER_START < 0 .or. iter < fpar%update%ITER_START) &
      call exit_MPI(worldrank, 'GLL model must be between ITER_START and M98')
    if (count(fpar%postproc%INV_TYPE) /= 1) &
      call exit_MPI(worldrank, 'GLL optimization requires exactly one POSTPROC.INV_TYPE')
    if (fpar%update%MODEL_TYPE /= 1) &
      call exit_MPI(worldrank, 'GLL optimization supports MODEL_TYPE: 1 (vp/vs/rho)')
    if (fpar%update%OPT_METHOD /= 1 .and. fpar%update%OPT_METHOD /= 2) &
      call exit_MPI(worldrank, 'GLL optimization supports OPT_METHOD: 1 (SD) or 2 (L-BFGS)')
    if (fpar%update%LBFGS_M_STORE < 1 .or. fpar%update%MAX_SLEN <= 0.0_cr) &
      call exit_MPI(worldrank, 'LBFGS_M_STORE and MAX_SLEN must be positive')
    if (fpar%update%VPVS_RATIO_RANGE(1) <= sqrt(4.0_cr/3.0_cr) .or. &
        fpar%update%VPVS_RATIO_RANGE(2) < fpar%update%VPVS_RATIO_RANGE(1)) &
      call exit_MPI(worldrank, 'Invalid VPVS_RATIO_RANGE for an elastic GLL model')
    if (NPROC /= worldsize .or. NUMBER_OF_SIMULTANEOUS_RUNS /= 1 .or. HDF5_IO_NODES /= 0) &
      call exit_MPI(worldrank, 'GLL optimization requires one MPI group with exactly NPROC ranks and no I/O ranks')
    if (ADIOS_ENABLED .or. HDF5_ENABLED .or. MESH_A_CHUNK_OF_THE_EARTH) &
      call exit_MPI(worldrank, 'GLL optimization requires binary databases and the standard internal mesher')
    if (IMODEL /= IMODEL_USER_EXTERNAL .and. IMODEL /= IMODEL_GLL) &
      call exit_MPI(worldrank, 'GLL optimization requires MODEL = external or gll')
    if (iter > 0 .and. IMODEL /= IMODEL_GLL) &
      call exit_MPI(worldrank, 'Set MODEL = gll for GLL updates after M00')
    if (ANISOTROPY .or. ANISOTROPIC_KL) &
      call exit_MPI(worldrank, 'GLL optimization requires isotropic model databases')
    if (fpar%update%DO_LS) then
      if (fpar%update%MAX_SUB_ITER < 1 .or. fpar%update%MAX_SHRINK <= 0.0_cr .or. &
          fpar%update%MAX_SHRINK >= 1.0_cr .or. fpar%update%C1 <= 0.0_cr .or. fpar%update%C1 >= 1.0_cr) &
        call exit_MPI(worldrank, 'Invalid GLL line-search parameters')
    endif

    do itype = 1, NUM_INV_TYPE
      if (fpar%postproc%INV_TYPE(itype)) simu_type = INV_TYPE_NAMES(itype)
    enddo
    local_path_backup = LOCAL_PATH
    call fpar%select_simu_type()
    database_path = local_path_fwat
    call get_kernel_names()
    call log%init('output_optimize_gll_'//trim(model_name)//'.log')

    ! Load the existing partition and its unscaled material arrays, without remeshing.
    SIMULATION_TYPE = 1
    SAVE_FORWARD = .true.
    call read_mesh_databases_for_init()
    if (.not. ELASTIC_SIMULATION .or. ACOUSTIC_SIMULATION .or. POROELASTIC_SIMULATION) &
      call exit_MPI(worldrank, 'GLL optimization requires a purely elastic mesh')
    call opt%init(iter)
    call opt%get_direction()
    if (opt%converged) then
      call log%write('GLL gradient is zero; no model update', .true.)
      call opt%finalize()
      return
    endif
    if (is_output_direction) &
      call write_gll_vector(trim(OPT_DIR)//'/DIRECTION_'//trim(model_name), parameter_names, opt%direction)

    step_len = fpar%update%MAX_SLEN
    MODEL = 'gll'
    IMODEL = IMODEL_GLL
    if (fpar%update%DO_LS) then
      call get_dat_type()
      call fpar%acqui%read()
      if (fpar%acqui%nevents < 1) call exit_MPI(worldrank, 'No events in the selected source list')
      current_misfit = 0.0_dp
      do ievt = 1, fpar%acqui%nevents
        current_misfit = current_misfit + read_evt_misfit(model_name, ievt)
      enddo
      call fpar%acqui%finalize()
      if (.not. ieee_is_finite(current_misfit)) call exit_MPI(worldrank, 'Non-finite current GLL misfit')
      call opt%run_linesearch(current_misfit)
    else
      call opt%update_trial()
    endif

    call write_gll_vector(gll_model_path(iter+1), parameter_names, opt%trial)
    call write_gll_vector(database_path, parameter_names, opt%trial)
    write(msg, '(A,ES12.4,A,I2.2,A)') &
      'Accepted step ', step_len, '; M', iter+1, ' saved to '//trim(database_path)
    call log%write(msg, .true.)
    call opt%finalize()
  end subroutine update_gll_model

end program fwat_optimize
