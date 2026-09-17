program fwat_optimize
  use optimize_grid, only: OptGridFlow
  use optimize_gll, only: OptGLLFlow, write_gll_vector, gll_model_path
  use fwat_mpi
  use config
  use input_params, fpar => fwat_par_global
  use argparse, only: parse_args_optimize
  use model_grid_data, only: write_grid_model
  use logger, only: log
  use param_check, only: check_model_name, check_inv_type, check_optimize_params, &
                        check_nevents, check_elastic_mesh, check_misfit

  implicit none

  integer :: iter

  call init_mpi()
  call init_mpi_fwat()

  call parse_args_optimize()

  call fpar%read(FWAT_PAR_FILE)
  call read_parameter_file(.true.)

  ! Reject an inconsistent setup before touching any model or history file.
  call check_model_name(iter)

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

    call check_optimize_params(iter)

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
    use common_lib, only: get_kernel_names, get_dat_type
    use kernel_io, only: read_mesh_databases_for_init
    use misfit_mod, only: read_evt_misfit
    use specfem_par, only: LOCAL_PATH, MODEL, IMODEL, SIMULATION_TYPE, SAVE_FORWARD
    type(OptGLLFlow) :: opt
    integer :: itype, ievt
    real(kind=dp) :: current_misfit
    character(len=MAX_STRING_LEN) :: database_path, msg

    ! A GLL update handles one data type, so the enabled one selects it.
    call check_inv_type(itype)
    simu_type = INV_TYPE_NAMES(itype)
    local_path_backup = LOCAL_PATH
    call fpar%select_simu_type()
    call check_optimize_params(iter)
    database_path = local_path_fwat
    call get_kernel_names()
    call log%init('output_optimize_gll_'//trim(model_name)//'.log')

    ! Load the existing partition and its unscaled material arrays, without remeshing.
    SIMULATION_TYPE = 1
    SAVE_FORWARD = .true.
    call read_mesh_databases_for_init()
    call check_elastic_mesh()
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
      call check_nevents()
      current_misfit = 0.0_dp
      do ievt = 1, fpar%acqui%nevents
        current_misfit = current_misfit + read_evt_misfit(model_name, ievt)
      enddo
      call fpar%acqui%finalize()
      call check_misfit(current_misfit)
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
