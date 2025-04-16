program fwat_optimize
  use optimize_grid, only: OptGridFlow
  use fwat_mpi
  use config
  use input_params, fpar => fwat_par_global
  use argparse, only: parse_args_optimize
  use specfem_par, only: LOCAL_PATH
  use opt_io, only: write_model
  use model_grid_data
  use logger

  implicit none

  type(OptGridFlow) :: fop

  call init_mpi()
  call init_mpi_fwat()

  call parse_args_optimize()

  call fpar%read(FWAT_PAR_FILE)
  call read_parameter_file(.true.)

  call fop%init(.true.)

  if (fpar%update%OPT_METHOD == 1) then
    call fop%get_SD_direction()
  elseif (fpar%update%OPT_METHOD == 2) then
    call fop%get_lbfgs_direction()
  else
    call exit_MPI(0, 'Unknown optimization method')
  endif

  ! if (fpar%update%DO_LS) then
  !   call fop%run_linesearch()
  ! endif
  call fop%model_update()

  call write_grid_model(fop%model_fname, fop%model)
  call write_grid_model(fop%output_model_path, fop%model)

  call synchronize_all()
  call log%write('*******************************************', .false.)
  call log%write('************ OPTIMIZATION DONE ************', .false.)
  call log%write('*******************************************', .false.)
  call finalize_MPI()

end program fwat_optimize