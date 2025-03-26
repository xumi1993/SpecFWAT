program fwat_optimize
  use optimize
  use fwat_mpi
  use config
  use input_params, fpar => fwat_par_global
  use argparse, only: parse_args_optimize
  use specfem_par, only: LOCAL_PATH
  use opt_io, only: write_model

  implicit none

  type(OptFlow) :: fop

  call init_mpi()
  call init_mpi_fwat()

  call parse_args_optimize()

  call fpar%read(FWAT_PAR_FILE)

  call fop%init(.true.)

  if (fpar%update%OPT_METHOD == 1) then
    call fop%get_SD_direction()
  elseif (fpar%update%OPT_METHOD == 2) then
    call fop%get_lbfgs_direction()
  else
    call exit_MPI(0, 'Unknown optimization method')
  endif

  if (fpar%update%DO_LS) then
    continue
  else
    call fop%model_update()
  endif

  call write_model(LOCAL_PATH, fop%model)
  call write_model(trim(OPT_DIR)//'/'//trim(fop%model_next), fop%model)

  
  call synchronize_all()
  call finalize_MPI()

end program fwat_optimize