program fwat_mesh_databases
  use generate_databases_subs
  use meshfem3D_subs
  use fwat_mpi
  use config
  use input_params, fpar => fwat_par_global
  use argparse, only: parse_args_mesh_databases
  use shared_parameters, only : ANISOTROPY

  ! init MPI
  call init_mpi()
  call init_mpi_fwat()

  ! parse command line arguments
  call parse_args_mesh_databases()

  ! read input parameters
  call fpar%read(FWAT_PAR_FILE)
  call read_parameter_file(.true.)
  if ( fpar%update%MODEL_TYPE > 1) then
    ANISOTROPY = .true.
  else
    ANISOTROPY = .false.
  end if

  ! select simulation type
  call fpar%select_simu_type()

  ! generate mesh and databases
  call meshfem3D_fwat(fpar%sim%mesh_par_file)
  call generate_databases_fwat(.true.)

  call finalize_mpi()
end program fwat_mesh_databases