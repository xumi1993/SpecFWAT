program fwat_mesh_databases
  use generate_databases_subs
  use meshfem3D_subs
  use fwat_mpi
  use config
  use input_params, fpar => fwat_par_global
  use argparse, only: parse_args_mesh_databases
  use shared_parameters, only : ANISOTROPY
  use common_lib, only: mkdir

  ! init MPI
  call init_mpi()
  call init_mpi_fwat()

  ! parse command line arguments
  call parse_args_mesh_databases()

  ! read input parameters
  call fpar%read(FWAT_PAR_FILE)
  call read_parameter_file(.true.)
  local_path_backup = trim(LOCAL_PATH)

  if ( fpar%update%MODEL_TYPE > 1) then
    ANISOTROPY = .true.
    ANISOTROPIC_KL = .true. 
  else
    ANISOTROPY = .false.
    ANISOTROPIC_KL = .false.
  end if

  ! PML parameters
  if (PML_CONDITIONS) then
    SIMULATION_TYPE = 1
    SAVE_FORWARD = .true.
  end if

  ! select simulation type
  call fpar%select_simu_type()

  ! set local path
  LOCAL_PATH = local_path_fwat
  call mkdir(LOCAL_PATH)
  call mkdir(OUTPUT_PATH)
  call synchronize_all()

  ! generate mesh and databases
  call meshfem3D_fwat(fpar%sim%mesh_par_file)
  call generate_databases_fwat()

  call finalize_mpi()
end program fwat_mesh_databases