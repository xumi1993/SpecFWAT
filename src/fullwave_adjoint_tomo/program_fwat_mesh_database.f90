program fwat_mesh_databse
  use fullwave_adjoint_tomo_par
  use fwat_input
  use fwat_utils
  use specfem_par, only: myrank, sizeprocs, LOCAL_PATH

  implicit none

  character(len=MAX_STRING_LEN) :: simutype, meshfile
  integer                       :: narg


  call init_mpi()
  call world_rank(myrank)
  call world_size(sizeprocs)

  call read_parameter_file(myrank,.true.)
  call read_fwat_par_file()

  if (myrank == 0) then
    call system('mkdir -p '//trim(OUTPUT_FILES_BASE)//'/')
    call system('mkdir -p '//trim(LOCAL_PATH)//'/')
  endif

  narg = command_argument_count()
  if (narg > 1) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xfwat_mesh_database simu_type'
      print *,'simu_type  --- simulation type: noise, tele, rf, telecd'
      stop 'Please check command line arguments'
    endif
  elseif (narg == 1) then
    call get_command_argument(1, simutype)
  else
    simutype = ''
  endif

  if (simutype == 'noise') then
    meshfile = noise_par%MESH_FILE_PATH
  elseif (index(simutype,'tele') /=0 .or. simutype == 'rf') then
    meshfile = tele_par%MESH_FILE_PATH
  elseif (simutype == 'leq') then
    meshfile = leq_par%MESH_FILE_PATH
  else
    meshfile = 'DATA/meshfem3D_files/Mesh_Par_file'
  endif
  call read_mesh_parameter_file_fwat(meshfile)
  call meshfem3d_fwat()
  call generate_database_fwat(USE_H5)

  ! MPI finish
  call finalize_mpi()
end program