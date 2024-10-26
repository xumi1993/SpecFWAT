program fwat_mesh_databse
  use fullwave_adjoint_tomo_par
  use fwat_input
  use fwat_utils
  use specfem_par, only: myrank, sizeprocs

  implicit none

  call init_mpi()
  call world_rank(myrank)
  call world_size(sizeprocs)

  ! BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,.true.)
  call read_mesh_parameter_file_fwat(get_mesh_file_path(0))
  call read_fwat_par_file()
  call synchronize_all()

  call meshfem3d_fwat()

  call generate_database_fwat(.true.)
  call generate_database_fwat(.false.)

  ! MPI finish
  call finalize_mpi()
end program