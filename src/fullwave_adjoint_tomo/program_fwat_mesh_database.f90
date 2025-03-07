program fwat_mesh_databse
  use fullwave_adjoint_tomo_par

  implicit none
  integer :: myrank, sizeprocs

  call init_mpi()
  call world_rank(myrank)
  call world_size(sizeprocs)

  ! BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,.true.)
  call read_mesh_parameter_file()
  call synchronize_all()

  call meshfem3d_fwat()

  call generate_database_fwat()



  ! MPI finish
  call finalize_mpi()
end program