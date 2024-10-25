set(FWAT_SPECFEM_INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/specfem_initialize_simulation_fwat.F90" 
                           "${CMAKE_CURRENT_SOURCE_DIR}/specfem_interface.F90"
                           "${CMAKE_CURRENT_SOURCE_DIR}/specfem_prepare_timerun_fwat.F90"
                           "${CMAKE_CURRENT_SOURCE_DIR}/specfem_save_adjoint_kernels_fwat.f90"
                           "${CMAKE_CURRENT_SOURCE_DIR}/specfem_setup_sources_receivers_fwat.f90"
                           )
set(FWAT_mesh_gen_INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/specfem_meshfem3D.f90"
                            "${CMAKE_CURRENT_SOURCE_DIR}/read_mesh_parameter_file_fwat.f90"
                            "${CMAKE_CURRENT_SOURCE_DIR}/create_regions_mesh_fwat.f90"
                            )
set(FWAT_read_mesh_INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/specfem_read_mesh_database.f90" )

include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
# Create an object library for the interface
add_library(specfem_interface_objs OBJECT ${FWAT_SPECFEM_INTERFACE} )
target_link_libraries(specfem_interface_objs PRIVATE fwat_shared_objs MPI::MPI_Fortran)

# Create an object library for the interface
add_library(mesh_gen_interface_objs OBJECT ${FWAT_mesh_gen_INTERFACE} )
target_link_libraries(mesh_gen_interface_objs PRIVATE meshfem3D_objs generate_databases_objs )

# Create an object library for the interface
add_library(read_mesh_interface_objs OBJECT ${FWAT_read_mesh_INTERFACE})
target_link_libraries(read_mesh_interface_objs PRIVATE generate_databases_objs )