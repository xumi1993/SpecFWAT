# object libraries of preproc
file(GLOB SOURCES_Fortran "*.f90")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/fwat.*$")
add_library(mesh_databases_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(mesh_databases_objs PRIVATE shared_objs
                      cuda_objs generate_databases_objs
                      meshfem3D_objs shared_fwat_objs
)

# executables
set(APP_SOURCES
    "${PROJECT_SOURCE_DIR}/src/mesh_database_api/fwat_mesh_databases.f90"
)
foreach(execsourcefile ${APP_SOURCES})
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)
  add_executable(${EXEC_NAME} ${execsourcefile}
    ${SOURCES_SUBS} 
    ${HEADERS}
  )

  set_target_properties(${EXEC_NAME} PROPERTIES PREFIX "x")
  target_link_libraries(${EXEC_NAME} PUBLIC
                        shared_objs mesh_databases_objs
                        generate_databases_objs cuda_objs
                        meshfem3D_objs external_values_objs
                        utils_objs hdf5_interface_objs config_objs
                        shared_fwat_objs specfem3D_objs fortran-yaml
                        ${MPI_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
  )

  if (USE_CUDA)
    target_link_libraries(${EXEC_NAME} PUBLIC cuda_kernels_objs)
  endif()
  set_target_properties(${EXEC_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
endforeach()