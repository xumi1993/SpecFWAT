# object libraries of preproc
file(GLOB SOURCES_Fortran "*.f90")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/fwat.*$")
if (USE_CUDA)
  file (GLOB SOURCES_CUDA "*.cu")
else()
  file (GLOB SOURCES_CUDA "*_stubs.cpp")
endif()
add_library(postproc_objs OBJECT ${SOURCES_Fortran} ${SOURCES_CUDA})
target_link_libraries(postproc_objs PRIVATE shared_objs shared_fwat_objs ${MPI_LIBRARIES}
                      cuda_objs grid_objs mesh_subs_objs databases_subs_objs)

# executables
set(APP_SOURCES
"${PROJECT_SOURCE_DIR}/src/postproc/fwat_post_proc.f90"
)
foreach(execsourcefile ${APP_SOURCES})
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)
  add_executable(${EXEC_NAME} ${execsourcefile}
    ${SOURCES_SUBS} 
    ${HEADERS}
  )

  set_target_properties(${EXEC_NAME} PROPERTIES PREFIX "x")
  target_link_libraries(${EXEC_NAME} PUBLIC shared_fwat_objs
                        shared_objs specfem3D_objs postproc_objs
                        cuda_objs grid_objs external_values_objs
                        config_objs hdf5_interface_objs utils_objs
                        meshfem3D_objs generate_databases_objs
                        mesh_subs_objs databases_subs_objs
                        specfem_api_objs window_chi_objs
  )
  if (USE_CUDA)
    target_link_libraries(${EXEC_NAME} PUBLIC cuda_kernels_objs)
  endif()
  target_link_libraries(${EXEC_NAME} PUBLIC ${MPI_LIBRARIES} fortran-yaml ${YAML_CPP_LIBRARIES})
  set_target_properties(${EXEC_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
endforeach()