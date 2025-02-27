include_directories(
  ${PROJECT_SOURCE_DIR}/include
)

# object libraries of specfem_interface
file(GLOB SOURCES_SPECFEM_API "specfem_interface/*.f90")
add_library(specfem_api_objs OBJECT ${SOURCES_SPECFEM_API})
target_link_libraries(specfem_api_objs PRIVATE specfem3D_objs)

# object libraries of preproc
file(GLOB SOURCES_Fortran "*.f90")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/fwat.*$")
add_library(preproc_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(preproc_objs PRIVATE shared_objs shared_fwat_objs cuda_objs)

# executables
set(APP_SOURCES
    "${PROJECT_SOURCE_DIR}/src/preproc/fwat_fwd_measure_adj.f90"
)
foreach(execsourcefile ${APP_SOURCES})
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)
  add_executable(${EXEC_NAME} ${execsourcefile}
    ${SOURCES_SUBS} 
    ${HEADERS}
  )

  set_target_properties(${EXEC_NAME} PROPERTIES PREFIX "x")
  target_link_libraries(${EXEC_NAME} PUBLIC preproc_objs shared_fwat_objs shared_objs cuda_objs specfem3D_objs specfem_api_objs)
  if (USE_CUDA)
    target_link_libraries(${EXEC_NAME} PUBLIC cuda_kernels_objs)
  endif()
  target_link_libraries(${EXEC_NAME} PUBLIC ${MPI_LIBRARIES} fortran-yaml ${YAML_CPP_LIBRARIES})
  set_target_properties(${EXEC_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
    
endforeach()
