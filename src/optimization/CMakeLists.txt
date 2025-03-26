# object libraries of preproc
file(GLOB SOURCES_Fortran "*.f90")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/fwat.*$")
add_library(opt_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(opt_objs PRIVATE shared_objs shared_fwat_objs
                      cuda_objs)

# executables
set(APP_SOURCES
"${PROJECT_SOURCE_DIR}/src/optimization/fwat_optimize.f90"
)
foreach(execsourcefile ${APP_SOURCES})
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)
  add_executable(${EXEC_NAME} ${execsourcefile}
    ${SOURCES_SUBS} 
    ${HEADERS}
  )

  set_target_properties(${EXEC_NAME} PROPERTIES PREFIX "x")
  target_link_libraries(${EXEC_NAME} PUBLIC shared_fwat_objs
                        shared_objs specfem3D_objs opt_objs
                        postproc_objs cuda_objs
  )
  if (USE_CUDA)
    target_link_libraries(${EXEC_NAME} PUBLIC cuda_kernels_objs)
  endif()
  target_link_libraries(${EXEC_NAME} PUBLIC ${MPI_LIBRARIES} fortran-yaml ${YAML_CPP_LIBRARIES})
  set_target_properties(${EXEC_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
endforeach()