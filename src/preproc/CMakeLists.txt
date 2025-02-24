include_directories(
  ${PROJECT_SOURCE_DIR}/include
)

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
  # target_link_libraries(${EXEC_NAME} PUBLIC "${EXEC_NAME}_objs" shared_objs cuda_objs ${MPI_LIBRARIES})
  target_link_libraries(${EXEC_NAME} PUBLIC shared_fwat_objs shared_objs cuda_objs)
  target_link_libraries(${EXEC_NAME} PUBLIC ${MPI_LIBRARIES} fortran-yaml)
  set_target_properties(${EXEC_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
    
endforeach()
