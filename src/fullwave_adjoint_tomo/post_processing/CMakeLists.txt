set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
find_package(MPI REQUIRED)


set(APP_SOURCES
  smooth_sem_pde.F90
  )

foreach(execsourcefile ${APP_SOURCES})
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)

  add_executable(${EXEC_NAME} ${execsourcefile})

  set_target_properties(${EXEC_NAME} PROPERTIES PREFIX "x")

  target_link_libraries(${EXEC_NAME} PUBLIC meshfem3D_objs generate_databases_objs specfem3D_objs shared_objs)

  target_link_libraries(${EXEC_NAME} PRIVATE MPI::MPI_Fortran)

  target_link_libraries(${EXEC_NAME} PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

endforeach()