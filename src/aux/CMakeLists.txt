# objs
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(GLOB SOURCES_Fortran "*.[fF]90")
file(GLOB SOURCES_C "*.c")

include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
find_package(MPI REQUIRED)

set(APP_SOURCES
  combine_vol_data_vtk.F95
  )
foreach(execsourcefile ${APP_SOURCES})
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)
  add_executable(${EXEC_NAME} ${execsourcefile} ${SOURCES_Fortran} ${SOURCES_C})
  set_target_properties( ${EXEC_NAME} PROPERTIES PREFIX "x")
  target_link_libraries(${EXEC_NAME} PRIVATE shared_objs)
  target_link_libraries(${EXEC_NAME} PRIVATE MPI::MPI_Fortran)

endforeach()


