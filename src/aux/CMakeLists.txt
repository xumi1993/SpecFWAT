# objs
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(GLOB SOURCES_C "*.c")

include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
find_package(MPI REQUIRED)

set(APP_SOURCES
  combine_vol_data_vtk.F95
  add_pert_gaus.f95
  add_pert_trigo.f95
  model_gll.f95
  )
foreach(execsourcefile ${APP_SOURCES})
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)
  if (${EXEC_NAME} STREQUAL "add_pert_trigo" )
    add_executable(${EXEC_NAME} ${execsourcefile} ${SOURCE_FWAT_DIR}/teleseis_stf/interpolation_mod.f90 npy.f90)
  elseif (${EXEC_NAME} STREQUAL "combine_vol_data_vtk")
    add_executable(${EXEC_NAME} ${execsourcefile} combine_vol_data_impl.f90 ${SOURCES_C})
  elseif(${EXEC_NAME} STREQUAL "model_gll")
    add_executable(${EXEC_NAME} ${execsourcefile} projection_on_FD_grid_mod_fwat.f90 npy.f90)
  elseif(${EXEC_NAME} STREQUAL "add_pert_gaus")
    add_executable(${EXEC_NAME} ${execsourcefile} npy.f90)
  endif()
  set_target_properties( ${EXEC_NAME} PROPERTIES PREFIX "x")
  target_link_libraries(${EXEC_NAME} PRIVATE shared_objs)
  target_link_libraries(${EXEC_NAME} PRIVATE MPI::MPI_Fortran)
  target_link_libraries(${EXEC_NAME} PRIVATE specfem3D_objs fwat_shared_objs)

endforeach()
