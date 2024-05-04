# objs
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(GLOB SOURCES_C "*.c")

include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
find_package(MPI REQUIRED)

add_library(aux_objs OBJECT
  "${CMAKE_CURRENT_SOURCE_DIR}/projection_on_FD_grid_mod_fwat.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/npy.f90"
)
target_link_libraries(aux_objs PRIVATE specfem3D_objs)

set(APP_SOURCES
  combine_vol_data_vtk.F90
  add_pert_gaus.f90
  add_pert_trigo.f90
  decompose_h5_gll.f90
  )
foreach(execsourcefile ${APP_SOURCES})
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)
  if (${EXEC_NAME} STREQUAL "add_pert_trigo" )
    add_executable(${EXEC_NAME} ${execsourcefile} ${SOURCE_FWAT_DIR}/teleseis_stf/interpolation_mod.f90)
  elseif (${EXEC_NAME} STREQUAL "combine_vol_data_vtk")
    add_executable(${EXEC_NAME} ${execsourcefile} combine_vol_data_impl.f90 ${SOURCES_C})
  elseif(${EXEC_NAME} STREQUAL "decompose_h5_gll")
    add_executable(${EXEC_NAME} ${execsourcefile})
  elseif(${EXEC_NAME} STREQUAL "add_pert_gaus")
    add_executable(${EXEC_NAME} ${execsourcefile} )
  endif()
  set_target_properties( ${EXEC_NAME} PROPERTIES PREFIX "x")
  target_link_libraries(${EXEC_NAME} PRIVATE shared_objs aux_objs)
  target_link_libraries(${EXEC_NAME} PRIVATE MPI::MPI_Fortran)
  target_link_libraries(${EXEC_NAME} PRIVATE specfem3D_objs fwat_shared_objs)
  target_link_libraries(${EXEC_NAME} PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

endforeach()
