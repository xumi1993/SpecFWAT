# objs
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
find_package(MPI REQUIRED)

file(GLOB SOURCES_Fortran "*.[Ff]90")
add_library(grid_interface_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(grid_interface_objs PRIVATE specfem3D_objs shared_objs MPI::MPI_Fortran)
target_link_libraries(grid_interface_objs PRIVATE ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
