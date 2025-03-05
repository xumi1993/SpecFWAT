# objs
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(GLOB SOURCES_Fortran "*.f90")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
add_library(measadj_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(measadj_objs PRIVATE MPI::MPI_Fortran)
target_link_libraries(measadj_objs PRIVATE fwat_shared_objs)

