# objs
file(GLOB SOURCES_Fortran "*.f90")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
add_library(measadj_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(measadj_objs PRIVATE MPI::MPI_Fortran)
target_link_libraries(measadj_objs PRIVATE shared_fwat_objs)

