# objs
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(GLOB SOURCES_Fortran "*.f90")
file(GLOB SOURCES_C "*.c")
file(GLOB SOURCES_CUDA "${SOURCE_DIR}/cuda/*.c*")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
find_package(MPI REQUIRED)
message(STATUS ${SOURCES_FWAT_SHARED})
add_library(measadj_objs OBJECT ${SOURCES_Fortran} ${SOURCES_C} )
target_link_libraries(measadj_objs PRIVATE MPI::MPI_Fortran)
target_link_libraries(measadj_objs PRIVATE fwat_shared_objs fwat_recfunc_objs fwat_telestf_objs)

