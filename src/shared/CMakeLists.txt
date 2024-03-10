# objs
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(GLOB SOURCES_Fortran "*.f90")
file(GLOB SOURCES_C "*.c")
file(GLOB SOURCES_CUDA "${SOURCE_DIR}/cuda/*.c*")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
find_package(MPI REQUIRED)
add_library(shared_objs OBJECT ${SOURCES_Fortran} ${SOURCES_C} ${SOURCES_CUDA})
target_link_libraries(shared_objs PRIVATE MPI::MPI_Fortran)

