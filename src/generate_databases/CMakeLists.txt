# objs
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(GLOB SOURCES_Fortran "*.[fF]90")
file(GLOB SOURCES_C "*.c")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
find_package(MPI REQUIRED)
add_library(generate_databases_objs OBJECT ${SOURCES_Fortran} ${SOURCES_C})
target_link_libraries(generate_databases_objs PRIVATE shared_objs)
target_link_libraries(generate_databases_objs PRIVATE MPI::MPI_Fortran)

