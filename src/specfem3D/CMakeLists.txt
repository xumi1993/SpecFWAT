# objs
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
file(GLOB SOURCES_Fortran "*.[fF]90")
file(GLOB SOURCES_C "*.c")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
set(FILTERED_SOURCE_FILES "")
foreach(FILE ${SOURCES_Fortran})
  if (NOT ${FILE} MATCHES "specfem3D.F90")
    list(APPEND FILTERED_SOURCE_FILES ${FILE})
  endif()
endforeach()
find_package(MPI REQUIRED)
add_library(specfem3D_objs OBJECT ${FILTERED_SOURCE_FILES} ${SOURCES_C})
target_link_libraries(specfem3D_objs PRIVATE shared_objs)
target_link_libraries(specfem3D_objs PRIVATE MPI::MPI_Fortran)

