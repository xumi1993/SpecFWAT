# objs
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(GB_SOURCE_DIR ${SOURCE_DIR}/generate_databases)
file(GLOB SOURCES_Fortran "${GB_SOURCE_DIR}/*.[fF]90")
file(GLOB SOURCES_C "${GB_SOURCE_DIR}/*.c")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
set(FILTERED_SOURCE_FILES "")
foreach(FILE ${SOURCES_Fortran})
  if (NOT ${FILE} MATCHES "${GB_SOURCE_DIR}/generate_databases.f90")
    list(APPEND FILTERED_SOURCE_FILES ${FILE})
  endif()
endforeach()
find_package(MPI REQUIRED)
add_library(generate_databases_objs OBJECT ${FILTERED_SOURCE_FILES} ${SOURCES_C})
target_link_libraries(generate_databases_objs PRIVATE shared_objs)
target_link_libraries(generate_databases_objs PRIVATE MPI::MPI_Fortran)

