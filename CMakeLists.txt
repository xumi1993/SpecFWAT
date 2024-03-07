cmake_minimum_required(VERSION 3.8)

# set the project name
set(CMAKE_PROJECT_NAME "specfwat")
project(${CMAKE_PROJECT_NAME} VERSION 1.3.0 LANGUAGES C Fortran)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "clang")
# endif()


# if ("${CMAKE_FORTRAN_COMPILER_ID}" MATCHES "GNU")
set(MPI_Fortran_COMPILE_FLAGS "-std=gnu -fimplicit-none -frange-check -fmax-errors=10 -pedantic -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wunderflow -ffpe-trap=invalid,zero,overflow -O3 -finline-functions -Wno-argument-mismatch -finline-functions")
# endif()

set(BUILD_TESTING OFF)
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

# setup mpi
find_package(MPI REQUIRED)
message(STATUS "MPI_FOUND: ${MPI_FOUND}")
message(STATUS "MPI_VERSION: ${MPI_VERSION}")


# setup source files
# file(GLOB SOURCES "${SOURCE_DIR}/*.f90")
file(GLOB SOURCES_SHARED "${SOURCE_DIR}/shared/*.[fFc]*")
file(GLOB HEADERS "${PROJECT_SOURCE_DIR}/include/*.h")
# add sub modules
# add_subdirectory("${PROJECT_SOURCE_DIR}/external_libs")

# fortran-yaml-cpp
# add_subdirectory("${SOURCE_DIR}/fortran-yaml")

# fortran-csv
# add_subdirectory("${SOURCE_DIR}/fortran-csv")

# setup include directory
include_directories(
  ${PROJECT_SOURCE_DIR}/include
  # ${CMAKE_CURRENT_BINARY_DIR}/src/fortran-yaml
  # ${CMAKE_CURRENT_BINARY_DIR}/src/fortran-csv
)

set(APP_SOURCES
  src/meshfem3D/meshfem3D.f95
  src/generate_databases/generate_databases.f95
  )
foreach(execsourcefile ${APP_SOURCES})
  # get exec name
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)

  file(GLOB SOURCES_SUBS "${SOURCE_DIR}/${EXEC_NAME}/*.[fFc]*")

  # add executable
  add_executable(${EXEC_NAME} ${execsourcefile}
    ${SOURCES_SHARED} ${SOURCES_SUBS} 
    ${HEADERS}
  )
  target_compile_options(${EXEC_NAME} PUBLIC ${MPI_Fortran_COMPILE_FLAGS} "-Wl,-no_compact_unwind")
  # target_compile_options(${EXEC_NAME} PUBLIC "-O3")

  set_target_properties(
    ${EXEC_NAME}
    PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_${EXEC_NAME}"
    PREFIX "x"
)
  # link MPI
  target_link_libraries(${EXEC_NAME} PUBLIC MPI::MPI_Fortran)
  
  # # link fortran stdlib
  # target_link_libraries(${EXEC_NAME} PUBLIC fortran_stdlib)

  # # link yaml-app:
  # target_link_libraries(${EXEC_NAME} PUBLIC fortran-yaml)

  # # link fortran-csv
  # target_link_libraries(${EXEC_NAME} PUBLIC fortran-csv)
  
  # # link h5fortran
  # target_link_libraries(${EXEC_NAME} PUBLIC h5fortran)
endforeach()

# if(BUILD_TESTING)
#   add_subdirectory(${PROJECT_SOURCE_DIR}/examples)
# endif()