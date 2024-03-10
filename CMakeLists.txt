cmake_minimum_required(VERSION 3.8)

# set the project name
set(CMAKE_PROJECT_NAME "specfwat")
project(${CMAKE_PROJECT_NAME} VERSION 1.3.0 LANGUAGES C Fortran)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)


if ("${CMAKE_FORTRAN_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-std=gnu -fimplicit-none -frange-check -fmax-errors=10 -pedantic -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wunderflow -ffpe-trap=invalid,zero,overflow -O3 -finline-functions -Wno-argument-mismatch")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -ggdb -fbacktrace -fbounds-check")
elseif("${CMAKE_FORTRAN_COMPILER_ID}" MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-std03 -O3 -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -diag-disable 6477 -implicitnone -gen-interfaces -warn all -warn nointerfaces")
  set(CMAKE_Fortran_FLAGS_DEBUG "-check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv")
elseif("${CMAKE_FORTRAN_COMPILER_ID}" MATCHES "Cray")
  set(CMAKE_Fortran_FLAGS "-O3 -Onoaggress -Oipa0 -hfp2 -Ovector3 -Oscalar3 -Ocache2 -Ounroll2 -Ofusion2 -M 1193 -M 1438")
  set(CMAKE_Fortran_FLAGS_DEBUG "eC -eD -ec -en -eI -ea -g -G0")
endif()

find_package(LAPACK REQUIRED)
message(STATUS "LAPACK_FOUND: ${LAPACK_FOUND}")
message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")

find_package(BLAS REQUIRED)
message(STATUS "BLAS_FOUND: ${BLAS_FOUND}")
message(STATUS "BLAS_LIBRARIES: ${BLAS_LIBRARIES}")

# set(BUILD_TESTING OFF)
# setup mpi
find_package(MPI REQUIRED)
message(STATUS "MPI_FOUND: ${MPI_FOUND}")
message(STATUS "MPI_VERSION: ${MPI_VERSION}")

# setup source files
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
set(SOURCE_FWAT_DIR "${PROJECT_SOURCE_DIR}/src/fullwave_adjoint_tomo")

add_subdirectory(${SOURCE_DIR}/shared)
add_subdirectory(${SOURCE_DIR}/specfem3D)
add_subdirectory(${SOURCE_DIR}/meshfem3D)
add_subdirectory(${SOURCE_DIR}/generate_databases)
add_subdirectory(${SOURCE_DIR}/fullwave_adjoint_tomo)

# cuda sources
file(GLOB SOURCES_CUDA "${SOURCE_DIR}/cuda/*.c*")

# setup include directory
include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${SOURCE_FWAT_DIR}/sac_inc
  ${CMAKE_BINARY_DIR}/modules
  # ${CMAKE_CURRENT_BINARY_DIR}/src/fortran-yaml
  # ${CMAKE_CURRENT_BINARY_DIR}/src/fortran-csv
)

set(APP_SOURCES
  src/meshfem3D/meshfem3D.F95
  src/generate_databases/generate_databases.F95
  src/specfem3D/specfem3D.F95
  # ${SOURCE_FWAT_DIR}/program_fwat_mesh_database.f95
  ${SOURCE_FWAT_DIR}/program_fwat0_forward_data.f95
  )
foreach(execsourcefile ${APP_SOURCES})
  # get exec name
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)

  if(EXEC_NAME STREQUAL "program_fwat_mesh_database")
    set(EXEC_NAME "fwat_mesh_database")
  elseif(EXEC_NAME STREQUAL "program_fwat0_forward_data")
    set(EXEC_NAME "fwat0_forward_data")
  endif()

  # add executable
  add_executable(${EXEC_NAME} ${execsourcefile}
    ${SOURCES_SUBS} 
    ${HEADERS}
  )

  # set prefix before executable name
  set_target_properties(${EXEC_NAME} PROPERTIES PREFIX "x")

  # link objects
  target_link_libraries(${EXEC_NAME} PUBLIC "${EXEC_NAME}_objs")
  if (${EXEC_NAME} STREQUAL "fwat0_forward_data")
    target_link_libraries(${EXEC_NAME} PUBLIC measadj_objs flexwin_objs 
                          specfem_interface_objs read_mesh_interface_objs
                          fwat_shared_objs fwat_recfunc_objs fwat_telestf_objs
                          specfem3D_objs meshfem3D_objs
                          )
  endif()
  

  # link shared objects
  target_link_libraries(${EXEC_NAME} PUBLIC shared_objs)
  
  # link MPI
  target_link_libraries(${EXEC_NAME} PUBLIC MPI::MPI_Fortran)

  # link LAPACK
  target_link_libraries(${EXEC_NAME} PUBLIC LAPACK::LAPACK)

  # link BLAS
  target_link_libraries(${EXEC_NAME} PUBLIC BLAS::BLAS)
  
  # link libtau
  target_link_libraries(${EXEC_NAME} PUBLIC tau)

endforeach()

add_subdirectory(${SOURCE_DIR}/aux)

# if(BUILD_TESTING)
#   add_subdirectory(${PROJECT_SOURCE_DIR}/examples)
# endif()
