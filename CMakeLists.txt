cmake_minimum_required(VERSION 3.8)

# set the project name
set(CMAKE_PROJECT_NAME "specfwat")
project(${CMAKE_PROJECT_NAME} VERSION 1.3.0 LANGUAGES C Fortran)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)


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
file(GLOB SOURCES_SHARED "${SOURCE_DIR}/shared/*.[fFc]*")
file(GLOB SOURCES_SPECFEM "${SOURCE_DIR}/specfem3D/*.[fF]90")
file(GLOB SOURCES_MESHFEM "${SOURCE_DIR}/meshfem3D/*.[fF]90")
file(GLOB SOURCES_GENERATE_DATABASES "${SOURCE_DIR}/generate_databases/*.[fF]90")
file(GLOB SOURCES_GENERATE_DATABASES_C "${SOURCE_DIR}/generate_databases/*.c")
list(APPEND SOURCES_GENERATE_DATABASES ${SOURCES_GENERATE_DATABASES_C})

# cuda sources
file(GLOB SOURCES_CUDA "${SOURCE_DIR}/cuda/*.c*")

# FWAT sources
set(SOURCES_FWAT_SHARED "${SOURCE_FWAT_DIR}/fullwave_adjoint_tomo_par.f90"
                       "${SOURCE_FWAT_DIR}/fwat_input.f90"
                       "${SOURCE_FWAT_DIR}/fwat_utils.f90")
file(GLOB SOURCES_RECFUNC "${SOURCE_FWAT_DIR}/preproc_recfunc/*.f90" )
file(GLOB SOURCES_FWAT_OPT "${SOURCE_FWAT_DIR}/optimization/*.[fF]90")
file(GLOB SOURCES_SPECFEMINTERFACE "${SOURCE_FWAT_DIR}/specfem_interface/*.[fF]*")


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
  src/meshfem3D/meshfem3D.F95
  src/generate_databases/generate_databases.F95
  src/specfem3D/specfem3D.F95
  src/aux/combine_vol_data_vtk.F95
  ${SOURCE_FWAT_DIR}/program_fwat_mesh_database.f95
  )
foreach(execsourcefile ${APP_SOURCES})
  # get exec name
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)

  file(GLOB SOURCES_SUBS "${SOURCE_DIR}/${EXEC_NAME}/*.[fFc]*")

  if (EXEC_NAME STREQUAL "specfem3D")
    list(APPEND SOURCES_SUBS ${SOURCES_CUDA})
  elseif(EXEC_NAME STREQUAL "combine_vol_data_vtk")
    file(GLOB SOURCES_AUX "${SOURCE_DIR}/aux/*.[Ff]90")
    list(APPEND SOURCES_SUBS ${SOURCES_AUX} "${SOURCE_DIR}/aux/vtk_writer.c")
  elseif(EXEC_NAME STREQUAL "program_fwat_mesh_database")
    set(EXEC_NAME "fwat_mesh_database")
    list(APPEND SOURCES_SUBS ${SOURCES_SPECFEM} ${SOURCES_MESHFEM}
         ${SOURCES_GENERATE_DATABASES} ${SOURCES_CUDA}
         ${SOURCES_SPECFEMINTERFACE} ${SOURCES_FWAT_SHARED})
  endif()
  # add executable
  add_executable(${EXEC_NAME} ${execsourcefile}
    ${SOURCES_SHARED} ${SOURCES_SUBS} 
    ${HEADERS}
  )

  set_target_properties(
    ${EXEC_NAME}
    PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod_${EXEC_NAME}"
    PREFIX "x"
)
  # link MPI
  target_link_libraries(${EXEC_NAME} PUBLIC MPI::MPI_Fortran)

  # link LAPACK
  target_link_libraries(${EXEC_NAME} PUBLIC LAPACK::LAPACK)

  # link BLAS
  target_link_libraries(${EXEC_NAME} PUBLIC BLAS::BLAS)
  
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