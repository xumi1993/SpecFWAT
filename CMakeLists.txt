cmake_minimum_required(VERSION 3.10)

set(CMAKE_PROJECT_NAME "SpecFWI")
project(${CMAKE_PROJECT_NAME}, VERSION 0.1.0 LANGUAGES C CXX Fortran)

execute_process(
    COMMAND git submodule update --init
)

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
elseif("${CMAKE_FORTRAN_COMPILER_ID}" MATCHES "Flang")
  set(CMAKE_Fortran_FLAGS "-O3 -fimplicit-none -ffpe-trap=invalid,zero,overflow -funroll-loops -Wall")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fbacktrace -fbounds-check")
endif()

# setup mpi
option(USE_SERIAL "Enable serial mode" OFF)
if (NOT USE_SERIAL)
  find_package(MPI REQUIRED)
  message(STATUS "MPI_FOUND: ${MPI_FOUND}")
  message(STATUS "MPI_VERSION: ${MPI_VERSION}")
else()
  message(STATUS "SpecFWI is compiled in serial mode, because -DUSE_SERIAL=True is defined")
endif()

# setup CUDA
if(USE_CUDA)
  message(STATUS "Running CMAKE FindCUDA.cmake...")
  enable_language(CUDA)
  find_package(CUDA)
else()
  message(STATUS "SpecFWI is compiled without cuda, because -DUSE_CUDA=True is not defined")
endif()

option(FORCE_DOWNLOAD_EXTERNAL_LIBS "Force download and use external libraries" OFF)
# yaml parser
find_package(yaml-cpp 0.8 QUIET)
if (yaml-cpp_FOUND AND NOT ${FORCE_DOWNLOAD_EXTERNAL_LIBS})
  message(STATUS "yaml-cpp found")
  message(STATUS "YAML_CPP_INCLUDE_DIR: ${YAML_CPP_INCLUDE_DIR}")
  message(STATUS "YAML_CPP_LIBRARIES: ${YAML_CPP_LIBRARIES}")
else()
  message(STATUS "yaml-cpp not found. Using external_libs/yaml-cpp ...")
  add_subdirectory(external_libs)
  set(YAML_CPP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/external_libs/yaml-cpp/include)
  set(YAML_CPP_LIBRARIES yaml-cpp)
endif()

# setup version
execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/specfem3d 
    RESULT_VARIABLE GIT_REV_PARSE_STATUS
    OUTPUT_VARIABLE GIT_COMMIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
file(READ "${PROJECT_SOURCE_DIR}/src/specfem3d/VERSION" GIT_PACKAGE_VERSION)
string(REGEX MATCH "^([^\n]*)" GIT_PACKAGE_VERSION "${GIT_PACKAGE_VERSION}")
execute_process(
    COMMAND git log -1 --format=%cd --date=short
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/specfem3d
    RESULT_VARIABLE GIT_LOG_STATUS
    OUTPUT_VARIABLE GIT_DATE_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
configure_file(
  ${PROJECT_SOURCE_DIR}/src/specfem3d/setup/version.fh.in
  ${PROJECT_SOURCE_DIR}/include/version.fh
)

# add -DENABLE_DOUBLE_PRECISION=ON to enable double precision
option(ENABLE_DOUBLE_PRECISION "Enable double precision for solver" OFF)
if (ENABLE_DOUBLE_PRECISION)
    set(CUSTOM_REAL "SIZE_DOUBLE")
    set(CUSTOM_MPI_TYPE "MPI_DOUBLE_PRECISION")
else()
    set(CUSTOM_REAL "SIZE_REAL")
    set(CUSTOM_MPI_TYPE "MPI_REAL")
endif()
if (NOT DEFINED LOCAL_PATH_IS_ALSO_GLOBAL)
    set(LOCAL_PATH_IS_ALSO_GLOBAL "true")
endif()
# setup constants
configure_file(
  ${PROJECT_SOURCE_DIR}/src/specfem3d/setup/constants.h.in
  ${PROJECT_SOURCE_DIR}/include/constants.h
)
# setup precision
configure_file(
  ${PROJECT_SOURCE_DIR}/src/specfem3d/setup/precision.h.in
  ${PROJECT_SOURCE_DIR}/include/precision.h
)

# config config.fh
set(PACKAGE_STRING "Specfem3D")
set(CONFIGURE_FLAGS "FC=${CMAKE_Fortran_COMPILER} FCFLAGS=${CMAKE_Fortran_FLAGS}")
configure_file(
  ${PROJECT_SOURCE_DIR}/src/specfem3d/setup/config.fh.in
  ${PROJECT_SOURCE_DIR}/include/config.fh
)

# setup adios
option(USE_ADIOS "Enable ADIOS support" OFF)
option(USE_ADIOS2 "Enable ADIOS2 support" OFF)
if (USE_ADIOS)
  add_definitions(-DUSE_ADIOS)
endif()
if (USE_ADIOS2)
  add_definitions(-DUSE_ADIOS2)
endif()

# setup asdf
option(USE_ASDF "Enable ASDF support" OFF)
if (USE_ASDF)
  add_definitions(-DUSE_ASDF)
endif()

# add specfem3d library and executable
add_subdirectory(src)