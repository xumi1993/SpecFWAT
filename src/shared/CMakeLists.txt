include_directories(
  ${PROJECT_SOURCE_DIR}/include
)

# Configuration file
set(config_src "${PROJECT_SOURCE_DIR}/src/shared/config.f90")
add_library(config_objs OBJECT ${config_src})
target_link_libraries(config_objs PRIVATE shared_objs)

# HDF5 interface
set(hdf5_interface_src "${PROJECT_SOURCE_DIR}/src/shared/hdf5_interface.f90")
add_library(hdf5_interface_objs OBJECT ${hdf5_interface_src})
target_link_libraries(hdf5_interface_objs PRIVATE
                      ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
)

# Utility functions
set(utils_src "${PROJECT_SOURCE_DIR}/src/shared/utils.f90")
add_library(utils_objs OBJECT ${utils_src})
target_link_libraries(utils_objs PRIVATE shared_objs config_objs)

# FFT package
set(FFT_SRC_DIR "${PROJECT_SOURCE_DIR}/src/shared/fftpack.f90")
add_library(fftpack_objs OBJECT ${FFT_SRC_DIR})

file(GLOB SOURCES_Fortran "*.f90")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/shared/utils.f90")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/shared/hdf5_interface.f90")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/shared/config.f90")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/shared/aniso.f90")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/shared/fftpack.f90")

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set_source_files_properties(sacio.f90 PROPERTIES COMPILE_OPTIONS -Wno-pedantic)
endif()

add_library(shared_fwat_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(shared_fwat_objs PUBLIC shared_objs ${MPI_LIBRARIES} 
                      ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
                      fortran-yaml ${YAML_CPP_LIBRARIES} config_objs
                      hdf5_interface_objs utils_objs
)
