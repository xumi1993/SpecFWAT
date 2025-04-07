file(GLOB SOURCES_Fortran "*.f90")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
add_library(shared_fwat_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(shared_fwat_objs PUBLIC shared_objs ${MPI_LIBRARIES} 
                      ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
                      fortran-yaml ${YAML_CPP_LIBRARIES}
)
