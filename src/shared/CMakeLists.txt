file(GLOB SOURCES_Fortran "*.f90")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
add_library(shared_fwat_objs OBJECT ${SOURCES_Fortran})
