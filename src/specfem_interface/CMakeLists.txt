# object libraries of specfem_interface
file(GLOB SOURCES_SPECFEM_API "*.f90")
add_library(specfem_api_objs OBJECT ${SOURCES_SPECFEM_API})
target_link_libraries(specfem_api_objs PRIVATE specfem3D_objs shared_fwat_objs)
