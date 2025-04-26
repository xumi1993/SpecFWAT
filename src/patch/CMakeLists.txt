file(GLOB SOURCES_Fortran "*.f90")
list(APPEND SOURCES_Fortran
    "${PROJECT_SOURCE_DIR}/src/specfem3d/src/generate_databases/generate_databases_par.F90"
)
list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/fwat.*$")
add_library(external_values_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(external_values_objs PRIVATE shared_objs
                      utils_objs hdf5_interface_objs config_objs
                                          
)
