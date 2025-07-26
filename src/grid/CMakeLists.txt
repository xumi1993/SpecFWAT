file(GLOB SOURCES_Fortran "*.f90")
add_library(grid_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(grid_objs PRIVATE shared_objs
                      cuda_objs specfem3D_objs 
                      generate_databases_objs
                      hdf5_interface_objs utils_objs  
)