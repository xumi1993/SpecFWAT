file(GLOB SOURCES_Fortran "*.f90")
add_library(grid_objs OBJECT ${SOURCES_Fortran})
list(APPEND SOURCES_Fortran
    "${PROJECT_SOURCE_DIR}/src/shared/hdf5_interface.f90"
    "${PROJECT_SOURCE_DIR}/src/shared/utils.f90"
)
target_link_libraries(grid_objs PRIVATE shared_objs
                      cuda_objs specfem3D_objs)