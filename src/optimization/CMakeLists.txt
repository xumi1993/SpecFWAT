# object libraries of preproc
file(GLOB SOURCES_Fortran "*.f90")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/fwat.*$")
add_library(opt_objs OBJECT ${SOURCES_Fortran})
target_link_libraries(opt_objs PRIVATE shared_objs shared_fwat_objs
                      cuda_objs)

                      