file(GLOB SOURCES "*.f90")
file(GLOB SOURCES_UTILS "utils/*.f90")

add_library(adjsrc_objs OBJECT ${SOURCES_UTILS} ${SOURCES})
target_link_libraries(adjsrc_objs PRIVATE shared_fwat_objs config_objs utils_objs tele_objs)