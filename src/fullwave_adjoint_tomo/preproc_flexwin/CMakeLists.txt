# objs
add_subdirectory(ttimes_mod)

file(GLOB SOURCES_Fortran "*.f90")
file(GLOB SOURCES_C "*.c")

include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${SOURCE_FWAT_DIR}/sac_inc
  ${CMAKE_BINARY_DIR}/modules
)
add_library(flexwin_objs OBJECT ${SOURCES_Fortran} ${SOURCES_C})
target_link_libraries(flexwin_objs PRIVATE tau fwat_shared_objs)
