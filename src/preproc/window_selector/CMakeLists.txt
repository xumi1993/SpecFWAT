add_subdirectory(ttimes_mod)

file(GLOB window_selector_sources "*.f90")
add_library(window_selector_objs OBJECT ${window_selector_sources})

target_link_libraries(window_selector_objs PRIVATE config_objs utils_objs
                      fftpack_objs shared_fwat_objs tau)