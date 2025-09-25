file(GLOB SOURCES_TELE "*.f90")

add_library(tele_objs OBJECT ${SOURCES_TELE})
target_link_libraries(tele_objs PRIVATE fftpack_objs LAPACK::LAPACK BLAS::BLAS)