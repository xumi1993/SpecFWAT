set(PROJECT_SOURCES "${PROJECT_SOURCE_DIR}/src/grid/projection_on_FD_grid_mod_fwat.f90")
add_library(project_lib OBJECT ${PROJECT_SOURCES})
target_link_libraries(project_lib PRIVATE specfem3D_objs shared_objs shared_fwat_objs)

set(APP_SOURCES
    "decompose_h5_gll.f90"
    "model_grid_cart_custom.f90"
)

foreach(execsourcefile ${APP_SOURCES})
    get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)
    if (EXEC_NAME MATCHES "decompose_h5_gll" AND NOT ${USE_HDF5})
        message(STATUS "decompose_h5_gll is not compiled because -DUSE_HDF5=True is not defined")
        continue()
    endif()
    add_executable(${EXEC_NAME} ${execsourcefile}
    )
    set_target_properties(${EXEC_NAME} PROPERTIES PREFIX "x")
    set_target_properties(${EXEC_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
    target_link_libraries(${EXEC_NAME} PRIVATE specfem3D_objs shared_objs shared_fwat_objs
                                        cuda_objs fortran-yaml project_lib
                                        ${MPI_LIBRARIES} ${YAML_CPP_LIBRARIES})
    if (USE_CUDA)
        target_link_libraries(${EXEC_NAME} PRIVATE cuda_kernels_objs)
    endif()
    if (USE_HDF5)
        target_link_libraries(${EXEC_NAME} PRIVATE ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    endif()
endforeach()