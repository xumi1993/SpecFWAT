# Compile the shared FK backends once so preproc and xfwat_fk use the same implementation.
set(FK_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/fk_coupling.f90)
# Enable CPU/GPU preprocessing in the directory owning the shared target.
set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/fk_coupling.f90
  TARGET_DIRECTORY shared_fwat_objs PROPERTIES Fortran_PREPROCESS ON)
include(reference.cmake)
target_sources(shared_fwat_objs PRIVATE ${FK_SOURCES})
# The shared backends import specfem_par; make its module available in clean parallel builds.
target_link_libraries(shared_fwat_objs PRIVATE specfem3D_objs)
if(USE_CUDA)
  target_sources(shared_fwat_objs PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/gpu/fk_gpu_backend.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/gpu/fk_cuda.cu)
  target_link_libraries(shared_fwat_objs PUBLIC CUDA::cufft CUDA::cudart)
  if(ENABLE_DOUBLE_PRECISION)
    target_compile_definitions(shared_fwat_objs PRIVATE FK_DOUBLE_PRECISION)
  endif()
endif()

add_executable(fwat_fk prepare_fk.f90 fwat_fk_compute.f90)
set_target_properties(fwat_fk PROPERTIES PREFIX "x" LINKER_LANGUAGE Fortran)
target_link_libraries(fwat_fk PRIVATE specfem_api_objs specfem3D_objs shared_objs shared_fwat_objs
  cuda_objs config_objs hdf5_interface_objs utils_objs fortran-yaml
  ${MPI_LIBRARIES} ${YAML_CPP_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
if(USE_CUDA)
  target_link_libraries(fwat_fk PRIVATE cuda_kernels_objs)
endif()
