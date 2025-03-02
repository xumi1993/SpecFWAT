# objs
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
set(SPECFEM_DIR ${PROJECT_SOURCE_DIR}/src/specfem3d)

# shared_objs
file(GLOB SOURCES_Fortran "${SPECFEM_DIR}/src/shared/*.[Ff]90")
file(GLOB SOURCES_C "${SPECFEM_DIR}/src/shared/*.c")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${SPECFEM_DIR}/src/shared/merge_sort.f90")
# if(USE_SERIAL)
  # list(FILTER SOURCES_Fortran EXCLUDE REGEX "${SPECFEM_DIR}/src/shared/parallel.f90")
# else()
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${SPECFEM_DIR}/src/shared/serial.f90")
set(MPI_LIBRARIES MPI::MPI_Fortran)
# endif()
if(NOT (USE_ADIOS AND USE_ADIOS2))
  list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/adios_.*$")
  list(FILTER SOURCES_C EXCLUDE REGEX "^.*/adios_.*$")
  list(APPEND SOURCES_Fortran "${SPECFEM_DIR}/src/shared/adios_manager.F90")
  list(APPEND SOURCES_C "${SPECFEM_DIR}/src/shared/adios_method_stubs.c")
endif()
if (USE_ASDF)
  list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/asdf_.*_stubs.*$")
  list(FILTER SOURCES_C EXCLUDE REGEX "^.*/asdf_.*_stubs.*$")
else()
  list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/asdf_manager.*$")
endif()
add_library(shared_objs OBJECT ${SOURCES_Fortran} ${SOURCES_C})
target_link_libraries(shared_objs PRIVATE ${MPI_LIBRARIES})
if (USE_HDF5)
  target_link_libraries(shared_objs PRIVATE hdf5::hdf5_fortran)
endif()

# CUDA
if (USE_CUDA)
  include_directories(
    ${SPECFEM_DIR}/src/gpu
    ${SPECFEM_DIR}/src/gpu/kernels
  )

  if (CUDA_ARCH EQUAL 12)
    set(CMAKE_CUDA_ARCHITECTURES 90)
  elseif (CUDA_ARCH EQUAL 11)
    set(CMAKE_CUDA_ARCHITECTURES 80)
  elseif (CUDA_ARCH EQUAL 10)
    set(CMAKE_CUDA_ARCHITECTURES 75)
  elseif (CUDA_ARCH EQUAL 9)
    set(CMAKE_CUDA_ARCHITECTURES 70)
  elseif (CUDA_ARCH EQUAL 8)
    set(CMAKE_CUDA_ARCHITECTURES 60)
  elseif (CUDA_ARCH EQUAL 7)
    set(CMAKE_CUDA_ARCHITECTURES 50)
  elseif (CUDA_ARCH EQUAL 6)
    set(CMAKE_CUDA_ARCHITECTURES 37)
  elseif (CUDA_ARCH EQUAL 5)
    set(CMAKE_CUDA_ARCHITECTURES 35)
  elseif (CUDA_ARCH EQUAL 4)
    set(CMAKE_CUDA_ARCHITECTURES 30)
  # else()
  #   set(CMAKE_CUDA_ARCHITECTURES 20)
  endif()
  file(GLOB SOURCES_CUDA "${SPECFEM_DIR}/src/gpu/*.cu")
  set(CUDA_KERNEL_DIR "${SPECFEM_DIR}/src/gpu/kernels")
  file(GLOB SOURCES_CUDA_KERNELS "${CUDA_KERNEL_DIR}/*.cu")
  add_library(cuda_objs OBJECT ${SOURCES_CUDA})
  add_library(cuda_kernels_objs OBJECT ${SOURCES_CUDA_KERNELS})
  target_compile_options(cuda_kernels_objs PRIVATE -include ${SPECFEM_DIR}/src/gpu/mesh_constants_gpu.h)
  message(STATUS "CUDA_ARCH: ${CUDA_ARCH}")
  if (${CUDA_ARCH} EQUAL 0)
    foreach(tgt cuda_objs cuda_kernels_objs)
      set_target_properties(${tgt} PROPERTIES CUDA_ARCHITECTURES native)
      target_link_libraries(${tgt} PRIVATE ${MPI_LIBRARIES})
    endforeach()
  endif()
else()
  set(SOURCES_CUDA "${SPECFEM_DIR}/src/gpu/specfem3D_gpu_cuda_method_stubs.c")
  add_library(cuda_objs OBJECT ${SOURCES_CUDA})
endif()

# meshfem3D_objs
file(GLOB SOURCES_Fortran "${SPECFEM_DIR}/src/meshfem3D/*.[Ff]90")
file(GLOB SOURCES_C "${SPECFEM_DIR}/src/meshfem3D/*.c")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${SPECFEM_DIR}/src/meshfem3D/meshfem3D.F90")
if(NOT (USE_ADIOS AND USE_ADIOS2))
  list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/.*adios.*$")
endif()
add_library(meshfem3D_objs OBJECT ${SOURCES_Fortran} ${SOURCES_C})
target_link_libraries(meshfem3D_objs PRIVATE shared_objs)

# generate_databases_objs
file(GLOB SOURCES_Fortran "${SPECFEM_DIR}/src/generate_databases/*.[Ff]90")
file(GLOB SOURCES_C "${SPECFEM_DIR}/src/generate_databases/*.c")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${SPECFEM_DIR}/src/generate_databases/generate_databases.f90")
# if(USE_SERIAL)
  # list(FILTER SOURCES_Fortran EXCLUDE REGEX "${SPECFEM_DIR}/src/generate_databases/model_sep.f90")
# else()
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${SPECFEM_DIR}/src/generate_databases/model_sep_nompi.f90")
# endif()
if(NOT (USE_ADIOS AND USE_ADIOS2))
  list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/.*adios.*$")
endif()
add_library(generate_databases_objs OBJECT ${SOURCES_Fortran} ${SOURCES_C})
target_link_libraries(generate_databases_objs PRIVATE shared_objs)

# specfem3D_objs
file(GLOB SOURCES_Fortran "${SPECFEM_DIR}/src/specfem3D/*.[Ff]90")
file(GLOB SOURCES_C "${SPECFEM_DIR}/src/specfem3D/*.c*")
list(FILTER SOURCES_Fortran EXCLUDE REGEX "${SPECFEM_DIR}/src/specfem3D/specfem3D.F90")
if(NOT (USE_ADIOS AND USE_ADIOS2))
  list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/.*adios.*$")
endif()
if (NOT USE_ASDF)
  list(FILTER SOURCES_Fortran EXCLUDE REGEX "^.*/.*ASDF.*$")
endif()
add_library(specfem3D_objs OBJECT ${SOURCES_Fortran} ${SOURCES_C})
target_link_libraries(specfem3D_objs PRIVATE shared_objs)

# ==============================================================================
# ----------------------------- Add executables -------------------------------
# ==============================================================================
set(APP_SOURCES
  "${SPECFEM_DIR}/src/meshfem3D/meshfem3D.F90"
  "${SPECFEM_DIR}/src/generate_databases/generate_databases.f90"
  "${SPECFEM_DIR}/src/specfem3D/specfem3D.F90"
) 
foreach(execsourcefile ${APP_SOURCES})
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)
  add_executable(${EXEC_NAME} ${execsourcefile}
    ${SOURCES_SUBS} 
    ${HEADERS}
  )
  set_target_properties(${EXEC_NAME} PROPERTIES PREFIX "x")
  target_link_libraries(${EXEC_NAME} PUBLIC "${EXEC_NAME}_objs" shared_objs cuda_objs MPI::MPI_Fortran)
  set_target_properties(${EXEC_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
  if (USE_CUDA)
    target_link_libraries(${EXEC_NAME} PUBLIC cuda_kernels_objs)
  endif()
  if (USE_HDF5)
    target_link_libraries(${EXEC_NAME} PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
  endif()
endforeach()

# ==============================================================================
# -------------------------- Add Specfem AUX -----------------------------------
# ==============================================================================
set(AUX_DIR ${SPECFEM_DIR}/src/auxiliaries)
set(APP_SOURCES
  "${AUX_DIR}/combine_vol_data_vtk_binary.F90"
) 
foreach(execsourcefile ${APP_SOURCES})
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)
  add_executable(${EXEC_NAME} ${execsourcefile}
    ${AUX_DIR}/combine_vol_data_impl.f90
    ${AUX_DIR}/vtk_writer.c
    ${HEADERS}
  )
  target_link_libraries(${EXEC_NAME} PUBLIC shared_objs)
  set_target_properties(${EXEC_NAME} PROPERTIES PREFIX "x")
  set_target_properties(${EXEC_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
endforeach()
# ==============================================================================
# ----------------------------- Add libraries --------------------------------
# ==============================================================================
# add yaml-fortran library
add_subdirectory(fortran-yaml)
add_subdirectory(shared)

# ==============================================================================
# ------------------------ Add FWAT executables --------------------------------
# ==============================================================================
add_subdirectory(preproc)
add_subdirectory(aux)