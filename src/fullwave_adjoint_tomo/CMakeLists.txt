file(GLOB SOURCES_FWAT_RECFUNC "${SOURCE_FWAT_DIR}/preproc_recfunc/*.f90" )
file(GLOB SOURCES_TELESTF "${SOURCE_FWAT_DIR}/teleseis_stf/telestf_mod.f90"
                          "${SOURCE_FWAT_DIR}/teleseis_stf/spanlib.f90"
                          "${SOURCE_FWAT_DIR}/teleseis_stf/fktimes.f90")
set(SOURCES_FWAT_SHARED "${SOURCE_FWAT_DIR}/fullwave_adjoint_tomo_par.f90"
                        "${SOURCE_FWAT_DIR}/fwat_input.f90"
                        "${SOURCE_FWAT_DIR}/fwat_utils.f90"
                        "${SOURCE_FWAT_DIR}/sacio.f90")
set(SOURCES_PRE_PROC "${SOURCE_FWAT_DIR}/preproc_subs.f90"
                     "${SOURCE_FWAT_DIR}/run_preprocessing.f90"
                     "${SOURCE_FWAT_DIR}/preproc_measure_adj_subs.f90"
                     "${SOURCE_FWAT_DIR}/preproc_collect_data.f90")
set(SOURCE_POST_PROC "${SOURCE_FWAT_DIR}/run_linesearch.f90"
                     "${SOURCE_FWAT_DIR}/postproc_subs.f90"
                     "${SOURCE_FWAT_DIR}/post_processing/taper_3D_mod.f90"
                     "${SOURCE_FWAT_DIR}/run_postprocessing.f90"
                     "${SOURCE_FWAT_DIR}/run_optimization.f90")
file(GLOB SOURCES_FWAT_OPT "${SOURCE_FWAT_DIR}/optimization/*.[fF]90")


include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/modules
)

# Add submodules
add_subdirectory(grid_interface)
add_subdirectory(preproc_measure_adj)
add_subdirectory(preproc_flexwin)
add_subdirectory(specfem_interface)


add_library(fwat_interpolation_objs OBJECT ${SOURCE_FWAT_DIR}/teleseis_stf/interpolation_mod.f90)
add_library(fwat_shared_objs OBJECT ${SOURCES_FWAT_SHARED})
target_link_libraries(fwat_shared_objs PRIVATE shared_objs specfem3D_objs MPI::MPI_Fortran)
add_library(fwat_recfunc_objs OBJECT ${SOURCES_FWAT_RECFUNC})
add_library(fwat_telestf_objs OBJECT ${SOURCES_TELESTF})
target_link_libraries(fwat_telestf_objs PRIVATE fwat_shared_objs fwat_interpolation_objs MPI::MPI_Fortran)
add_library(fwat_opt_objs OBJECT ${SOURCES_FWAT_OPT})
target_link_libraries(fwat_opt_objs PRIVATE fwat_shared_objs shared_objs MPI::MPI_Fortran)
add_library(fwat_post_proc_objs OBJECT ${SOURCE_POST_PROC})
target_link_libraries(fwat_post_proc_objs PRIVATE fwat_shared_objs shared_objs fwat_interpolation_objs MPI::MPI_Fortran)
add_library(fwat_pre_proc_objs OBJECT ${SOURCES_PRE_PROC})
target_link_libraries(fwat_pre_proc_objs PUBLIC fwat_shared_objs shared_objs
                      measadj_objs fwat_recfunc_objs flexwin_objs specfem_interface_objs
                      fwat_interpolation_objs fwat_telestf_objs fwat_opt_objs MPI::MPI_Fortran)

# Define fwat0 library
set(SOURCES_fwat0 "run_fwat0_forward_data.f90"
                  "run_semd2sac.f90")
add_library(fwat0_forward_data_objs OBJECT ${SOURCES_fwat0})
target_link_libraries(fwat0_forward_data_objs PUBLIC measadj_objs flexwin_objs 
                      specfem_interface_objs read_mesh_interface_objs
                      fwat_shared_objs fwat_recfunc_objs fwat_interpolation_objs
                      MPI::MPI_Fortran)

# Define fwat1 library
set(SOURCES_fwat1 "run_fwat1_fwd_measure_adj.f90"
                  "run_precond.f90")
add_library(fwat1_fwd_measure_adj_objs OBJECT ${SOURCES_fwat1})
target_link_libraries(fwat1_fwd_measure_adj_objs PUBLIC measadj_objs flexwin_objs 
                      specfem_interface_objs read_mesh_interface_objs fwat_telestf_objs
                      fwat_shared_objs fwat_recfunc_objs fwat_pre_proc_objs
                      MPI::MPI_Fortran)

# Define fwat2 library
set(SOURCES_fwat2 "run_fwat2_postproc_opt.f90")
add_library(fwat2_postproc_opt_objs OBJECT ${SOURCES_fwat2})
target_link_libraries(fwat2_postproc_opt_objs PUBLIC
                      read_mesh_interface_objs fwat_opt_objs measadj_objs flexwin_objs 
                      fwat_post_proc_objs fwat_interpolation_objs fwat_telestf_objs
                      fwat_pre_proc_objs MPI::MPI_Fortran)

set(SOURCES_fwat3 "run_fwat3_optim.f90")
add_library(fwat3_linesearch_objs OBJECT ${SOURCES_fwat3})
target_link_libraries(fwat3_linesearch_objs PUBLIC fwat_mesh_database_objs
                      fwat_opt_objs fwat_post_proc_objs fwat_telestf_objs
                      fwat_pre_proc_objs  MPI::MPI_Fortran)

# Define fwat_main library
set(SOURCES_fwat_main "fullwave_adjoint_tomo_main.f90")
add_library(fullwave_adjoint_tomo_objs OBJECT ${SOURCES_fwat_main})
target_link_libraries(fullwave_adjoint_tomo_objs PUBLIC fwat_mesh_database_objs
                      fwat1_fwd_measure_adj_objs fwat2_postproc_opt_objs fwat3_linesearch_objs)