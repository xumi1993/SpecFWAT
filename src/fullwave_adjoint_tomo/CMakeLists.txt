file(GLOB SOURCES_FWAT_RECFUNC "${SOURCE_FWAT_DIR}/preproc_recfunc/*.f90" )
file(GLOB SOURCES_TELESTF "${SOURCE_FWAT_DIR}/teleseis_stf/*.f90")
set(SOURCES_FWAT_SHARED "${SOURCE_FWAT_DIR}/fullwave_adjoint_tomo_par.f90"
                        "${SOURCE_FWAT_DIR}/fwat_input.f90"
                        "${SOURCE_FWAT_DIR}/fwat_utils.f90"
                        "${SOURCE_FWAT_DIR}/sacio.f90")
file(GLOB SOURCES_FWAT_OPT "${SOURCE_FWAT_DIR}/optimization/*.[fF]90")


add_library(fwat_shared_objs OBJECT ${SOURCES_FWAT_SHARED})
target_link_libraries(fwat_shared_objs PRIVATE shared_objs specfem3D_objs)
add_library(fwat_recfunc_objs OBJECT ${SOURCES_FWAT_RECFUNC})
add_library(fwat_telestf_objs OBJECT ${SOURCES_TELESTF})
target_link_libraries(fwat_telestf_objs PRIVATE fwat_shared_objs)

# Add submodules
add_subdirectory(preproc_measure_adj)
add_subdirectory(preproc_flexwin)
add_subdirectory(specfem_interface)


# Define fwat0 library
set(SOURCES_fwat0 "run_fwat0_forward_data.f90"
                  "run_semd2sac.f90")
add_library(fwat0_forward_data_objs OBJECT ${SOURCES_fwat0})
target_link_libraries(fwat0_forward_data_objs PUBLIC measadj_objs flexwin_objs 
                      specfem_interface_objs read_mesh_interface_objs
                      fwat_shared_objs fwat_recfunc_objs)

# Define fwat1 library
set(SOURCES_fwat1 "run_fwat1_fwd_measure_adj.f90"
                  "preproc_subs.f90"
                  "run_preprocessing.f90"
                  "preproc_measure_adj_subs.f90"
                  "preproc_collect_data.f90"
                  "run_precond.f90")
add_library(fwat1_fwd_measure_adj_objs OBJECT ${SOURCES_fwat1})
target_link_libraries(fwat1_fwd_measure_adj_objs PUBLIC measadj_objs flexwin_objs 
                      specfem_interface_objs read_mesh_interface_objs
                      fwat_shared_objs fwat_recfunc_objs)