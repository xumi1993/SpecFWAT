#=====================================================================
#
#               Full Waveform Adjoint Tomography -v1.1
#               ---------------------------------------
#
#     Main historical authors: Kai Wang
#                              Macquarie Uni, Australia
#                            & University of Toronto, Canada
#                           (c) Martch 2020
#
#=====================================================================
#
### compilation directories
S := ${S_TOP}/src/fullwave_adjoint_tomo
$(fullwave_adjoint_tomo_OBJECTS): S = ${S_TOP}/src/fullwave_adjoint_tomo

#######################################
# solver objects - no statically allocated arrays anymore

####
#### targets
####

# default targets for the pure Fortran version
fullwave_adjoint_tomo_TARGETS = \
	$E/xfullwave_adjoint_tomo \
	$E/xfwat1_fwd_measure_adj \
	$E/xfwat2_postproc_opt \
	$E/xdynmbat_sele_ctrlgrp \
	$E/xfwat3_linesearch \
	$E/xfwat0_forward_data \
	$E/xfwat4_cmt3d_fwd \
	$E/xfwat5_rtm \
    $E/xmeasure_adj \
    $E/xflexwin \
    $E/xdynmbat_sele_minibatch \
	$E/xgrid_and_combine_vol_data \
	$E/xgrid_and_combine_vol_data_cert \
	$E/xdecomp_grid2gll \
	$E/xadd_pert_gaus \
	$E/xadd_pert_trigo \
	$E/xfwat_mesh_database \
	$(EMPTY_MACRO)


## source folder objects
fwat_top_OBJECTS = \
	$O/fwat_utils.fwat_utils.o \
	$O/run_fwat1_fwd_measure_adj.fwat_fwdadj.o \
	$O/fwat_input.fwat_readpar.o \
	$O/specfem_initialize_simulation_fwat.fwat_spec.o \
	$O/preproc_collect_data.fwat_preproc.o \
	$O/preproc_subs.fwat_preproc.o \
	$O/preproc_measure_adj_subs.fwat_preproc.o \
	$O/fftpack.fwat_deconit.o \
	$O/deconit.fwat_deconit.o \
	$O/run_preprocessing.fwat_preproc.o \
	$O/run_precond.fwat_precond.o \
	$O/run_semd2sac.fwat_semd2sac.o \
	$O/interpolation_mod.fwat_interp.o \
	$O/spanlib.fwat_pca.o \
	$O/telestf_mod.fwat_telestf.o \
	$O/specfem_read_mesh_database.fwat_sr.o \
	$O/specfem_setup_sources_receivers_fwat.fwat_sr.o \
	$O/specfem_save_adjoint_kernels_fwat.fwat_savekern.o \
	$O/specfem_interface.fwat_spec.o \
	$O/specfem_prepare_timerun_fwat.fwat_spec.o \
	$O/fullwave_adjoint_tomo_par.fwat_par.o \
	$(EMPTY_MACRO)

preproc_measure_adj_OBJECTS = \
	$O/xapiir_sub.fwat_preproc_meas.o \
	$O/sacio.fwat_preproc_meas.o \
	$O/ascii_rw.fwat_preproc_meas.o \
	$O/ma_sub.fwat_preproc_meas.o \
	$O/ma_sub2.fwat_preproc_meas.o \
	$O/measure_adj.fwat_preproc_meas.o \
	$(EMPTY_MACRO)

preproc_flexwin_OBJECTS = \
	$O/dcops.fwat_preproc_flex.o \
	$O/nextpow2.fwat_preproc_flex.o \
	$O/cfft.fwat_preproc_flex.o \
	$O/hilenv.fwat_preproc_flex.o\
	$O/hec_fortran_wrapper.fwat_preproc_flex.o \
	$O/distaz.fwat_preproc_flex.o \
	$O/seismo_subs.fwat_preproc_flex.o \
	$O/user_functions.fwat_preproc_flex.o \
	$O/maxima.fwat_preproc_flex.o \
	$O/travel_times.fwat_preproc_flex.o \
	$O/io_subs.fwat_preproc_flex.o \
	$O/select_windows_stalta2.fwat_preproc_flex.o \
	$(EMPTY_MACRO)


fullwave_adjoint_tomo_OBJECTS = \
	$(fwat_top_OBJECTS) \
	$(preproc_measure_adj_OBJECTS) \
	$(preproc_flexwin_OBJECTS) \
	$(EMPTY_MACRO)

## objects from other source directories
fullwave_adjoint_tomo_OBJECTS += \
	$O/specfem3D_par.spec_module.o \
	$O/asdf_data.spec_module.o \
	$O/assemble_MPI_vector.spec.o \
	$O/check_stability.spec.o \
	$O/comp_source_time_function.spec.o \
	$O/compute_add_sources_acoustic.spec.o \
	$O/compute_add_sources_viscoelastic.spec.o \
	$O/compute_add_sources_poroelastic.spec.o \
	$O/compute_adj_source_frechet.spec.o \
	$O/compute_arrays_source.spec.o \
	$O/compute_boundary_kernel.spec.o \
	$O/compute_coupling_acoustic_el.spec.o \
	$O/compute_coupling_acoustic_po.spec.o \
	$O/compute_coupling_viscoelastic_ac.spec.o \
	$O/compute_coupling_viscoelastic_po.spec.o \
	$O/compute_coupling_poroelastic_ac.spec.o \
	$O/compute_coupling_poroelastic_el.spec.o \
	$O/compute_forces_acoustic_calling_routine.spec.o \
	$O/compute_forces_acoustic_NGLL5_fast.spec.o \
	$O/compute_forces_acoustic_NGLLnot5_generic_slow.spec.o \
	$O/compute_forces_viscoelastic_calling_routine.spec.o \
	$O/compute_forces_viscoelastic.spec.o \
	$O/compute_element_att_memory.spec.o \
	$O/compute_forces_poro_fluid_part.spec.o \
	$O/compute_forces_poroelastic_calling_routine.spec.o \
	$O/compute_forces_poro_solid_part.spec.o \
	$O/compute_gradient_in_acoustic.spec.o \
	$O/compute_interpolated_dva.spec.o \
	$O/compute_kernels.spec.o \
	$O/compute_seismograms.spec.o \
	$O/compute_stacey_acoustic.spec.o \
	$O/compute_stacey_viscoelastic.spec.o \
	$O/compute_stacey_poroelastic.spec.o \
	$O/compute_energy.spec.o \
	$O/convert_time.spec.o \
	$O/couple_with_injection.spec.o \
	$O/calendar.spec.o \
	$O/create_color_image.spec.o \
	$O/detect_mesh_surfaces.spec.o \
	$O/fault_solver_common.spec.o \
	$O/fault_solver_dynamic.spec.o \
	$O/fault_solver_kinematic.spec.o \
	$O/finalize_simulation.spec.o \
	$O/get_cmt.spec.o \
	$O/get_elevation.spec.o \
	$O/get_force.spec.o \
	$O/gravity_perturbation.spec.o \
	$O/initialize_simulation.spec.o \
	$O/iterate_time.spec.o \
	$O/locate_MPI_slice.spec.o \
	$O/locate_point.spec.o \
	$O/locate_receivers.spec.o \
	$O/locate_source.spec.o \
	$O/make_gravity.spec.o \
	$O/noise_tomography.spec.o \
	$O/pml_allocate_arrays.spec.o \
	$O/pml_output_VTKs.spec.o \
	$O/pml_compute_accel_contribution.spec.o \
	$O/pml_compute_memory_variables.spec.o \
	$O/pml_par.spec.o \
	$O/prepare_attenuation.spec.o \
	$O/prepare_gpu.spec.o \
	$O/prepare_gravity.spec.o \
	$O/prepare_noise.spec.o \
	$O/prepare_timerun.spec.o \
	$O/prepare_wavefields.spec.o \
	$O/print_stf_file.spec.o \
	$O/read_external_stf.spec.o \
	$O/read_mesh_databases.spec.o \
	$O/save_adjoint_kernels.spec.o \
	$O/setup_GLL_points.spec.o \
	$O/setup_movie_meshes.spec.o \
	$O/setup_sources_receivers.spec.o \
	$O/station_filter.spec.o \
	$O/surface_or_volume_integral.spec.o \
	$O/update_displacement_scheme.spec.o \
	$O/update_displacement_LDDRK.spec.o \
	$O/write_movie_output.spec.o \
	$O/write_output_ASCII_or_binary.spec.o \
	$O/write_output_SU.spec.o \
	$O/write_seismograms.spec.o \
	$(EMPTY_MACRO)


fullwave_adjoint_tomo_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/meshfem3D_par.mesh.o \
	$O/assemble_MPI_scalar.shared.o \
	$O/check_mesh_resolution.shared.o \
	$O/create_name_database.shared.o \
	$O/define_derivation_matrices.shared.o \
	$O/detect_surface.shared.o \
	$O/exit_mpi.shared.o \
	$O/force_ftz.cc.o \
	$O/get_attenuation_model.shared.o \
	$O/get_element_face.shared.o \
	$O/get_jacobian_boundaries.shared.o \
	$O/get_shape3D.shared.o \
	$O/gll_library.shared.o \
	$O/heap_sort.shared.o \
	$O/hex_nodes.shared.o \
	$O/lagrange_poly.shared.o \
	$O/netlib_specfun_erf.shared.o \
	$O/param_reader.cc.o \
	$O/prepare_assemble_MPI.shared.o \
	$O/read_topo_bathy_file.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/save_header_file.shared.o \
	$O/search_kdtree.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/utm_geo.shared.o \
	$O/write_VTK_data.shared.o \
	$O/write_c_binary.cc.o \
	$(EMPTY_MACRO)

fwat_meshfem3D_OBJECTS = \
	$O/calc_gll_points.mesh.o \
	$O/create_interfaces_mesh.mesh.o \
	$O/compute_parameters.mesh.o \
	$O/create_meshfem_mesh.mesh.o \
	$O/chunk_earth_mesh_mod.mesh.o \
	$O/earth_chunk.mesh.o \
	$O/read_value_mesh_parameters.mesh.o \
	$O/determine_cavity.mesh.o \
	$O/create_CPML_regions.mesh.o \
	$O/create_visual_files.mesh.o \
	$O/store_coords.mesh.o \
	$O/store_boundaries.mesh.o \
	$O/check_mesh_quality.mesh.o \
	$O/save_databases.mesh.o \
	$O/meshfem3D_adios_stubs.mesh_noadios.o \
	$O/define_subregions.mesh.o \
	$O/define_subregions_heuristic.mesh.o \
	$O/define_superbrick.mesh.o \
	$O/get_flags_boundaries.mesh.o \
	$O/get_global.shared.o \
	$O/read_mesh_parameter_file.mesh.o \
	$(EMPTY_MACRO)

fwat_generate_database_OBJECTS = \
	$O/generate_databases_par.gen_mod.o \
	$O/calc_jacobian.gen.o \
	$O/fault_generate_databases.gen.o \
	$O/create_mass_matrices.gen.o \
	$O/create_regions_mesh.gen.o \
	$O/finalize_databases.gen.o \
	$O/get_absorbing_boundary.gen.o \
	$O/get_coupling_surfaces.gen.o \
	$O/get_model.gen.o \
	$O/get_MPI.gen.o \
	$O/get_perm_color.gen.o \
	$O/model_1d_cascadia.gen.o \
	$O/model_1d_prem.gen.o \
	$O/model_1d_socal.gen.o \
	$O/model_aniso.gen.o \
	$O/model_coupled.gen.o \
	$O/model_default.gen.o \
	$O/model_external_values.gen.o \
	$O/model_ipati.gen.o \
	$O/parse_sep.genc.o \
	$O/model_gll.gen.o \
	$O/model_salton_trough.gen.o \
	$O/model_tomography.gen.o \
	$O/pml_set_local_dampingcoeff.gen.o \
	$O/read_parameters.gen.o \
	$O/read_partition_files.gen.o \
	$O/save_arrays_solver.gen.o \
	$O/setup_color_perm.gen.o \
	$O/setup_mesh.gen.o \
	$O/memory_eval.gen.o \
	$O/model_sep.mpi_gen.o \
	$O/get_shape2D.shared.o \
	$(EMPTY_MACRO)


fullwave_adjoint_tomo_MODULES = \
	$(FC_MODDIR)/fullwave_adjoint_tomo_par.$(FC_MODEXT) \
	$(FC_MODDIR)/fwat_input.$(FC_MODEXT) \
	$(FC_MODDIR)/fwat_utils.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_interface.$(FC_MODEXT) \
	$(FC_MODDIR)/interpolation_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/telestf_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/ascii_rw.$(FC_MODEXT) \
	$(FC_MODDIR)/ma_sub.$(FC_MODEXT) \
	$(FC_MODDIR)/ma_sub2.$(FC_MODEXT) \
	$(FC_MODDIR)/measure_adj.$(FC_MODEXT) \
	$(FC_MODDIR)/taper_3D_mod.$(FC_MODEXT) \
	$(EMPTY_MACRO)


###
### MPI
###
fullwave_adjoint_tomo_SHARED_OBJECTS += $(COND_MPI_OBJECTS)

###
### OPENMP
###
fullwave_adjoint_tomo_SHARED_OBJECTS += $(COND_OPENMP_OBJECTS)

###
### CUDA
###
cuda_fullwave_adjoint_tomo_OBJECTS = \
	$O/assemble_MPI_scalar_cuda.cuda.o \
	$O/assemble_MPI_vector_cuda.cuda.o \
	$O/check_fields_cuda.cuda.o \
	$O/compute_add_sources_acoustic_cuda.cuda.o \
	$O/compute_add_sources_viscoelastic_cuda.cuda.o \
	$O/compute_coupling_cuda.cuda.o \
	$O/compute_forces_acoustic_cuda.cuda.o \
	$O/compute_forces_viscoelastic_cuda.cuda.o \
	$O/compute_kernels_cuda.cuda.o \
	$O/compute_stacey_acoustic_cuda.cuda.o \
	$O/compute_stacey_viscoelastic_cuda.cuda.o \
	$O/initialize_cuda.cuda.o \
	$O/noise_tomography_cuda.cuda.o \
	$O/prepare_mesh_constants_cuda.cuda.o \
	$O/save_and_compare_cpu_vs_gpu.cudacc.o \
	$O/transfer_fields_cuda.cuda.o \
	$O/update_displacement_cuda.cuda.o \
	$O/write_seismograms_cuda.cuda.o \
	$O/fault_solver_dynamics.cuda.o \
	$(EMPTY_MACRO)

cuda_fullwave_adjoint_tomo_STUBS = \
	$O/specfem3D_gpu_cuda_method_stubs.cudacc.o \
	$(EMPTY_MACRO)

cuda_fullwave_adjoint_tomo_DEVICE_OBJ = \
	$O/cuda_device_obj.o \
	$(EMPTY_MACRO)

ifeq ($(CUDA),yes)
fullwave_adjoint_tomo_OBJECTS += $(cuda_fullwave_adjoint_tomo_OBJECTS)
ifeq ($(CUDA_PLUS),yes)
fullwave_adjoint_tomo_OBJECTS += $(cuda_fullwave_adjoint_tomo_DEVICE_OBJ)
endif
else
fullwave_adjoint_tomo_OBJECTS += $(cuda_fullwave_adjoint_tomo_STUBS)
endif

###
### ADIOS
###

# using ADIOS files

adios_fullwave_adjoint_tomo_OBJECTS= \
	$O/read_mesh_databases_adios.spec_adios.o \
	$O/save_forward_arrays_adios.spec_adios.o \
	$O/read_forward_arrays_adios.spec_adios.o \
	$O/save_kernels_adios.spec_adios.o \
	$O/read_partition_files_adios.gen_adios.o \
	$O/save_arrays_solver_adios.gen_adios.o \
	$O/save_moho_adios.gen_adios.o \
	$O/model_gll_adios.gen_adios.o \
	$O/model_ipati_adios.gen_adios.o \
	$(EMPTY_MACRO)

adios_fullwave_adjoint_tomo_PREOBJECTS = \
	$O/adios_manager.shared_adios.o \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o \
	$(EMPTY_MACRO)


adios_fullwave_adjoint_tomo_STUBS = \
	$O/specfem3D_adios_stubs.spec_noadios.o \
	$O/generate_databases_adios_stubs.gen_noadios.o

adios_fullwave_adjoint_tomo_PRESTUBS = \
	$O/adios_manager_stubs.shared_noadios.o

# conditional adios linking
ifeq ($(ADIOS),no)
adios_fullwave_adjoint_tomo_OBJECTS = $(adios_fullwave_adjoint_tomo_STUBS)
adios_fullwave_adjoint_tomo_PREOBJECTS = $(adios_fullwave_adjoint_tomo_PRESTUBS)
endif
fullwave_adjoint_tomo_OBJECTS += $(adios_fullwave_adjoint_tomo_OBJECTS)
fullwave_adjoint_tomo_SHARED_OBJECTS += $(adios_fullwave_adjoint_tomo_PREOBJECTS)

###
### ASDF
###

asdf_fullwave_adjoint_tomo_OBJECTS = \
        $O/write_output_ASDF.spec.o \
        $O/read_adjoint_sources_ASDF.spec.o \
        $(EMPTY_MACRO)

asdf_fullwave_adjoint_tomo_SHARED_OBJECTS = \
        $O/asdf_manager.shared_asdf.o \
        $(EMPTY_MACRO)

asdf_fullwave_adjoint_tomo_SHARED_STUBS = \
        $O/asdf_method_stubs.cc.o \
        $O/asdf_manager_stubs.shared_asdf.o \
        $(EMPTY_MACRO)

# conditional asdf linking
ifeq ($(ASDF),yes)
SPECFEM_LINK_FLAGS += $(ASDF_LIBS) -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lstdc++
fullwave_adjoint_tomo_OBJECTS += $(asdf_fullwave_adjoint_tomo_OBJECTS)
fullwave_adjoint_tomo_SHARED_OBJECTS += $(asdf_fullwave_adjoint_tomo_SHARED_OBJECTS)
else
fullwave_adjoint_tomo_OBJECTS += $(asdf_fullwave_adjoint_tomo_STUBS)
fullwave_adjoint_tomo_SHARED_OBJECTS += $(asdf_fullwave_adjoint_tomo_SHARED_STUBS)
endif
#

###
### VTK
###

ifeq ($(VTK),yes)
fullwave_adjoint_tomo_OBJECTS += \
	$O/vtk_window.spec.o \
	$O/vtk_helper.visualcc.o \
	$(EMPTY_MACRO)
fullwave_adjoint_tomo_MODULES += \
	$(FC_MODDIR)/vtk_window_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)
endif


##
## xfullwave_adjoint_tomo
##
xfullwave_adjoint_tomo_OBJECTS = \
	$O/specfem_meshfem3D.fwat_sr.o \
	$(fwat_meshfem3D_OBJECTS) \
	$(fwat_generate_database_OBJECTS) \
	$(fullwave_adjoint_tomo_OBJECTS) \
	$O/fullwave_adjoint_tomo_main.fwat_main.o \
	$O/program_fullwave_adjoint_tomo.fwat.o \
	$(EMPTY_MACRO)
xfwat1_fwd_measure_adj_OBJECTS = \
	$O/program_fwat1_fwd_measure_adj.fwat_stage1.o \
	$(fullwave_adjoint_tomo_OBJECTS) \
	$(EMPTY_MACRO)
xfwat3_linesearch_OBJECTS = \
	$O/run_fwat3_linesearch.fwat_ls.o \
	$O/program_fwat3_linesearch.fwat_stage3.o \
	$(fullwave_adjoint_tomo_OBJECTS) \
	$(EMPTY_MACRO)
xfwat0_forward_data_OBJECTS = \
	$O/run_fwat0_forward_data.fwat_fwd.o \
	$O/program_fwat0_forward_data.fwat_stage0.o \
	$(fullwave_adjoint_tomo_OBJECTS) \
	$(EMPTY_MACRO)
xfwat4_cmt3d_fwd_OBJECTS = \
	$O/run_fwat4_cmt3d_fwd.fwat_fwd.o \
	$O/program_fwat4_cmt3d_fwd.fwat_stage4.o \
	$(fullwave_adjoint_tomo_OBJECTS) \
	$(EMPTY_MACRO)
xfwat5_rtm_OBJECTS = \
	$O/tomography_par.opt_module.o \
	$O/taper_3D_mod.fwat_postproc.o \
	$O/run_postprocessing.fwat_postproc.o \
	$O/postprocess_par.fwat_postproc_module.o \
	$O/parse_kernel_names.fwat_postproc.o \
	$O/run_fwat5_RTM.fwat_rtm.o \
	$O/program_fwat5_RTM.fwat_stage5.o \
	$(fullwave_adjoint_tomo_OBJECTS) \
	$(EMPTY_MACRO)

xfwat_mesh_database_OBJECTS = \
	$O/specfem_meshfem3D.fwat_sr.o \
	$O/program_fwat_mesh_database.fwat_mesh_database.o \
	$(fwat_meshfem3D_OBJECTS) \
	$(fwat_generate_database_OBJECTS) \
	$(adios_fullwave_adjoint_tomo_OBJECTS) \
	$(fullwave_adjoint_tomo_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

##
## xfwat2_postproc_opt
##
xfwat2_postproc_opt_top_OBJECTS = \
	$O/tomography_par.opt_module.o \
	$O/postprocess_par.fwat_postproc_module.o \
	$O/parse_kernel_names.fwat_postproc.o \
	$O/taper_3D_mod.fwat_postproc.o \
	$O/run_postprocessing.fwat_postproc.o \
	$O/postproc_subs.fwat_postproc.o \
	$O/run_optimization.fwat_opt.o \
	$O/compute_kernel_integral.opt.o \
	$O/get_cg_direction.opt.o \
	$O/get_lbfgs_direction.opt.o \
	$O/get_sd_direction.opt.o \
	$O/read_kernels.opt.o \
	$O/read_kernels_cg.opt.o \
	$O/read_model.opt.o \
	$O/read_parameters_tomo.opt.o \
	$O/read_parameters_lbfgs.opt.o \
	$O/write_gradient.opt.o \
	$O/write_new_model.opt.o \
	$O/write_new_model_perturbations.opt.o \
	$O/run_fwat2_postproc_opt.fwat_postopt.o \
	$(EMPTY_MACRO)
xfwat2_postproc_opt_OBJECTS = \
	$O/program_fwat2_postproc_opt.fwat_stage2.o \
	$O/interpolation_mod.fwat_interp.o \
	$(xfwat2_postproc_opt_top_OBJECTS) \
	$(EMPTY_MACRO)
xfwat2_postproc_opt_SHARED_OBJECTS = \
	$O/fullwave_adjoint_tomo_par.fwat_par.o \
	$O/fwat_input.fwat_readpar.o \
	$O/fwat_utils.fwat_utils.o \
	$(xsmooth_sem_sph_pde_SHARED_OBJECTS) \
	$(EMPTY_MACRO)
xdynmbat_sele_ctrlgrp_top_OBJECTS = \
	$O/tomography_par.opt_module.o \
	$O/postprocess_par.fwat_postproc_module.o \
	$O/parse_kernel_names.fwat_postproc.o \
	$O/interpolation_mod.fwat_interp.o \
	$O/taper_3D_mod.fwat_postproc.o \
	$O/run_postprocessing.fwat_postproc.o \
	$O/run_optimization.fwat_opt.o \
	$O/compute_kernel_integral.opt.o \
	$O/get_cg_direction.opt.o \
	$O/get_lbfgs_direction.opt.o \
	$O/get_sd_direction.opt.o \
	$O/read_kernels.opt.o \
	$O/read_kernels_cg.opt.o \
	$O/read_model.opt.o \
	$O/read_parameters_tomo.opt.o \
	$O/read_parameters_lbfgs.opt.o \
	$O/write_gradient.opt.o \
	$O/write_new_model.opt.o \
	$O/write_new_model_perturbations.opt.o \
	$O/run_sele_ctrlgrp.dynmbat.o \
	$(EMPTY_MACRO)
xdynmbat_sele_ctrlgrp_OBJECTS = \
	$O/program_dynmbat_sele_ctrlgrp.dynmbat.o \
	$(xdynmbat_sele_ctrlgrp_top_OBJECTS) \
	$(EMPTY_MACRO)
xdynmbat_sele_ctrlgrp_SHARED_OBJECTS = \
	$O/fullwave_adjoint_tomo_par.fwat_par.o \
	$O/fwat_input.fwat_readpar.o \
	$(xsmooth_sem_sph_pde_SHARED_OBJECTS) \
	$(EMPTY_MACRO)
grid_and_combine_vol_data_OBJECTS = \
	$O/fullwave_adjoint_tomo_par.fwat_par.o \
	$O/npy.fwat_aux.o \
	$O/specfem3D_par.spec_module.o \
	$O/projection_on_FD_grid_mod_fwat.fwat_aux.o \
	$O/model_grid.fwat_aux.o \
	$(EMPTY_MACRO)
grid_and_combine_vol_data_cert_OBJECTS = \
	$O/fullwave_adjoint_tomo_par.fwat_par.o \
	$O/npy.fwat_aux.o \
	$O/specfem3D_par.spec_module.o \
	$O/projection_on_FD_grid_mod_fwat.fwat_aux.o \
	$O/model_grid_cert.fwat_aux.o \
	$(EMPTY_MACRO)
decomp_grid2gll_OBJECTS = \
	$O/fullwave_adjoint_tomo_par.fwat_par.o \
	$O/npy.fwat_aux.o \
	$O/specfem3D_par.spec_module.o \
	$O/projection_on_FD_grid_mod_fwat.fwat_aux.o \
	$O/model_gll.fwat_aux.o \
	$(EMPTY_MACRO)
add_pert_gaus_OBJECTS = \
	$O/fullwave_adjoint_tomo_par.fwat_par.o \
	$O/specfem3D_par.spec_module.o \
	$O/add_pert_gaus.fwat_aux.o \
	$(EMPTY_MACRO)
add_pert_trigo_OBJECTS = \
	$O/interpolation_mod.fwat_interp.o \
	$O/specfem3D_par.spec_module.o \
	$O/add_pert_trigo.fwat_aux.o \
	$(EMPTY_MACRO)

	
ifeq ($(CUDA),yes)
## cuda version
xfwat2_postproc_opt_OBJECTS += $(cuda_smooth_sem_OBJECTS)
xdynmbat_sele_ctrlgrp_OBJECTS += $(cuda_smooth_sem_OBJECTS)
ifeq ($(CUDA_PLUS),yes)
xfwat2_postproc_opt_OBJECTS += $(cuda_smooth_sem_DEVICE_OBJ)
xdynmbat_sele_ctrlgrp_OBJECTS += $(cuda_smooth_sem_DEVICE_OBJ)
endif
## libs
#xsmooth_sem_LIBS = $(MPILIBS) $(CUDA_LINK)
#INFO_CUDA_SEM="building xsmooth_sem with CUDA support"
else
## non-cuda version
xfwat2_postproc_opt_OBJECTS += $(cuda_smooth_sem_STUBS)
xdynmbat_sele_ctrlgrp_OBJECTS += $(cuda_smooth_sem_STUBS)
## libs
#xsmooth_sem_LIBS = $(MPILIBS)
#INFO_CUDA_SEM="building xsmooth_sem without CUDA support"
endif

##
## xmeasure_adj
##
xmeasure_adj_OBJECTS = \
	$(preproc_measure_adj_OBJECTS) \
	$O/fftpack.fwat_deconit.o \
	$O/deconit.fwat_deconit.o \
	$O/interpolation_mod.fwat_interp.o \
	$O/fullwave_adjoint_tomo_par.fwat_par.o \
	$O/exit_mpi.shared.o \
	$O/parallel.sharedmpi.o  \
	$O/param_reader.cc.o \
	$O/read_value_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/specfem3D_par.spec_module.o \
	$O/shared_par.shared_module.o \
	$O/program_measure_adj.fwat_xmeas.o \
	$(EMPTY_MACRO)

##
## xflexwin
##
xflexwin_OBJECTS = \
	$(preproc_flexwin_OBJECTS) \
	$O/xapiir_sub.fwat_preproc_meas.o \
	$O/sacio.fwat_preproc_meas.o \
	$O/fullwave_adjoint_tomo_par.fwat_par.o \
	$O/exit_mpi.shared.o \
	$O/parallel.sharedmpi.o \
	$O/read_value_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/param_reader.cc.o \
	$O/specfem3D_par.spec_module.o \
	$O/shared_par.shared_module.o \
	$O/program_flexwin.fwat_xflex.o \
	$(EMPTY_MACRO)

##
## xdynmbat_sele_minibatch
##
xdynmbat_sele_minibatch_OBJECTS = \
	$O/program_dynmbat_sele_minibatch.dynmbat.o \
	$O/distaz.fwat_preproc_flex.o \
	$(EMPTY_MACRO)

#######################################

####
#### rules for executables
####


ifeq ($(CUDA),yes)
## cuda version
ifeq ($(CUDA_PLUS),yes)
## cuda 5x & 6x version
INFO_CUDA_INVERSE_PROBLEM="building xfullwave_adjoint_tomo with CUDA support"
else
## cuda 4 version
INFO_CUDA_INVERSE_PROBLEM="building xfullwave_adjoint_tomo with CUDA 4 support"
endif

${E}/xfullwave_adjoint_tomo: $(xfullwave_adjoint_tomo_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA_INVERSE_PROBLEM)
	@echo ""
	${FCLINK} -o ${E}/xfullwave_adjoint_tomo $(xfullwave_adjoint_tomo_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS) $(CUDA_LINK)
	@echo ""
${E}/xfwat1_fwd_measure_adj: $(xfwat1_fwd_measure_adj_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA_INVERSE_PROBLEM)
	@echo ""
	${FCLINK} -o ${E}/xfwat1_fwd_measure_adj $(xfwat1_fwd_measure_adj_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS) $(CUDA_LINK)
	@echo ""
${E}/xfwat3_linesearch: $(xfwat3_linesearch_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA_INVERSE_PROBLEM)
	@echo ""
	${FCLINK} -o ${E}/xfwat3_linesearch $(xfwat3_linesearch_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS) $(CUDA_LINK)
	@echo ""
${E}/xfwat0_forward_data: $(xfwat0_forward_data_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA_INVERSE_PROBLEM)
	@echo ""
	${FCLINK} -o ${E}/xfwat0_forward_data $(xfwat0_forward_data_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS) $(CUDA_LINK)
	@echo ""
${E}/xfwat4_cmt3d_fwd: $(xfwat4_cmt3d_fwd_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA_INVERSE_PROBLEM)
	@echo ""
	${FCLINK} -o ${E}/xfwat4_cmt3d_fwd $(xfwat4_cmt3d_fwd_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS) $(CUDA_LINK)
	@echo ""
${E}/xfwat5_rtm: $(xfwat5_rtm_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo "building xfwat5_rtm"
	@echo ""
	${FCLINK} -o ${E}/xfwat5_rtm $(xfwat5_rtm_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS) $(CUDA_LINK)
	@echo ""
else

## non-cuda version
${E}/xfullwave_adjoint_tomo: $(xfwat2_postproc_opt_top_OBJECTS)  $(xfullwave_adjoint_tomo_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) 
	@echo ""
	@echo "building xfullwave_adjoint_tomo"
	@echo ""
	${FCLINK} -o ${E}/xfullwave_adjoint_tomo  $(xfwat2_postproc_opt_top_OBJECTS)  $(xfullwave_adjoint_tomo_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""
${E}/xfwat1_fwd_measure_adj: $(xfwat1_fwd_measure_adj_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo "building xfwat1_fwd_measure_adj"
	@echo ""
	${FCLINK} -o ${E}/xfwat1_fwd_measure_adj $(xfwat1_fwd_measure_adj_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""
${E}/xfwat3_linesearch: $(xfwat3_linesearch_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo "building xfwat3_linesearch"
	@echo ""
	${FCLINK} -o ${E}/xfwat3_linesearch $(xfwat3_linesearch_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""
${E}/xfwat0_forward_data: $(xfwat0_forward_data_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo "building xfwat0_forward_data"
	@echo ""
	${FCLINK} -o ${E}/xfwat0_forward_data $(xfwat0_forward_data_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""
${E}/xfwat4_cmt3d_fwd: $(xfwat4_cmt3d_fwd_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo "building xfwat4_cmt3d_fwd"
	@echo ""
	${FCLINK} -o ${E}/xfwat4_cmt3d_fwd $(xfwat4_cmt3d_fwd_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""

${E}/xfwat5_rtm: $(xfwat5_rtm_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS)
	@echo ""
	@echo "building xfwat5_rtm"
	@echo ""
	${FCLINK} -o ${E}/xfwat5_rtm $(xfwat5_rtm_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""
${E}/xgrid_and_combine_vol_data: $(grid_and_combine_vol_data_OBJECTS)
	@echo ""
	@echo "building xgrid_and_combine_vol_data"
	@echo ""
	${FCLINK} -o ${E}/xgrid_and_combine_vol_data $(grid_and_combine_vol_data_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""
${E}/xgrid_and_combine_vol_data_cert: $(grid_and_combine_vol_data_cert_OBJECTS)
	@echo ""
	@echo "building xgrid_and_combine_vol_data_cert"
	@echo ""
	${FCLINK} -o ${E}/xgrid_and_combine_vol_data_cert $(grid_and_combine_vol_data_cert_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""
${E}/xdecomp_grid2gll: $(decomp_grid2gll_OBJECTS)
	@echo ""
	@echo "building xdecomp_grid2gll"
	@echo ""
	${FCLINK} -o ${E}/xdecomp_grid2gll $(decomp_grid2gll_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""
${E}/xadd_pert_gaus: $(add_pert_gaus_OBJECTS)
	@echo ""
	@echo "building xadd_pert_gaus"
	@echo ""
	${FCLINK} -o ${E}/xadd_pert_gaus $(add_pert_gaus_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""
${E}/xadd_pert_trigo: $(add_pert_trigo_OBJECTS)
	@echo ""
	@echo "building xadd_pert_trigo"
	@echo ""
	${FCLINK} -o ${E}/xadd_pert_trigo $(add_pert_trigo_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(MPILIBS)
	@echo ""
endif

${E}/xfwat2_postproc_opt: $(xfwat2_postproc_opt_OBJECTS) $(xfwat2_postproc_opt_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xfwat2_postproc_opt"
	@echo ""
	${FCLINK} -o $@ $(xfwat2_postproc_opt_OBJECTS) $(xfwat2_postproc_opt_SHARED_OBJECTS) $(COND_MPI_OBJECTS) $(xsmooth_sem_LIBS)
	@echo ""

${E}/xdynmbat_sele_ctrlgrp: $(xdynmbat_sele_ctrlgrp_OBJECTS) $(xdynmbat_sele_ctrlgrp_SHARED_OBJECTS)  $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xdynmbat_sele_ctrlgrp"
	@echo ""
	${FCLINK} -o $@ $(xdynmbat_sele_ctrlgrp_OBJECTS) $(xdynmbat_sele_ctrlgrp_SHARED_OBJECTS)  $(COND_MPI_OBJECTS) $(xsmooth_sem_LIBS)
	@echo ""

$E/xmeasure_adj: $(xmeasure_adj_OBJECTS)
	@echo ""
	@echo "building xmeasure_adj"
	@echo ""
	${FCLINK} -o ${E}/xmeasure_adj $(xmeasure_adj_OBJECTS)  $(MPILIBS)
	@echo ""

$E/xflexwin: $(xflexwin_OBJECTS)
	@echo ""
	@echo "building xflexwin"
	@echo ""
	${FCLINK} -o ${E}/xflexwin $(xflexwin_OBJECTS)  $(MPILIBS)
	@echo ""

$E/xdynmbat_sele_minibatch: $(xdynmbat_sele_minibatch_OBJECTS)
	@echo ""
	@echo "building xdynmbat_sele_minibatch"
	@echo ""
	${FCLINK} -o ${E}/xdynmbat_sele_minibatch $(xdynmbat_sele_minibatch_OBJECTS) 
	@echo ""

$E/xfwat_mesh_database: $(xfwat_mesh_database_OBJECTS) $(fullwave_adjoint_tomo_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xfwat_mesh_database"
	@echo ""
	${FCLINK} -o ${E}/xfwat_mesh_database $(xfwat_mesh_database_OBJECTS) $(MPILIBS)
	@echo ""
#######################################

###
### Module dependencies
###

$O/fullwave_adjoint_tomo_main.fwat.o: \
	$O/run_fwat1_fwd_measure_adj.fwat_fwdadj.o

$O/run_fwat1_fwd_measure_adj.fwat_fwdadj.o: \
	$O/specfem_initialize_simulation_fwat.fwat_spec.o \
	$O/specfem_interface.fwat_spec.o \
	$O/preproc_collect_data.fwat_preproc.o \
	$O/run_preprocessing.fwat_preproc.o \
	$O/run_postprocessing.fwat_postproc.o \
	$O/fwat_input.fwat_readpar.o
$O/program_dynmbat_sele_ctrlgrp.dynmbat.o: \
	$O/run_sele_ctrlgrp.dynmbat.o	

$O/preproc_collect_data.fwat_preproc.o: \
	$(preproc_measure_adj_OBJECTS)
$O/preproc_subs.fwat_preproc.o : \
	$O/spanlib.fwat_pca.o \
	$O/telestf_mod.fwat_telestf.o

$O/run_preprocessing.fwat_preproc.o: \
	$(preproc_flexwin_OBJECTS) \
	$O/preproc_collect_data.fwat_preproc.o \
	$O/preproc_subs.fwat_preproc.o  \
	$O/preproc_measure_adj_subs.fwat_preproc.o \
	$O/interpolation_mod.fwat_interp.o \
	$O/spanlib.fwat_pca.o
$O/measure_adj.fwat_preproc_meas.o: \
	$O/fftpack.fwat_deconit.o \
	$O/deconit.fwat_deconit.o \
	$O/interpolation_mod.fwat_interp.o \
	$O/xapiir_sub.fwat_preproc_meas.o \
	$O/sacio.fwat_preproc_meas.o \
	$O/ma_sub.fwat_preproc_meas.o \
	$O/ma_sub2.fwat_preproc_meas.o \
	$O/ascii_rw.fwat_preproc_meas.o
$O/ma_sub.fwat_preproc_meas.o: \
	$O/xapiir_sub.fwat_preproc_meas.o \
	$O/sacio.fwat_preproc_meas.o \
	$O/ma_sub2.fwat_preproc_meas.o \
	$O/ascii_rw.fwat_preproc_meas.o

$O/taper_3D_mod.fwat_postproc.o:\
	$O/interpolation_mod.fwat_interp.o 

$O/run_postprocessing.fwat_postproc.o: \
	$O/taper_3D_mod.fwat_postproc.o

####
#### rule to build each .o file
####

## main module
$O/%.fwat_par.o: $S/%.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o $O/specfem3D_par.spec_module.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

## file object rules
$O/%.fwat.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_main.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_stage1.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_stage2.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.dynmbat.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_stage3.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_stage0.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_stage4.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_stage5.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_mesh_database.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_xmeas.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_utils.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_xflex.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_rtm.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_postopt.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_fwdadj.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_ls.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_fwd.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_readpar.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_preproc.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_precond.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_semd2sac.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_postproc.o: $S/%.F90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_opt.o: $S/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_preproc_meas.o: $S/preproc_measure_adj/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o  $@ $<
$O/%.fwat_preproc_flex.o: $S/preproc_flexwin/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o  $@ $<
$O/%.fwat_preproc_flex.o: $S/preproc_flexwin/%.f
	${FC} ${FCFLAGS} -c -o  $@ $<
$O/%.fwat_preproc_flex.o: $S/preproc_flexwin/%.c
	${CC} ${CFLAGS} -lm -I./src/fullwave_adjoint_tomo/sac_inc -c -o $@ $<
$O/%.fwat_savekern.o: $S/specfem_interface/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_sr.o: $S/specfem_interface/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_spec.o: $S/specfem_interface/%.F90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_interp.o: $S/teleseis_stf/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_pca.o: $S/teleseis_stf/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_telestf.o: $S/teleseis_stf/%.f90 ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_deconit.o: $S/preproc_recfunc/%.f90 
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
$O/%.fwat_aux.o: $S/auxiliaries/%.f90  ${SETUP}/constants.h $O/fullwave_adjoint_tomo_par.fwat_par.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<
