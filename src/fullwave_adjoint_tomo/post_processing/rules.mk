#=====================================================================
#
#               S p e c f e m 3 D  V e r s i o n  3 . 0
#               ---------------------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#                and CNRS / University of Marseille, France
#                 (there are currently many more authors!)
# (c) Princeton University and CNRS / University of Marseille, July 2012
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================
# WK change directory: fullwave_adjoint_tomo/post_processing
# WK change target: fwat_postproc
## compilation directories
S := ${S_TOP}/src/fullwave_adjoint_tomo/post_processing
$(fullwave_adjoint_tomo/post_processing_OBJECTS): S := ${S_TOP}/src/fullwave_adjoint_tomo/post_processing

#######################################

fullwave_adjoint_tomo/post_processing_TARGETS = \
	$E/xclip_sem \
	$E/xcombine_sem \
	$E/xsmooth_sem \
	$E/xsmooth_sem_glob \
	$E/xsmooth_sem_sph_pde \
	$E/xsmooth_sem_pde \
	$(EMPTY_MACRO)

fullwave_adjoint_tomo/post_processing_OBJECTS = \
	$(xclip_sem_OBJECTS) \
	$(xcombine_sem_OBJECTS) \
	$(xsmooth_sem_OBJECTS) \
	$(xsmooth_sem_sph_pde_OBJECTS) \
	$(xsmooth_sem_pde_OBJECTS) \
	$(EMPTY_MACRO)

fullwave_adjoint_tomo/post_processing_SHARED_OBJECTS = \
	$(xclip_sem_SHARED_OBJECTS) \
	$(xcombine_sem_SHARED_OBJECTS) \
	$(xsmooth_sem_SHARED_OBJECTS) \
	$(xsmooth_sem_sph_pde_SHARED_OBJECTS) \
	$(xsmooth_sem_pde_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

fullwave_adjoint_tomo/post_processing_MODULES = \
	$(FC_MODDIR)/postprocess_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

####
#### rules for executables
####

.PHONY: fwat_postproc


fwat_postproc: $(fullwave_adjoint_tomo/post_processing_TARGETS)

post_processing: fwat_postproc

fullwave_adjoint_tomo/post_processing: fwat_postproc

### single targets

clip_sem: xclip_sem
xclip_sem: $E/xclip_sem

combine_sem: xcombine_sem
xcombine_sem: $E/xcombine_sem

smooth_sem: xsmooth_sem
xsmooth_sem: $E/xsmooth_sem

smooth_sem_sph_pde: xsmooth_sem_sph_pde
xsmooth_sem_sph_pde: $E/xsmooth_sem_sph_pde

smooth_sem_pde: xsmooth_sem_pde
xsmooth_sem_pde: $E/xsmooth_sem_pde

#######################################

####
#### rules for each program follow
####

#######################################


##
## xclip_sem
##
xclip_sem_OBJECTS = \
	$O/postprocess_par.fwat_postproc_module.o \
	$O/parse_kernel_names.fwat_postproc.o \
	$O/clip_sem.fwat_postproc.o \
	$(EMPTY_MACRO)

xclip_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

${E}/xclip_sem: $(xclip_sem_OBJECTS) $(xclip_sem_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xclip_sem"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""


##
## xcombine_sem
##
xcombine_sem_OBJECTS = \
	$O/postprocess_par.fwat_postproc_module.o \
	$O/parse_kernel_names.fwat_postproc.o \
	$O/combine_sem.fwat_postproc.o \
	$(EMPTY_MACRO)

xcombine_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

${E}/xcombine_sem: $(xcombine_sem_OBJECTS) $(xcombine_sem_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xcombine_sem"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""


##
## xsmooth_sem
##
xsmooth_sem_OBJECTS = \
	$O/postprocess_par.fwat_postproc_module.o \
	$O/parse_kernel_names.fwat_postproc.o \
	$O/smooth_sem.fwat_postproc.o \
	$(EMPTY_MACRO)

xsmooth_sem_glob_OBJECTS = \
	$O/postprocess_par.fwat_postproc_module.o \
	$O/parse_kernel_names.fwat_postproc.o \
	$O/smooth_sem_glob.fwat_postproc.o \
	$O/smooth_weights_vec.fwat_postproc.o \
	$(EMPTY_MACRO)

xsmooth_sem_sph_pde_OBJECTS = \
	$O/postprocess_par.fwat_postproc_module.o \
	$O/parse_kernel_names.fwat_postproc.o \
	$O/smooth_sem_sph_pde.fwat_postproc.o \
	$(EMPTY_MACRO)

xsmooth_sem_pde_OBJECTS = \
	$O/postprocess_par.fwat_postproc_module.o \
	$O/parse_kernel_names.fwat_postproc.o \
	$O/smooth_sem_pde.fwat_postproc.o \
	$(EMPTY_MACRO)

xsmooth_sem_SHARED_OBJECTS = \
	$O/specfem3D_par.spec_module.o \
	$O/pml_par.spec.o \
	$O/read_mesh_databases.spec.o \
	$O/shared_par.shared_module.o \
	$O/check_mesh_resolution.shared.o \
	$O/create_name_database.shared.o \
	$O/exit_mpi.shared.o \
	$O/gll_library.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/write_VTK_data.shared.o \
	$(EMPTY_MACRO)

xsmooth_sem_sph_pde_SHARED_OBJECTS = \
	$(xsmooth_sem_SHARED_OBJECTS) \
	$O/assemble_MPI_scalar.shared.o \
	$O/define_derivation_matrices.shared.o \
	$O/detect_surface.shared.o \
	$O/force_ftz.cc.o \
	$O/get_attenuation_model.shared.o \
	$O/get_element_face.shared.o \
	$O/get_jacobian_boundaries.shared.o \
	$O/get_shape3D.shared.o \
	$O/heap_sort.shared.o \
	$O/hex_nodes.shared.o \
	$O/lagrange_poly.shared.o \
	$O/netlib_specfun_erf.shared.o \
	$O/prepare_assemble_MPI.shared.o \
	$O/read_topo_bathy_file.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/save_header_file.shared.o \
	$O/search_kdtree.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/utm_geo.shared.o \
	$O/write_c_binary.cc.o

xsmooth_sem_pde_SHARED_OBJECTS = \
	$(xsmooth_sem_sph_pde_SHARED_OBJECTS)

cuda_smooth_sem_STUBS = \
	$O/smooth_sem_cuda_stubs.fwat_postproc.o \
	$(EMPTY_MACRO)

cuda_smooth_sem_OBJECTS = \
	$O/smooth_cuda.fwat_postproc.cuda.o \
	$O/check_fields_cuda.cuda.o \
	$O/initialize_cuda.cuda.o \
	$(EMPTY_MACRO)

cuda_smooth_sem_DEVICE_OBJ = \
	$O/cuda_device_smooth_obj.o \
	$(EMPTY_MACRO)

ifeq ($(CUDA),yes)
## cuda version
xsmooth_sem_OBJECTS += $(cuda_smooth_sem_OBJECTS)
ifeq ($(CUDA_PLUS),yes)
xsmooth_sem_OBJECTS += $(cuda_smooth_sem_DEVICE_OBJ)
endif
## libs
xsmooth_sem_LIBS = $(MPILIBS) $(CUDA_LINK)
INFO_CUDA_SEM="building xsmooth_sem with CUDA support"
else
## non-cuda version
xsmooth_sem_OBJECTS += $(cuda_smooth_sem_STUBS)
xsmooth_sem_glob_OBJECTS += $(cuda_smooth_sem_STUBS)
xsmooth_sem_sph_pde_OBJECTS += $(cuda_smooth_sem_STUBS)
xsmooth_sem_pde_OBJECTS += $(cuda_smooth_sem_STUBS)
## libs
xsmooth_sem_LIBS = $(MPILIBS)
INFO_CUDA_SEM="building xsmooth_sem without CUDA support"
endif

# extra dependencies
$O/smooth_sem.fwat_postproc.o: $O/specfem3D_par.spec_module.o $O/postprocess_par.fwat_postproc_module.o
$O/smooth_sem_sph_pde.fwat_postproc.o: $O/specfem3D_par.spec_module.o $O/postprocess_par.fwat_postproc_module.o

${E}/xsmooth_sem: $(xsmooth_sem_OBJECTS) $(xsmooth_sem_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA_SEM)
	@echo ""
	${FCLINK} -o $@ $+ $(xsmooth_sem_LIBS)
	@echo ""

${E}/xsmooth_sem_glob: $(xsmooth_sem_glob_OBJECTS) $(xsmooth_sem_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA_SEM)
	@echo ""
	${FCLINK} -o $@ $+ $(xsmooth_sem_LIBS)
	@echo ""

${E}/xsmooth_sem_sph_pde: $(xsmooth_sem_sph_pde_OBJECTS) $(xsmooth_sem_sph_pde_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xsmooth_sem_sph_pde"
	@echo ""
	${FCLINK} -o $@ $+ $(xsmooth_sem_LIBS)
	@echo ""

${E}/xsmooth_sem_pde: $(xsmooth_sem_pde_OBJECTS) $(xsmooth_sem_pde_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xsmooth_sem_pde"
	@echo ""
	${FCLINK} -o $@ $+ $(xsmooth_sem_LIBS)
	@echo ""

#######################################

###
### Module dependencies
###
$O/postprocess_par.fwat_postproc_module.o: $O/shared_par.shared_module.o

####
#### rule for each .o file below
####

$O/%.fwat_postproc_module.o: $S/%.f90 ${SETUP}/constants_tomography.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_postproc.o: $S/%.f90 ${SETUP}/constants_tomography.h $O/postprocess_par.fwat_postproc_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_postproc.o: $S/%.F90 ${SETUP}/constants_tomography.h $O/postprocess_par.fwat_postproc_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.fwat_postproc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<

###
### CUDA
###
$O/%.fwat_postproc.cuda.o: $S/%.cu ${SETUP}/config.h $S/smooth_cuda.h
	${NVCC} -c $< -o $@ $(NVCC_FLAGS)

$(cuda_smooth_sem_DEVICE_OBJ): $(cuda_smooth_sem_OBJECTS)
	${NVCCLINK} -o $(cuda_smooth_sem_DEVICE_OBJ) $(cuda_smooth_sem_OBJECTS)
