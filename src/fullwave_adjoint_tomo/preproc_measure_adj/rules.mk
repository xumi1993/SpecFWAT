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

## compilation directories
S := ${S_TOP}/fullwave_adjoint_tomo/prepoc_measure_adj
$(fullwave_adjoint_tomo/prepoc_measure_adj_OBJECTS): S := ${S_TOP}/src/fullwave_adjoint_tomo/prepoc_measure_adj

#######################################

fullwave_adjoint_tomo/prepoc_measure_adj_TARGETS = \
	$E/xmeasure_adj \
	$(EMPTY_MACRO)

fullwave_adjoint_tomo/prepoc_measure_adj_OBJECTS = \
	$(xmeasure_adj_OBJECTS) \
	$(EMPTY_MACRO)

fullwave_adjoint_tomo/prepoc_measure_adj_SHARED_OBJECTS = \
	$(xmeasure_adj_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

fullwave_adjoint_tomo/prepoc_measure_adj_MODULES = \
	$(FC_MODDIR)/fullwave_adjoint_tomo_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

####
#### rules for executables
####

.PHONY: measure_adj


measure_adj: $(fullwave_adjoint_tomo/prepoc_measure_adj_TARGETS)

prepoc_measure_adj: measure_adj

fullwave_adjoint_tomo/prepoc_measure_adj: measure_adj

### single targets

xmeasure_adj: $E/xmeasure_adj



#######################################

####
#### rules for each program follow
####

#######################################


##
## xmeasure_adj
##
xmeasure_adj_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/clip_sem.postprocess.o \
	$(EMPTY_MACRO)

xmeasure_adj_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

${E}/xmeasure_adj: $(xmeasure_adj_OBJECTS) $(xmeasure_adj_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xmeasure_adj"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""


#######################################

###
### Module dependencies
###

####
#### rule for each .o file below
####

$O/%.postprocess_module.o: $S/%.f90 ${SETUP}/constants_fullwave_adjoint_tomo.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.f90 ${SETUP}/constants_fullwave_adjoint_tomo.h $O/postprocess_par.postprocess_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.F90 ${SETUP}/constants_fullwave_adjoint_tomo.h $O/postprocess_par.postprocess_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<

###
### CUDA
###
$O/%.postprocess.cuda.o: $S/%.cu ${SETUP}/config.h $S/smooth_cuda.h
	${NVCC} -c $< -o $@ $(NVCC_FLAGS)

$(cuda_smooth_sem_DEVICE_OBJ): $(cuda_smooth_sem_OBJECTS)
	${NVCCLINK} -o $(cuda_smooth_sem_DEVICE_OBJ) $(cuda_smooth_sem_OBJECTS)
