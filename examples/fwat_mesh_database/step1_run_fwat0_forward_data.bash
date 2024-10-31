#!/bin/bash -l

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

setpar_fwat DATA/FWAT.PAR INITIAL_MODEL_PATH ./DATA/target_model.h5
setpar_fwat DATA/FWAT.PAR USE_H5 .true.

# mpirun -np $NPROC ../../bin/xfwat_mesh_database rf
# mpirun -np $NPROC ../../bin/xfwat0_forward_data M00 set1 rf

mpirun -np $NPROC ../../bin/xfwat_mesh_database noise
# mpirun -np $NPROC ../../bin/xfwat0_forward_data M00 set0 noise
