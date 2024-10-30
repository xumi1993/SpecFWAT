#!/bin/bash -l

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`


setpar_fwat DATA/FWAT.PAR INITIAL_MODEL_PATH ./DATA/initial_model.h5
setpar_fwat DATA/FWAT.PAR USE_H5 .true.

# Run the fullwaveform adjoint tomography with 2 iterations
mpirun -np $NPROC ../../bin/xfullwave_adjoint_tomo 0 1

