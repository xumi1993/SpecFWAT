#!/bin/bash -l

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

mkdir -p model_initial
setpar_fwat DATA/Par_file MODEL gll

# Decompose the initial model
for name in vs vp rho; do
  mpirun -np $NPROC ../../bin/xdecompose_h5_gll ${name} ./initial_model.h5 model_initial/
done

cp model_initial/* DATABASES_MPI/

# Run the fullwaveform adjoint tomography with 2 iterations
mpirun -np $NPROC ../../bin/xfullwave_adjoint_tomo 0 2

