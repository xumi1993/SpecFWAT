#!/bin/bash
set -e

rm -rf OUTPUT_FILES

mkdir -p DATABASES_MPI
mkdir -p OUTPUT_FILES

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`

setpar_fwat DATA/Par_file MODEL default
mpirun -np $NPROC ../../bin/xmeshfem3D
mpirun -np $NPROC ../../bin/xgenerate_databases

setpar_fwat DATA/Par_file MODEL gll
mkdir -p model_target
for key in vp vs rho; do
    mpirun -np $NPROC ../../bin/xdecompose_h5_gll $key ./target_model.h5 ./model_target/ 
done
cp model_target/* DATABASES_MPI/

# # Run Forward to generate data
mpirun -np $NPROC ../../bin/xmeshfem3D
mpirun -np $NPROC ../../bin/xgenerate_databases
mpirun -np $NPROC ../../bin/xfwat_fwd_measure_adj -m M00 -s noise -r 1