#!/bin/bash

# rm -rf DATABASES_MPI
rm -rf OUTPUT_FILES

mkdir -p DATABASES_MPI
mkdir -p OUTPUT_FILES

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`

# # Run Forward to generate data
# cp model_target/* DATABASES_MPI/
# mpirun --oversubscribe -np $NPROC ../../bin/xmeshfem3D
# mpirun --oversubscribe -np $NPROC ../../bin/xgenerate_databases
# mpirun --oversubscribe -np $NPROC ../../bin/xfwat_fwd_measure_adj -m M00 -d tele -r 1

# Run Inversion
mkdir -p model_initial
# for key in vp vs rho; do
#     mpirun --oversubscribe -np $NPROC ../../bin/xdecompose_h5_gll $key ./initial_model.h5 ./model_initial/ 
# done
cp model_initial/* DATABASES_MPI/
mpirun --oversubscribe -np $NPROC ../../bin/xmeshfem3D
mpirun --oversubscribe -np $NPROC ../../bin/xgenerate_databases
mpirun --oversubscribe -np $NPROC ../../bin/xfwat_fwd_measure_adj -m M00 -d tele -r 3
