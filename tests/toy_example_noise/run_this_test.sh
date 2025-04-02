#!/bin/bash
set -e

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`

# Run Inversion
mkdir -p model_initial
for key in vp vs rho; do
    mpirun -np $NPROC ../../bin/xdecompose_h5_gll $key ./initial_model.h5 ./model_initial/ 
done
cp model_initial/* DATABASES_MPI/

for iter in `seq 0 0`; do
    MODEL=`printf "M%02d" $iter`
    mpirun -np $NPROC ../../bin/xmeshfem3D
    mpirun -np $NPROC ../../bin/xgenerate_databases
    mpirun -np $NPROC ../../bin/xfwat_fwd_measure_adj -m $MODEL -s noise -r 3
    mpirun -np $NPROC ../../bin/xfwat_post_proc -m $MODEL
    # mpirun -np $NPROC ../../bin/xfwat_optimize -m $MODEL
done