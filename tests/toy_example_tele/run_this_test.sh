#!/bin/bash

# rm -rf DATABASES_MPI
rm -rf OUTPUT_FILES

# mkdir -p DATABASES_MPI
mkdir -p OUTPUT_FILES


# copy model to DATABASES_MPI
cp model_target/* DATABASES_MPI/

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`
mpirun --oversubscribe -np $NPROC ../../bin/xmeshfem3D
mpirun --oversubscribe -np $NPROC ../../bin/xgenerate_databases
mpirun --oversubscribe -np $NPROC ../../bin/xfwat_fwd_measure_adj -m M00 -d tele -r 2 -e 1

