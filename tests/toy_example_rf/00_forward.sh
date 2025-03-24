#!/bin/bash
set -e

rm -rf OUTPUT_FILES

mkdir -p DATABASES_MPI
mkdir -p OUTPUT_FILES

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`

mkdir -p model_target

# # Run Forward to generate data
cp model_target/* DATABASES_MPI/
mpirun --oversubscribe -np $NPROC ../../bin/xmeshfem3D
mpirun --oversubscribe -np $NPROC ../../bin/xgenerate_databases
mpirun --oversubscribe -np $NPROC ../../bin/xfwat_fwd_measure_adj -m M00 -d rf -r 1