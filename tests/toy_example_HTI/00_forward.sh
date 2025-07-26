#!/bin/bash
set -e

rm -rf OUTPUT_FILES

mkdir -p DATABASES_MPI
mkdir -p OUTPUT_FILES

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`

cp target_model.h5 DATA/tomo_files/tomography_model.h5
# # Run Forward to generate data
mpirun -np $NPROC ../../bin/xfwat_mesh_databases -s noise
mpirun -np $NPROC ../../bin/xfwat_fwd_measure_adj -m M00 -s noise -r 1