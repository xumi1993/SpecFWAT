#!/bin/bash
set -e


rm -rf DATABASES_MPI
rm -rf OUTPUT_FILES

mkdir -p DATABASES_MPI
mkdir -p OUTPUT_FILES

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`

for iter in `seq 0 0`; do
    MODEL=`printf "M%02d" $iter`
    if [ $iter -eq 0 ]; then
        cp initial_model.h5 DATA/tomo_files/tomography_model.h5
    fi
    mpirun -np $NPROC ../../bin/xfwat_mesh_databases -s tele
    mpirun -np $NPROC ../../bin/xfwat_fwd_measure_adj -m $MODEL -s tele -r 3
    mpirun -np $NPROC ../../bin/xfwat_post_proc -m $MODEL
    mpirun -np $NPROC ../../bin/xfwat_optimize -m $MODEL
done