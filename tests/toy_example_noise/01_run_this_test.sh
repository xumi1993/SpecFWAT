#!/bin/bash
set -e

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`

# Run Inversion

for it in `seq 0 14`; do
    MODEL=`printf "M%02d" $it`
    if [ $it -eq 0 ]; then
        cp initial_model.h5 DATA/tomo_files/tomography_model.h5
    fi
    mpirun -np $NPROC ../../bin/xfwat_mesh_databases -s noise
    mpirun -np $NPROC ../../bin/xfwat_fwd_measure_adj -m $MODEL -s noise -r 3
    mpirun -np $NPROC ../../bin/xfwat_post_proc -m $MODEL
    mpirun -np $NPROC ../../bin/xfwat_optimize -m $MODEL
done