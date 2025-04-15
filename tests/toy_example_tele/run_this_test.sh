#!/bin/bash
set -e

# rm -rf DATABASES_MPI
rm -rf OUTPUT_FILES

mkdir -p DATABASES_MPI
mkdir -p OUTPUT_FILES

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`


# Run Inversion
for it in `seq 0 0`; do
    model=`printf "M%02d" $it`
    if [ $it -eq 0 ]; then
        cp initial_model.h5 DATA/tomo_files/tomography_model.h5
    else
        cp optimize/model_${model}.h5 DATA/tomo_files/tomography_model.h5
    fi
    mpirun -np $NPROC ../../bin/xfwat_mesh_databases -s tele
    mpirun -np $NPROC ../../bin/xfwat_fwd_measure_adj -m $model -s tele -r 3
    mpirun -np $NPROC ../../bin/xfwat_post_proc -m $model
    mpirun -np $NPROC ../../bin/xfwat_optimize -m $model
done
