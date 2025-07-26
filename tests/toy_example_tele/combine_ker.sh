#!/bin/bash

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

mpirun -np $NPROC ../../bin/xmodel_grid_cart_custom beta_kernel optimize/SUM_KERNELS_M00/ ./ M00 false
