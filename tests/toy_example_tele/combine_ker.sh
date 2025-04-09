#!/bin/bash
cat > grid.txt << eof
833950 -44274.0 -80000
4000 4000 2000
51 23 41
eof
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

echo ${NPROC}
mpirun -np $NPROC ../../bin/xmodel_grid_cart_custom grid.txt beta_kernel_smooth optimize/SUM_KERNELS_M00/ ./ M00 false

rm grid.txt