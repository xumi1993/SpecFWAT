#!/bin/bash -l

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

SET=set1

mpirun -np $NPROC ../../bin/xfwat_mesh_database rf
mpirun -np $NPROC ../../bin/xfwat0_forward_data M00 $SET rf
