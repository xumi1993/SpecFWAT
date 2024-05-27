#!/bin/bash -l

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

mkdir -p model_initial
setpar_fwat DATA/Par_file MODEL gll

for name in vs vp rho; do
  mpirun -np $NPROC ../../bin/xdecompose_h5_gll ${name} ./initial_model.h5 model_initial/
done

mpirun -np $NPROC ../../bin/xfullwave_adjoint_tomo 0 2

# SET=set0
# for i in  `seq 0 1 2`;do
#   MODEL=`echo $i|awk '{printf "M%02d","'$i'"}'`
#   if [ $MODEL == 'M00' ]; then
#     cp model_initial/* ./DATABASES_MPI
#   fi
#   mpirun -np $NPROC ../../bin/xmeshfem3D
#   mpirun -np $NPROC ../../bin/xgenerate_databases
#   mpirun -np $NPROC ../../bin/xfwat1_fwd_measure_adj $MODEL $SET rf 3
#   mpirun -np $NPROC ../../bin/xfwat2_postproc_opt $MODEL
# done

