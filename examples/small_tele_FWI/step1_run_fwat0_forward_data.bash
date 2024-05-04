#!/bin/bash -l
#PBS -l nodes=5:ppn=16
#PBS -l walltime=00:05:00
#PBS -N FWAT.ITER
#PBS -q hpq

#cd $PBS_O_WORKDIR

#module load intel/intel-18 openmpi/3.0.0-intel-18

NPROC=4
cp model_target/* ./DATABASES_MPI
sed -i "/MODEL                           =/c\MODEL                           = gll" DATA/Par_file

for SET in set0;do
  echo "Running SET: " $SET
  mpirun -np $NPROC ../../bin/xmeshfem3D
  mpirun -np $NPROC ../../bin/xgenerate_databases
  # mpirun -np $NPROC ../../bin/xfwat0_forward_data M00 $SET rf
done
