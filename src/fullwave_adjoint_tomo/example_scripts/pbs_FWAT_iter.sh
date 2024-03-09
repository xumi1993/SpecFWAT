#!/bin/bash -l
#PBS -l nodes=5:ppn=16
#PBS -l walltime=09:30:00
#PBS -N FWAT.ITER
#PBS -q hpq

cd $PBS_O_WORKDIR

module load intel/intel-18 openmpi/3.0.0-intel-18




SET=set0

for MODEL in M00 M01 M02 M03;do
#for MODEL in M04 M05 M06 M07 M08 M09 M10;do
#for MODEL in M07 M08 M09 M10;do
#for MODEL in M10 M11 M12 M13 M14;do
  mpirun -np 80 ./bin/xmeshfem3D
  mpirun -np 80 ./bin/xgenerate_databases
  mpirun -np 80 ./bin/xfwat1_fwd_measure_adj $MODEL $SET noise
  mpirun -np 80 ./bin/xfwat2_postproc_opt $MODEL $SET $SET
done

