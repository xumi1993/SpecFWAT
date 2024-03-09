#!/bin/bash -l
#PBS -l nodes=5:ppn=16
#PBS -l walltime=01:30:00
#PBS -N M06.LS0.050
#PBS -q hpq

cd $PBS_O_WORKDIR

module load intel/intel-18 openmpi/3.0.0-intel-18




MODEL=M06_step0.050
SET=ls
SIMU_TYPE=noise

mpirun -np 80 ./bin/xfwat3_linesearch $MODEL $SET $SIMU_TYPE

