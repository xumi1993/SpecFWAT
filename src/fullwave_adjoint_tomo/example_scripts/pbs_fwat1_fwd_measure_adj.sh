#!/bin/bash -l
#PBS -l nodes=5:ppn=16
#PBS -l walltime=01:30:00
#PBS -N M07.FWD18
#PBS -q hpq

cd $PBS_O_WORKDIR

module load intel/intel-18 openmpi/3.0.0-intel-18



#mpirun -np 80 ./bin/xmeshfem3D
#mpirun -np 80 ./bin/xgenerate_databases
#exit

MODEL=M07
SET=set18
set_num=`echo $SET |awk -Fset '{printf"%d\n",$2}'`
if [ $set_num -le 11 ];then
  mpirun -np 80 ./bin/xfwat1_fwd_measure_adj $MODEL $SET noise
else
  mpirun -np 80 ./bin/xfwat1_fwd_measure_adj $MODEL $SET tele
fi  
#mpirun -np 80 ./bin/xfwat0_forward_data $MODEL $SET noise

