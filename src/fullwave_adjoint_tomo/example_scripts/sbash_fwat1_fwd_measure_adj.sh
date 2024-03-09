#!/bin/bash
#SBATCH --nodes=3  
#SBATCH --ntasks=120
#SBATCH --time=00:30:00
#SBATCH --job-name FWD1
#SBATCH --output=FWD1_%j.txt


# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load intel/15.0.2 openmpi/intel/1.6.4
module load intel openmpi
#=====
#cd $PBS_O_WORKDIR
cd $SLURM_SUBMIT_DIR
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/intel/lib/intel64:/usr/local/openmpi/lib
#=====



#mpirun -np 120 ./bin/xmeshfem3D
#mpirun -np 120 ./bin/xgenerate_databases

MODEL=M00
SET=set1
TYPE=tele
mpirun -np 120 ./bin/xfwat1_fwd_measure_adj $MODEL $SET $TYPE

