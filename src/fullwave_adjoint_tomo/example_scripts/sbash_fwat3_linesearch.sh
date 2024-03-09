#!/bin/bash
#SBATCH --nodes=3  
#SBATCH --ntasks=120
#SBATCH --time=00:30:00
#SBATCH --job-name LS0.050
#SBATCH --output=LS0.050_%j.txt


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

MODEL=M10_step0.050
SET=ls
mpirun -np 120 ./bin/xfwat3_linesearch $MODEL $SET tele

