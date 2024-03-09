#!/bin/bash -l
#PBS -l nodes=5:ppn=16
#PBS -l walltime=00:16:00
#PBS -N sem_mesh
#PBS -q hpq

cd $PBS_O_WORKDIR

ulimit -S -s unlimited
module load intel/intel-18 openmpi/3.0.0-intel-18

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
sumkern_dir=model_1d/1d_ak135/
mkdir -p $sumkern_dir
for knm in rho vp vs ;do  ##### change rho
   mpirun -np 80  ./bin/xsmooth_sem  5000  10000 $knm $sumkern_dir $sumkern_dir FALSE # not GPU mode
done

mkdir model_1d/1d_ak135_smooth_h5km_v10km
cd model_1d/1d_ak135_smooth_h5km_v10km
mv ../1d_ak135/*_smooth.bin .
rename rho_smooth.bin rho.bin *rho_smooth.bin
rename vp_smooth.bin vp.bin *vp_smooth.bin
rename vs_smooth.bin vs.bin *vs_smooth.bin
cd ../..
