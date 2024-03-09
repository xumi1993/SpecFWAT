#!/bin/bash -l
#PBS -l nodes=5:ppn=16
#PBS -l walltime=00:16:00
#PBS -N addpert
#PBS -q hpq

cd $PBS_O_WORKDIR

ulimit -S -s unlimited
module load intel/intel-18 openmpi/3.0.0-intel-18

prog=~/progs/FWAT_tools/model_addpert

model_dir=model_1d/1d_ak135_smooth_h5km_v10km

topo_dir=./OUTPUT_FILES/DATABASES_MPI
# xyz_infile topo_dir model_dir data_name gmt_outfile
# Note: set the right NSPEC_AB when compile the code
out_dir=model_1d/1d_ak135sm_chb3Dx2_rotate
mkdir -p $out_dir

for tag in rho vp vs;do
   mpirun -np 80 $prog/sem_model_checkboard checkboard_grid.dat $topo_dir $model_dir $out_dir $tag 10000
   #mpirun $prog/sem_model_checkboard gaus01.dat $topo_dir $model_dir $out_dir $tag 10000
done
