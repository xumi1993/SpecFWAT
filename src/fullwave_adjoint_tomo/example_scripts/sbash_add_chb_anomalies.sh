#!/bin/bash -l
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --time=10:00:00
#SBATCH --job-name fwi_inv
#SBATCH --output=fwi_inv_%j.txt

ulimit -S -s unlimited

module load intel openmpi

cd $SLURM_SUBMIT_DIR

prog=/home/l/liuqy/kai/progs/SEM_tools/model_addpert

model_dir=model_1d/1d_ak135_smooth_h5km_v10km

topo_dir=./OUTPUT_FILES/DATABASES_MPI
# xyz_infile topo_dir model_dir data_name gmt_outfile
# Note: set the right NSPEC_AB when compile the code
out_dir=model_1d/1d_ak135sm_chb_rotate
mkdir -p $out_dir

for tag in rho vp vs;do
   mpirun $prog/sem_model_checkboard checkboard_grid.dat $topo_dir $model_dir $out_dir $tag 10000
   #mpirun $prog/sem_model_checkboard gaus01.dat $topo_dir $model_dir $out_dir $tag 10000
done
