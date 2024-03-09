#!/bin/bash
#SBATCH --nodes=3  
#SBATCH --ntasks=120
#SBATCH --time=00:20:00
#SBATCH --job-name MODEL
#SBATCH --output=M18_%j.txt

# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#=====
#module load intel/15.0.2 openmpi/intel/1.6.4
#cd $PBS_O_WORKDIR
module load intel openmpi
cd $SLURM_SUBMIT_DIR
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/intel/lib/intel64:/usr/local/openmpi/lib
#=====
# number of processes
NPROC=120

prog=~/progs/SEM_tools/model_slice/sem_model_slice

mod=M05
oldmod=`echo $mod |awk -FM '{printf"M%d",$2-1}'`
#model_dir=/scratch/l/liuqy/kai/Cecal_ANAT_Rayl/optimize/SUM_KERNELS_M00/
#model_dir=../solver/${mod}/5/EKERNEL
#model_dir=../solver/${mod}/GRADIENT/
#model_dir=../optimize/SUM_KERNELS_${mod}/
model_dir=../optimize/MODEL_${mod}/


topo_dir=../OUTPUT_FILES/DATABASES_MPI
#topo_dir=OUTPUT_FILES/DATABASES_MPI

# xyz_infile topo_dir model_dir data_name gmt_outfile
# Note: set the right NSPEC_AB when compile the code

#for kern in vpv vph rho;do
#for kern in alpha_kernel beta_kernel rhop_kernel;do #r_aniso voigt_vs;do
#for kern in alpha_kernel beta_kernel hess_kernel;do #r_aniso voigt_vs;do
#for kern in alpha_kernel_smooth beta_kernel_smooth rho_kernel_smooth;do #r_aniso voigt_vs;do
for kern in vp vs rho;do
#mpirun -np $NPROC $prog model_regrid.xyz $topo_dir $model_dir ${kern}  grad.${kern}.${mod}.5.regridded.xyz 
#mpirun -np $NPROC $prog model_regrid.xyz $topo_dir $model_dir ${kern}  grad.${kern}.${mod}.regridded.xyz 
mpirun -np $NPROC $prog model_regrid.xyz $topo_dir $model_dir ${kern}   model.${kern}.${mod}.regridded.xyz 
done
