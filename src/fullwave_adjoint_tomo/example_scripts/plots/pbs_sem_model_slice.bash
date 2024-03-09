#!/bin/bash -l
#PBS -l nodes=5:ppn=16
#PBS -l walltime=01:00:00
#PBS -N sem_slice
#PBS -q hpq

cd $PBS_O_WORKDIR

ulimit -S -s unlimited
module load intel/intel-18 openmpi/3.0.0-intel-18

#=====
# number of processes
NPROC=80

prog=~/progs/FWAT_tools/model_slice/sem_model_slice_opt_Cecal

mod=M03
oldmod=`echo $mod |awk -FM '{printf"M%d",$2-1}'`
#model_dir=/scratch/l/liuqy/kai/Cecal_ANAT_Rayl/optimize/SUM_KERNELS_M00/
#model_dir=../solver/${mod}/5/EKERNEL
#model_dir=../solver/${mod}/GRADIENT/
#model_dir=../optimize/SUM_KERNELS_${mod}/
#model_dir=../model_1d/1d_ak135sm_chb3Dx1_rotate
model_dir=../optimize/MODEL_${mod}

topo_dir=../OUTPUT_FILES/DATABASES_MPI
#topo_dir=OUTPUT_FILES/DATABASES_MPI

# xyz_infile topo_dir model_dir data_name gmt_outfile
# Note: set the right NSPEC_AB when compile the code

#for kern in vpv vph rho;do
#for kern in alpha_kernel beta_kernel rhop_kernel;do #r_aniso voigt_vs;do
#for kern in alpha_kernel beta_kernel hess_kernel;do #r_aniso voigt_vs;do
#for kern in alpha_kernel_smooth beta_kernel_smooth rho_kernel_smooth;do #r_aniso voigt_vs;do
for kern in rho vp vs;do
#mpirun -np $NPROC $prog model_regrid.xyz $topo_dir $model_dir ${kern}  grad.${kern}.${mod}.5.regridded.xyz 
#mpirun -np $NPROC $prog model_regrid.xyz $topo_dir $model_dir ${kern}  grad.${kern}.${mod}.regridded.xyz 
mpirun -np $NPROC $prog model_regrid.xyz $topo_dir $model_dir ${kern}   model.${kern}.${mod}.regridded.xyz 
done
