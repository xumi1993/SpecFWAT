#!/bin/bash
#SBATCH --nodes=3  
#SBATCH --ntasks=120
#SBATCH --time=00:30:00
#SBATCH --job-name FWD10
#SBATCH --output=FWD10_%j.txt


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



MODEL=M10
SETB=set1
SETE=set4
mpirun -np 120 ./bin/xfwat2_postproc_opt $MODEL $SETB $SETE

for step in `grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`;do
  echo "======================================="
  echo  Meshing for model $MODEL step:$step
  echo "========================================"
  #=====
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${MODEL}_step${step}" DATA/Par_file
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${MODEL}_step${step}" DATA/meshfem3D_files/Mesh_Par_file
  sed -i '/SAVE_MESH_FILES                 =/c\SAVE_MESH_FILES                 = .false.' DATA/Par_file
 
  mpirun -np 120 ./bin/xmeshfem3D
  mpirun -np 120 ./bin/xgenerate_databases
done



