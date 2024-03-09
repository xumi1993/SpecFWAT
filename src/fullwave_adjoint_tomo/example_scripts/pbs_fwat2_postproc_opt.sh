#!/bin/bash -l
#PBS -l nodes=5:ppn=16
#PBS -l walltime=00:30:00
#PBS -N M00.OPT
#PBS -q sandyq

cd $PBS_O_WORKDIR

module load intel/intel-18 openmpi/3.0.0-intel-18




MODEL=M06
SETB=set1
SETE=set11
mpirun -np 80 ./bin/xfwat2_postproc_opt $MODEL $SETB $SETE

for step in `grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`;do
  echo "======================================="
  echo  Meshing for model $MODEL step:$step
  echo "========================================"
  #=====
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${MODEL}_step${step}" DATA/Par_file
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${MODEL}_step${step}" DATA/meshfem3D_files/Mesh_Par_file
  sed -i '/SAVE_MESH_FILES                 =/c\SAVE_MESH_FILES                 = .false.' DATA/Par_file
 
  mpirun -np 80 ./bin/xmeshfem3D
  mpirun -np 80 ./bin/xgenerate_databases
done



