#!/bin/bash
mod=$1
setb=$2
sete=$3
env=pbs # slurm or pbs

###
if [ $# -ne 3 ];then
    echo " Usage: ./submit_job_fwat2.bash M?? setb sete "
    echo "narg: " $#
    exit
fi
###

if [ $env == 'slurm' ];then 
  script=sbash_fwat2_postproc_opt.sh
elif [ $env == 'pbs' ];then
  script=pbs_fwat2_postproc_opt.sh
fi

if [ $mod == "M00" ];then
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI" DATA/Par_file
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI" DATA/meshfem3D_files/Mesh_Par_file
else
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${mod}" DATA/Par_file
  sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${mod}" DATA/meshfem3D_files/Mesh_Par_file
fi 


sed -i "/MODEL=/c\MODEL=${mod}" $script
sed -i "/SETB=/c\SETB=set${setb}" $script
sed -i "/SETE=/c\SETE=set${sete}" $script

if [ $env == 'slurm' ];then 
  sbatch $script
elif [ $env == 'pbs' ];then
  qsub $script
fi
##################
