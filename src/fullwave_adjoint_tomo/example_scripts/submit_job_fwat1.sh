#!/bin/bash
mod=$1
simu_type=$2
setb=$3
sete=$4
env=pbs # slurm or pbs

###
if [ $# -ne 4 ];then
    echo " Usage: ./submit_job_fwat1.bash M?? noise/tele setb sete "
    echo "narg: " $#
    exit
fi
###

if [ $simu_type == 'noise' ];then
  echo "copying params for noise"
  cp DATA/Par_file.noise DATA/Par_file
  cp fwat_params/FWAT.PAR.noise fwat_params/FWAT.PAR
  cp fwat_params/MEASUREMENT.PAR.noise fwat_params/MEASUREMENT.PAR
  cp src_rec/sources_ls.dat.noise src_rec/sources_ls.dat
elif [ $simu_type == 'tele' ];then
  echo "copying params for tele"
  cp DATA/Par_file.tele DATA/Par_file
  cp fwat_params/FWAT.PAR.tele fwat_params/FWAT.PAR
  cp fwat_params/MEASUREMENT.PAR.tele fwat_params/MEASUREMENT.PAR
  cp src_rec/sources_ls.dat.tele src_rec/sources_ls.dat
fi 

if [ $env == 'slurm' ];then
   fwd=sbash_fwat1_fwd_measure_adj.sh
elif [ $env == 'pbs' ];then
   fwd=pbs_fwat1_fwd_measure_adj.sh
fi

for ipart in `seq $setb 1 $sete`;do
  if [ $mod == "M00" ];then
    echo "======================================="
    echo Model:$mod Part:$ipart
    sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI" DATA/Par_file
    sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./OUTPUT_FILES/DATABASES_MPI" DATA/meshfem3D_files/Mesh_Par_file
  else
    echo "======================================="
    echo Model:$mod Part:$ipart
    sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${mod}" DATA/Par_file
    sed -i "/LOCAL_PATH                      =/c\LOCAL_PATH                      = ./optimize/MODEL_${mod}" DATA/meshfem3D_files/Mesh_Par_file
  fi 

  echo "========================================"

  if [ $env == 'slurm' ];then
  #######   for slurm ######
    sed -i "/#SBATCH --job-name/c\#SBATCH --job-name FWD${ipart}" $fwd 
    sed -i "/#SBATCH --time/c\#SBATCH --time=00:30:00" $fwd
    sed -i "/#SBATCH --output/c\#SBATCH --output=FWD${ipart}_\%j.txt" $fwd
    sed -i "/MODEL=/c\MODEL=${mod}" $fwd
    sed -i "/SET=/c\SET=set${ipart}" $fwd
    sbatch $fwd
  elif [ $env == 'pbs' ];then
  #######   for PBS ######
    sed -i "/#PBS -N/c\#PBS -N ${mod}.FWD${ipart}" $fwd 
    sed -i "/#PBS -l walltime/c\#PBS -l walltime=01:30:00" $fwd
    sed -i "/MODEL=/c\MODEL=${mod}" $fwd
    sed -i "/SET=/c\SET=set${ipart}" $fwd
    qsub $fwd
    sleep 3
  fi
done
##################
