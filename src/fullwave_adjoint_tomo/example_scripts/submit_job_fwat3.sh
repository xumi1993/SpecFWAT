#!/bin/bash
mod=$1
simu_type=$2
env=pbs # slurm or pbs

###
if [ $# -ne 2 ];then
    echo " Usage: ./submit_fwat3.bash M?? noise/tele "
    echo "narg: " $#
    exit
fi
###


if [ $env == 'slurm' ];then
   fwd=sbash_fwat3_linesearch.sh 
elif [ $env == 'pbs' ];then
   fwd=pbs_fwat3_linesearch.sh 
fi



for ipart in `grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`;do
#for ipart in 0.020;do
  echo "======================================="
  echo Model:$mod Part:$ipart
  echo "========================================"
  if [ $env == 'slurm' ];then
  #######   for slurm ######
    sed -i "/#SBATCH --job-name/c\#SBATCH --job-name LS${ipart}" $fwd 
    sed -i "/#SBATCH --time/c\#SBATCH --time=00:30:00" $fwd
    sed -i "/#SBATCH --output/c\#SBATCH --output=LS${ipart}_\%j.txt" $fwd
    sed -i "/MODEL=/c\MODEL=${mod}_step${ipart}" $fwd
    sed -i "/SET=/c\SET=ls" $fwd
    sbatch $fwd
  elif [ $env == 'pbs' ];then
  #######   for PBS ######
    sed -i "/#PBS -N/c\#PBS -N ${mod}.LS${ipart}" $fwd 
    sed -i "/#PBS -l walltime/c\#PBS -l walltime=01:30:00" $fwd
    sed -i "/MODEL=/c\MODEL=${mod}_step${ipart}" $fwd
    sed -i "/SET=/c\SET=ls" $fwd
    sed -i "/SIMU_TYPE=/c\SIMU_TYPE=${simu_type}" $fwd
    qsub $fwd
    sleep 3
  fi
done
##################
