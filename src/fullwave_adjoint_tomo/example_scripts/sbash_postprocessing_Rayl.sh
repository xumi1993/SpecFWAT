#!/bin/bash
#SBATCH --nodes=3  
#SBATCH --ntasks=120
#SBATCH --time=00:16:00
#SBATCH --job-name POST
#SBATCH --output=POST_%j.txt





# script runs mesher,database generation and solver
# using this example setup
#
###################################################
#module load intel/15.0.2 openmpi/intel/1.6.4
module load intel/2018.2 #intelmpi/2018.2
module load openmpi/3.1.0rc3
#=====
#cd $PBS_O_WORKDIR
cd $SLURM_SUBMIT_DIR
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/intel/lib/intel64:/usr/local/openmpi/lib
#=====

cdir=`pwd`
MODEL=M00
is_ls=false
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
#######################################

oldmod=`echo $MODEL |awk -FM '{printf"M%02d",$2-1}'`
newmod=`echo $MODEL |awk -FM '{printf"M%02d",$2+1}'`
sumkern_dir=SUM_KERNELS_$MODEL/
is_sumkern=true
is_sd_update=true
is_cg_update=false
is_LBFGS_update=false




#######################################################
if $is_sumkern;then
echo "Post processing"
cd optimize_tomo
#mpirun -np $NPROC ./bin/xsum_kernels kernels_list_$MODEL.txt $sumkern_dir
mpirun -np $NPROC ./bin/xsum_preconditioned_kernels kernels_list_$MODEL.txt $sumkern_dir
for knm in beta_kernel alpha_kernel rho_kernel;do  ##### change to rhop_kernel for any density inversion
   mpirun -np $NPROC ./bin/xsmooth_sem  10000  10000 $knm $sumkern_dir $sumkern_dir FALSE # not GPU mode
done
echo "Post processing done!"
cd ..
fi
#######################################################
if $is_sd_update;then
date
echo "Steep Descent (SD) update ..."
cd optimize_tomo
if [ ! -d MODEL_M00 ];then
   ln -s ../initial_model MODEL_M00
fi

for step in `cat step_len.dat`;do
if $is_ls;then
   mkdir -p MODEL_${newmod}_slen${step}
   mpirun -np $NPROC ./bin/xadd_model_iso $step MODEL_$MODEL MODEL_${newmod}_slen${step} SUM_KERNELS_$oldmod SUM_KERNELS_$MODEL
else
   #mpirun -np $NPROC ./bin/xadd_model_iso $step MODEL_$MODEL MODEL_${newmod} SUM_KERNELS_$oldmod SUM_KERNELS_$MODEL
   mpirun -np $NPROC ./bin/xadd_model_iso $step MODEL_$MODEL MODEL_${newmod} none SUM_KERNELS_$MODEL
fi
done
cd ../
date
fi
#######################################################
if $is_cg_update;then
date
echo "Conjugate Gradient (CG) update ..."
cd optimize_tomo
if [ ! -d MODEL_M00 ];then
   ln -sf ../initial_model MODEL_M00
fi

for step in `cat step_len.dat`;do
if $is_ls;then
   mkdir -p MODEL_${newmod}_slen${step}
   mpirun -np $NPROC ./bin/xadd_model_iso_cg $step MODEL_$MODEL MODEL_${newmod}_slen${step} SUM_KERNELS_$oldmod SUM_KERNELS_$MODEL
else
   mpirun -np $NPROC ./bin/xadd_model_iso_cg $step MODEL_$MODEL MODEL_${newmod} SUM_KERNELS_$oldmod SUM_KERNELS_$MODEL
fi
done
cd ..
date
fi
#######################################################
if $is_LBFGS_update;then
date
echo "Limited memory BFGS  update ..."
cd optimize_tomo
iter=`echo $mod |awk -FM '{print $2}'`
for step in `cat step_len.dat`;do

if [ ! -d MODEL_M00 ];then
   ln -s ../initial_model MODEL_M00
fi
### using previous gradients and models together with current gradient and model to obtain LBFGS direction
mkdir -p MODEL_${newmod}_slen${step}
mpirun -np $NPROC ./bin/xadd_model_iso_lbfgs $step 0 $iter MODEL_$MODEL MODEL_${newmod} SUM_KERNELS_$oldmod SUM_KERNELS_$MODEL 

done
cd ..
date
fi
###
