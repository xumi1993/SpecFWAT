#!/bin/bash


ln -sf ../bin .
ln -sf ../DATA .
ln -sf ../OUTPUT_FILES .
seq 0 1 99 >slices.lst

in_path=../OUTPUT_FILES/DATABASES_MPI/
out_path=VTKmodel_M00

if [ ! -d $out_path ];then
mkdir $out_path
for tag in rho vp vs; do
#for tag in rho_smooth vp_smooth vs_smooth; do
    ../bin/xcombine_vol_data_vtk slices.lst $tag $in_path $out_path 0
done
else
  echo "$out_path already exist!!!"
fi

