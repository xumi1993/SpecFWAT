#!/bin/bash
cube2sph=true
cube2sph_dir=/scratch/l/liuqy/kai/cube2sph/cube2sph_utils/bin
# format of stations_all.lst
#net.stnm     lon           lat          ele
#T1.KCD01    103.91800     30.96120      589.0
#T1.KCD02    103.75800     30.90980      630.0
#T1.KCD03    103.85000     30.71890      537.0
#T1.KCD04    103.51900     30.57460      527.0

if false;then
mkdir -p src_rec
cat stations_all.lst |awk '{printf"%s_Z %f %f %f 0.0\n",$1,$3,$2,$4}' >src_rec/sources.dat
fi

#~/progs/specfem3d_for_fwat_devel/bin/xdynmbat_sele_minibatch src_rec/sources.dat 20 src_rec/sources_batch.M00.dat

cat src_rec/sources.dat |sed -n '593,$p' |
while read evtfile;do
   eid=`echo $evtfile |awk '{print $1}'`
   stnm=`echo $eid |awk -F'[.,_]' '{print $2}'`
   net=`echo $eid |awk -F'[.,_]' '{print $1}'`
   evla=`echo $evtfile |awk '{print $2}'`
   evlo=`echo $evtfile |awk '{print $3}'`
   echo $eid: $stnm $evla $evlo 

   if $cub2sph;then
      module load intel openmpi

      if [ ! -f src_rec/FORCESOLUTION_${eid} ];then
      cat >DATA/FORCESOLUTION_geo <<eof
FORCE  001
time shift:     0.0000
half duration:  0.2 !Half duration (s) for Step function, frequency (Hz) for Ricker
latorUTM:       $evla
longorUTM:      $evlo
depth:          0.0
source time function: 1 ! 0=Step function, 1=Ricker wavelet
factor force source:             1.d15
component dir vect source E:     0.d0
component dir vect source N:     0.d0
component dir vect source Z_UP:  1.d0
eof
      $cube2sph_dir/write_force_solution_file
      mv DATA/FORCESOLUTION_cartesian src_rec/FORCESOLUTION_${eid} 
     else
      echo "    src_rec/FORCESOLUTION_${eid} already exists!"   
     fi

     if [ ! -f src_rec/STATIONS_${eid} ];then
      cat /dev/null >DATA/STATIONS_geo
      for sacf in `ls data/$eid/*BXZ.sac`;do
         stnm=`echo $sacf |awk -Fdata/$eid/ '{print $2}' |awk -F. '{print $2}'` 
         grep $net.$stnm stations_all.lst |awk '{print substr($1,4),substr($1,1,2),$3,$2,$4,0.0}' >>DATA/STATIONS
      done
      #grep -v "^${net}.${stnm}" stations_all.lst |awk '{printf"%s %s %.5f %.5f %.1f 0.0\n",substr($1,4),substr($1,1,2),$3,$2,$4}' >DATA/STATIONS

      $cube2sph_dir/write_stations_file
      mv DATA/STATIONS_cartesian  src_rec/STATIONS_${eid}
     else
      echo "    src_rec/STATIONS_${eid} already exists!"
     fi

   else
      cat >src_rec/FORCESOLUTION_${eid} <<eof
FORCE  001
time shift:     0.0000
f0:             1.0
latorUTM:       $evla
longorUTM:      $evlo
depth:          0.0000
factor force source:             1.d15
component dir vect source E:     0.d0
component dir vect source N:     0.d0
component dir vect source Z_UP:  1.d0
eof

      cat /dev/null >src_rec/STATIONS_${eid}
      for sacf in `ls data/$eid/*BXZ.sac`;do
         stnm=`echo $sacf |awk -Fdata/$eid/ '{print $2}' |awk -F. '{print $2}'` 
         grep $stnm stations_all.lst |awk '{print substr($1,4),substr($1,1,2),$3,$2,$4,$5}' >>src_rec/STATIONS_${eid}
      done
   fi
done


