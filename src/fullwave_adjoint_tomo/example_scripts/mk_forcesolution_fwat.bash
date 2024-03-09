#!/bin/bash
# format of station_used.lst
#net.stnm     lon           lat          ele
#T1.KCD01    103.91800     30.96120      589.0
#T1.KCD02    103.75800     30.90980      630.0
#T1.KCD03    103.85000     30.71890      537.0
#T1.KCD04    103.51900     30.57460      527.0
mkdir -p src_rec
cat station_used.lst |awk '{printf"%s_Z %f %f %f 0.0\n",$1,$3,$2,$4}' >src_rec/sources.dat
if false;then
ntot=`cat src_rec/sources.dat |wc -l`
ntot1=`echo $ntot |awk '{print $1-1}'`
nlen=5
npart=11
i=1
for i in `seq $npart`;do
   ib=`echo $i |awk '{print ($1-1)*a+1}' a=$nlen`
   ie=`echo $i |awk '{print $1*a}' a=$nlen`
   if [ $ie -gt $ntot1 ];then
      ie=$ntot1
   fi
   echo $i $ib $ie
   cat src_rec/sources.dat |sed -n "$ib,${ie}p" >src_rec/sources_set$i.dat
done
fi

cat src_rec/sources.dat |
while read evtfile;do
   eid=`echo $evtfile |awk '{print $1}'`
   stnm=`echo $eid |awk -F'[.,_]' '{print $2}'`
   evla=`echo $evtfile |awk '{print $2}'`
   evlo=`echo $evtfile |awk '{print $3}'`
   echo $eid: $stnm $evla $evlo 

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

#grep -v "^${eid}" src_rec/sources.dat |awk '{print substr($1,4,4),substr($1,1,2),$2,$3,$4,$5}' >src_rec/STATIONS_${eid}
cat /dev/null >src_rec/STATIONS_${eid}
for sacf in `ls data/$eid/*BXZ.sac`;do
  stnm=`echo $sacf |awk -Fdata/$eid/ '{print $2}' |awk -F. '{print $2}'` 
  grep $stnm src_rec/sources.dat |awk '{print substr($1,4,4),substr($1,1,2),$2,$3,$4,$5}' >>src_rec/STATIONS_${eid}
done
done


