#!/bin/bash

# For mini-batch
#~/progs/minibatch_opt/sampling_minibatch src_rec/sources.dat 6 src_rec/sources_M00.dat

rm -f src_rec/sources_set*.dat
srclist=src_rec/sources_M00.dat

ntot=`cat $srclist |wc -l`
ntot1=`echo $ntot |awk '{print $1-1}'`
nlen=1
npart=6
i=1
for i in `seq $npart`;do
   ib=`echo $i |awk '{print ($1-1)*a+1}' a=$nlen`
   ie=`echo $i |awk '{print $1*a}' a=$nlen`
   if [ $ie -gt $ntot ];then
      ie=$ntot
   fi
   echo $i $ib $ie
   cat $srclist |sed -n "$ib,${ie}p" >src_rec/sources_set$i.dat
done

