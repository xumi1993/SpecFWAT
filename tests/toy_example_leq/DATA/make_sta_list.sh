#!/bin/bash

cat /dev/null > STATIONS
i=1
for y in `seq -30 4 30`;do
   for x in `seq -30 4 30`;do
   stnm=`echo $i |awk '{printf"S%03d",$1}'`
   echo $stnm WK $y $x |awk '{printf"%s %s %10.2f %10.2f 0.0 0.0\n",$1,$2,$3*1000,$4*1000}' >> STATIONS
   let i=i+1
   done
done

