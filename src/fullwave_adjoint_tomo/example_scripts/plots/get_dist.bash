#!/bin/bash

cat >Geo_loc.txt <<eof
-120.67 36.07 SAF
-120.33 36.13  GV-L
-119.14 36.34  GV-R
-118.09 36.50 SNB-R
-117.5 36.63 Ap
eof

evlo=-121.56
evla=35.93
cat Geo_loc.txt |while read loc;do
   stlo=`echo $loc |awk '{print $1}'`
   stla=`echo $loc |awk '{print $2}'`
   stnm=`echo $loc |awk '{print $3}'`
   echo $stnm $evlo $evla $stlo $stla
   ~/progs/sactools_c/distbaz $evlo $evla $stlo $stla
done
