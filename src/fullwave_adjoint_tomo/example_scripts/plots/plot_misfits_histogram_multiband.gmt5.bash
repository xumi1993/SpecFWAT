#!/bin/bash
module load gcc hdf5/1.8.20 netcdf/4.6.1
#        ! KEY: write misfit function values to file (two for each window)
#        ! Here are the 20 columns of the vector window_chi
#        !  1: MT-TT chi,    2: MT-dlnA chi,    3: XC-TT chi,    4: XC-dlnA chi
#        !  5: MT-TT meas,   6: MT-dlnA meas,   7: XC-TT meas,   8: XC-dlnA meas
#        !  9: MT-TT error, 10: MT-dlnA error, 11: XC-TT error, 12: XC-dlnA error
#        ! WINDOW     : 13: data power, 14: syn power, 15: (data-syn) power, 16: window duration
#        ! FULL RECORD: 17: data power, 18: syn power, 19: (data-syn) power, 20: record duration
#        ! Example of a reduced file: awk '{print $2,$3,$4,$5,$6,$31,$32}' window_chi > window_chi_sub
gmtset ANNOT_FONT_SIZE_PRIMARY 12p
 
mod=$1
cmp=$2
scratch=/scratch/l/liuqy/kai/Cecal_noise_fwat

output=dT_${mod}_${cmp}_multiband.ps

i=1
for band in T010_T020 T015_T030 T020_T040;do
input=dT_${mod}_${band}.dat
#cat ../../output/misfits/${mod}_${band}_*_chi |awk '{print $15,$19 }' >misfits_${mod}_${band}.dat
#cat window_chi |awk '{print a,$1,$15,$19,$13,$17}' a=$eid >> ../output/misfits/misfits_${mod}_${T}_BXZ.dat
cat /dev/null >$input
cat /dev/null >$input.all
if [ $mod == "M15" ];then
  ls $scratch/misfits/${mod}_${band}_*_chi >chi.dat
else
  ls $scratch/misfits/${mod}.set*_${band}_*_chi >chi.dat
fi
for file in `cat chi.dat`;do
   #echo $file
   evt=`echo $file |awk -F$scratch/misfits/ '{print $2}' |awk -F_ '{print $4 }'`

   cat $file |grep $cmp |awk '{print a,$1,$13,$17}' a=$evt >>$input.all
### MTTT
   cat $file |grep $cmp |awk '{if($13!=0) print a,$1,$13,$17}' a=$evt >>$input
### XCTT
#   cat $file |awk '{print a,$1,$15,$19}' a=$evt >>$input.all
    cat $file |grep $cmp |awk '{if($13==0&&$15!=0) print a,$1,$15,$19}' a=$evt >>$input
done 
nwin=`cat $input.all |wc -l`
nmeas=`cat $input |wc -l`
if [ $band == "T020_T050" ] ;then
   if [ $i -eq 1 ];then
     gmt psbasemap -JX4.8 -R-10/10/0/4000 -Ba4f1/a1000f500WSen -K -X2 -Y22 -P >$output
   else
     gmt  psbasemap -JX4.8 -R-10/10/0/4000 -Ba4f1/a1000f500WSen -K -O -X7 >>$output
   fi
   ty=1400
   ty1=800
   ty2=500
else
   if [ $i -eq 1 ];then
     gmt psbasemap -JX3.2 -R-5/5/0/300 -Ba4f1/a50f25WSen -K -X2 -Y22 -P >$output
   else
     gmt psbasemap -JX3.2 -R-5/5/0/300 -Ba4f1/a50f25WSen -K -O -X4.8 >>$output
   fi
   ty=350
   ty1=280
   ty2=260
fi
echo $band : $ty
a=`cat $input |awk '{print $3}' |awk 'BEGIN { sum=0;nn=0 } {sum = sum + $1; nn = nn+1} END { printf"%-8.2f\n", sum/nn }'`
b=`cat $input |awk '{print $3}' |awk -v avg=$a 'BEGIN { sum=0;nn=0 } {sum = sum + ($1-avg)*($1-avg); nn = nn+1} END { printf"%-8.2f\n", sqrt(sum/nn) }'`
#cat $input |awk '{print $3}' |pshistogram -JX -R -S -W1.0 -L4.5p/green -F -G255/255/255 -O -K >>$output
#cat $input |awk '{print $3}' |gmt pshistogram -JX -R -S -W1.0 -L1p,0/100/0 -N0+p0.5p,0/100/0 -F -O -K >>$output
#cat $input |awk '{print $3}' |gmt pshistogram -JX -R -S -W0.5 -L1p,100/200/200 -N0+p0.5p,100/200/200 -F -O -K >>$output
cat $input |awk '{print $3}' |gmt pshistogram -JX -R -S -W0.5 -L1p,green -N0+p0.5p,green -F -O -K >>$output

gmt pstext -J -R -O -K -N <<eof >>$output
0.0 $ty 12 0.0 0 CM  ${mod} ${band}
-4.8 $ty1 8 0.0 0 LM  nmeas=$nmeas
-4.8 $ty2 8 0.0 0 LM  Mean:$a std:$b
eof
rm $input $input.all
let i=$i+1
done ## end of band

pwd |gmt psxy -J -R -Sc0.1 -O >>$output
gs $output
