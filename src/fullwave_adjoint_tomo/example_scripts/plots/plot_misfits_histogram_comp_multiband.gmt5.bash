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
#gmtset ANNOT_FONT_SIZE_PRIMARY 12p
#gmtset LABEL_FONT_SIZE 12p
#gmtset HEADER_FONT_SIZE 36p
#gmtset ANNOT_FONT_SIZE_SECONDARY 12p
gmt set FONT_ANNOT_PRIMARY 12p
gmt set FONT_LABEL 14p

mod=$1
mod1=$2
cmp=$3
scratch=/mnt/scratch-lustre/wangkai/CeCal_fwat_alter_realdata

output=dT_${mod}_${mod1}_${cmp}_multiband.ps

misf=dT_${mod}_all.dat
misf1=dT_${mod1}_all.dat

cat /dev/null >$misf
cat /dev/null >$misf1
i=1
#for band in T006_T015 T010_T020 T015_T030 T020_T040;do
for band in T006_T015 T010_T020 T015_T035;do
####################  for mod ========================
input=dT_${mod}_${band}_${cmp}.dat
cat /dev/null >$input
cat /dev/null >$input.all
ls $scratch/misfits/${mod}.set*_${band}_*_chi >chi.dat
for file in `cat chi.dat`;do
   echo $file
   evt=`echo $file |awk -F$scratch/misfits/ '{print $2}' |awk -F_ '{print $4 }'`

   cat $file |grep $cmp |awk '{print a,$1,$13,$17}' a=$evt >>$input.all
### MTTT
   cat $file |grep $cmp |awk '{if($13!=0) print a,$1,$13,$17}' a=$evt >>$input
### XCTT
#   cat $file |awk '{print a,$1,$15,$19}' a=$evt >>$input.all
    cat $file |grep $cmp  |awk '{if($13==0&&$15!=0) print a,$1,$15,$19}' a=$evt >>$input
done 
nwin=`cat $input.all |wc -l`
nmeas=`cat $input |wc -l`
cat $input >>$misf
###################### for mod1 #############################
input1=dT_${mod1}_${band}_${cmp}.dat
cat /dev/null >$input1
cat /dev/null >$input1.all
if [ $mod1 == "M18" ]  && [ $band == "T005_T050" ];then
  ls $scratch/misfits/$mod1/${mod1}.set*_${band}_*_chi >chi.dat
else
  ls $scratch/misfits/${mod1}.set*_${band}_*_chi >chi.dat
fi
for file in `cat chi.dat`;do
   echo $file
   evt=`echo $file |awk -F$scratch/misfits/ '{print $2}' |awk -F_ '{print $4 }'`

   cat $file |grep $cmp |awk '{print a,$1,$13,$17}' a=$evt >>$input1.all
### MTTT
   cat $file |grep $cmp |awk '{if($13!=0) print a,$1,$13,$17}' a=$evt >>$input1
### XCTT
#   cat $file |awk '{print a,$1,$15,$19}' a=$evt >>$input.all
    cat $file |grep $cmp |awk '{if($13==0&&$15!=0) print a,$1,$15,$19}' a=$evt >>$input1
done 
nwin1=`cat $input1.all |wc -l`
nmeas1=`cat $input1 |wc -l`
cat $input1 >>$misf1
##############

if [ $band == "T020_T050" ] ;then
   if [ $i -eq 1 ];then
     gmt psbasemap -JX4.8 -R-10/10/0/4000 -Ba4f1/a1000f500:"Num. of Windows":WSen -K -X2 -Y22 -P >$output
   else
     gmt psbasemap -JX4.8 -R-10/10/0/4000 -Ba4f1/a1000f500WSen -K -O -X7 >>$output
   fi
   ty=4400
   ty1=3800
   ty2=3500
else
    if [ $i -eq 1 ];then
     gmt psbasemap -JX3.2 -R-6/6/0/600 -Ba4f1:"dT (sec)":/a100f50:"Num. of Windows":WSen -K -X2 -Y22 -P >$output
   else
     gmt psbasemap -JX3.2 -R-6/6/0/600 -Ba4f1:"dT (sec)":/a100f50wSen -K -O -X4.0 >>$output
   fi
   ty=660
   ty1=540
   ty2=480
fi
### plot 1
a=`cat $input |awk '{print $3}' |awk 'BEGIN { sum=0;nn=0 } {sum = sum + $1; nn = nn+1} END { printf"%-4.2f\n", sum/nn }'`
b=`cat $input |awk '{print $3}' |awk -v avg=$a 'BEGIN { sum=0;nn=0 } {sum = sum + ($1-avg)*($1-avg); nn = nn+1} END { printf"%-8.2f\n", sqrt(sum/nn) }'`
#cat $input |awk '{print $3}' |pshistogram -JX -R -S -W1.0 -L2.5p/green -F -G255/255/255 -O -K >>$output
cat $input |awk '{print $3}' |gmt pshistogram -JX -R -S -W0.5 -L1p,green -N1+p0.5p,green  -F -O -K >>$output

gmt pstext -J -R -O -K -N <<eof >>$output
0.0 $ty 12 0.0 0 CM  $band
eof
gmt pstext -J -R -Ggreen -O -K -N <<eof >>$output
-5.8 $ty1 8 0.0 0 LM $mod: $nmeas,${a}\261${b}
eof
### plot 2
a=`cat $input1 |awk '{print $3}' |awk 'BEGIN { sum=0;nn=0 } {sum = sum + $1; nn = nn+1} END { printf"%-4.2f\n", sum/nn }'`
b=`cat $input1 |awk '{print $3}' |awk -v avg=$a 'BEGIN { sum=0;nn=0 } {sum = sum + ($1-avg)*($1-avg); nn = nn+1} END { printf"%-8.2f\n", sqrt(sum/nn) }'`
cat $input1 |awk '{print $3}' |gmt pshistogram -JX -R -S -W0.5 -L1p,red -N1+p0.5p,red -F -G255/255/255 -O -K >>$output
#cat $input1 |awk '{print $3}' |pshistogram -JX -R  -W1.0 -L0.5p/red -F -G255/255/255 -O -K >>$output

gmt pstext -J -R -O -K -Gred -N <<eof >>$output
-5.8 $ty2 8 0.0 0 LM $mod1: $nmeas1,$a\261$b
eof

rm chi.dat $input $input1 $input.all $input1.all

let i=$i+1
done ## end of band
####
nmeas=`cat $misf |wc -l`
nmeas1=`cat $misf1 |wc -l`
ty=1650
ty1=1350
ty2=1200
gmt psbasemap -JX -R-6/6/0/1500 -Ba4f1:"dT (sec)":/200f100wSEn -K -O -X4.0 >>$output
a=`cat $misf |awk '{print $3}' |awk 'BEGIN { sum=0;nn=0 } {sum = sum + $1; nn = nn+1} END { printf"%-4.2f\n", sum/nn }'`
b=`cat $misf |awk '{print $3}' |awk -v avg=$a 'BEGIN { sum=0;nn=0 } {sum = sum + ($1-avg)*($1-avg); nn = nn+1} END { printf"%-8.2f\n", sqrt(sum/nn) }'`
cat $misf |awk '{print $3}' |gmt pshistogram -JX -R -S -W0.5 -L1p,green -N1+p0.5p,green -F -G255/255/255 -O -K >>$output
gmt pstext -J -R -O -K -N <<eof >>$output
0.0 $ty 12 0.0 0 CM  Over all
eof
gmt pstext -J -R -O -K -Ggreen -N <<eof >>$output
-5.8 $ty1 8 0.0 0 LM ${mod}: $nmeas,$a\261$b
eof
a=`cat $misf1 |awk '{print $3}' |awk 'BEGIN { sum=0;nn=0 } {sum = sum + $1; nn = nn+1} END { printf"%-4.2f\n", sum/nn }'`
b=`cat $misf1 |awk '{print $3}' |awk -v avg=$a 'BEGIN { sum=0;nn=0 } {sum = sum + ($1-avg)*($1-avg); nn = nn+1} END { printf"%-8.2f\n", sqrt(sum/nn) }'`
cat $misf1 |awk '{print $3}' |gmt pshistogram -JX -R -S -W0.5 -L1p,red -N1+p0.5p,red -F -G255/255/255 -O -K >>$output
gmt pstext -J -R -O -K -Gred -N <<eof >>$output
-5.8 $ty2 8 0.0 0 LM $mod1: $nmeas1,$a\261$b
eof

####

pwd |gmt psxy -J -R -Sc0.1 -O >>$output
gs $output
