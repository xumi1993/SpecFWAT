#!/bin/bash
#        ! KEY: write misfit function values to file (two for each window)
#        ! Here are the 20 columns of the vector window_chi
#        !  1: MT-TT chi,    2: MT-dlnA chi,    3: XC-TT chi,    4: XC-dlnA chi
#        !  5: MT-TT meas,   6: MT-dlnA meas,   7: XC-TT meas,   8: XC-dlnA meas
#        !  9: MT-TT error, 10: MT-dlnA error, 11: XC-TT error, 12: XC-dlnA error
#        ! WINDOW     : 13: data power, 14: syn power, 15: (data-syn) power, 16: window duration
#        ! FULL RECORD: 17: data power, 18: syn power, 19: (data-syn) power, 20: record duration
#        ! Example of a reduced file: awk '{print $2,$3,$4,$5,$6,$31,$32}' window_chi > window_chi_sub
rm .gmtdefaults4 
#gmtset ANOT_FONT_SIZE_PRIMARY 10
gmtset LABEL_FONT_SIZE 18

scratch=..
cmp=BX
out=chi_iter_alter_TeleFWI_${cmp}.ps
range=-R-2./36/0.3/0.7
xtick=a4f2
ytick=a0.2
xc=`echo $range |awk '{print substr($1,3)}' |awk -F/ '{print ($2-$1)/2}'`
yc=`echo $range |awk '{print substr($1,3)}' |awk -F/ '{print $4*1.05}'`
if true;then
#for band in T006_T015 T010_T020 T015_T030 T020_T040;do
#for band in T006_T015 T010_T020 T015_T030;do
for band in T005_T050;do
    echo $band
    cat /dev/null >chi_$band.$cmp.dat
    for itnm in `seq 1 2 33`;do
       mod=`echo $itnm |awk '{printf"M%02d",$1}'`
       if [ $itnm -eq 0 ];then
          misf=`cat $scratch/misfits/${mod}.set*_${band}_*chi |grep $cmp |awk 'BEGIN{sum=0;n=0} {if($29!=0) {sum=sum+$29;n=n+1}} END{print sum}'`
          misn=`cat $scratch/misfits/${mod}.set*_${band}_*chi |grep $cmp |awk 'BEGIN{sum=0;n=0} {if($29!=0) {sum=sum+$29;n=n+1}} END{print n}'`
          #misf=`cat ../../output/misfits/$mod/${mod}_${band}_*_chi |awk 'BEGIN{sum=0;n=0} {if($13!=0) {sum=sum+$13*$13/$17/$17;n=n+1}} END{print 0.5*sum/n}'`
          #misf=`cat ../../output/misfits/$mod/${mod}_${band}_*_chi |awk 'BEGIN{sum=0;n=0} {if($15!=0) {sum=sum+$15*$15/$19/$19;n=n+1}} END{print 0.5*sum/n}'`
       else
          misf=`cat $scratch/misfits/${mod}.set*_${band}_*chi |grep $cmp |awk 'BEGIN{sum=0;n=0} {if($29!=0) {sum=sum+$29;n=n+1}} END{print sum}'`
          misn=`cat $scratch/misfits/${mod}.set*_${band}_*chi |grep $cmp |awk 'BEGIN{sum=0;n=0} {if($29!=0) {sum=sum+$29;n=n+1}} END{print n}'`
          #misf=`cat ../../output/misfits/$mod/${mod}.set*_${band}_*_chi |awk 'BEGIN{sum=0;n=0} {if($13!=0) {sum=sum+$13*$13/$17/$17;n=n+1}} END{print 0.5*sum/n}'`
          #misf=`cat ../../output/misfits/$mod/${mod}.set*_${band}_*_chi |awk 'BEGIN{sum=0;n=0} {if($15!=0) {sum=sum+$15*$15/$19/$19;n=n+1}} END{print 0.5*sum/n}'`
       fi
    echo $itnm $misf $misn  >>chi_$band.${cmp}.dat 
    done
done
fi
psbasemap -JX12/4 $range -B"$xtick":"Iter number":/"$ytick":"Misfits":WSne -K -P -X4 -Y16 >$out
cat chi_T005_T050.${cmp}.dat  |awk '{printf"%d %f\n",$1,$2/$3}' |psxy -J -R -W1p,red -O -K >>$out
cat chi_T005_T050.${cmp}.dat  |awk '{printf"%d %f\n",$1,$2/$3}' |psxy -J -R -Sa0.08i -Gred -Wblack -O -K >>$out
cat chi_TeleFWI.dat |awk '{printf"%d %f\n",$1*2+1,$2/$3}' |psxy -J -R -W1p,blue,- -O -K >>$out
cat chi_TeleFWI.dat |awk '{printf"%d %f\n",$1*2+1,$2/$3}' |psxy -J -R -St0.08i -Gblue -Wblack -O -K >>$out
### plot legend
pslegend -J -R -D24/0.638/3.1/2.0/LT -O -K <<eof >>$out
S 0.45c - 0.3i - 0.5p,red 0.6i
S 0.45c - 0.3i - 0.5p,- 0.6i
S 0.45c - 0.3i - 0.5p,- 0.6i
eof

pslegend -J -R -D24/0.638/3.1/2.0/LT -F -O <<eof >>$out
S 0.45c a 0.08i red black 1.c Alter
S 0.45c t 0.08i blue black 1.c TeleFWI
S 0.45c s 0.08i blue black 1.c ANAT
eof

