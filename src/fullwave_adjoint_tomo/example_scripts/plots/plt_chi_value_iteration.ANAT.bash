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
out=chi_iter_alter_ANAT_${cmp}.ps
range=-R-1./36/0.65/1.5
xtick=a4f2
ytick=a0.2
xc=`echo $range |awk '{print substr($1,3)}' |awk -F/ '{print ($2-$1)/2}'`
yc=`echo $range |awk '{print substr($1,3)}' |awk -F/ '{print $4*1.05}'`
if true;then
#for band in T006_T015 T010_T020 T015_T030 T020_T040;do
for band in T006_T015 T010_T020 T015_T035;do
    echo $band
    cat /dev/null >chi_$band.$cmp.dat
    for itnm in `seq 0 2 32`;do
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
#########
cat /dev/null  >chi_overall_band.$cmp.dat
l=0
for i in `seq 0 2 32`;do
  let l=$l+1
  #chi1=`cat chi_T020_T040.$cmp.dat |sed -n "${l}p" |awk '{print $2}'`
  #nc1=`cat chi_T020_T040.$cmp.dat |sed -n "${l}p" |awk '{print $3}'`
  chi2=`cat chi_T015_T035.$cmp.dat |sed -n "${l}p" |awk '{print $2}'`
  nc2=`cat chi_T015_T035.$cmp.dat |sed -n "${l}p" |awk '{print $3}'`
  chi3=`cat chi_T010_T020.$cmp.dat |sed -n "${l}p" |awk '{print $2}'`
  nc3=`cat chi_T010_T020.$cmp.dat |sed -n "${l}p" |awk '{print $3}'`
  chi4=`cat chi_T006_T015.$cmp.dat |sed -n "${l}p" |awk '{print $2}'`
  nc4=`cat chi_T006_T015.$cmp.dat |sed -n "${l}p" |awk '{print $3}'`
  #echo $i $chi1 $chi2 $chi3 $chi4 $nc1 $nc2 $nc3 $nc4 |awk '{print $1,($2+$3+$4+$5)/($6+$7+$8+$9)}' >>chi_overall_band.$cmp.dat
  echo $i $chi2 $chi3 $chi4 $nc2 $nc3 $nc4 |awk '{print $1,($2+$3+$4)/($5+$6+$7)}' >>chi_overall_band.$cmp.dat
  #echo $i $chi1 $chi2 $chi3 $chi4 |awk '{print $1,($2+$3+$4)/3}' >>chi_overall_band.$cmp.dat
done
psbasemap -JX12/4 $range -B"$xtick":"Iter number":/"$ytick":"Misfits":WSne -K -P -X4 -Y16 >$out
cat chi_overall_band.$cmp.dat |psxy -J -R -W1p,red -O -K >>$out
cat chi_overall_band.$cmp.dat |psxy -J -R -Sa0.08i -Gred -Wblack -O -K >>$out
cat chi_ANAT.dat |awk '{print $1*2,$2}' |psxy -J -R -W1p,blue,- -O -K >>$out
cat chi_ANAT.dat |awk '{print $1*2,$2}' |psxy -J -R -Ss0.08i -Gblue -Wblack -O >>$out
