#!/bin/bash
module load gcc/8.3.0  hdf5/1.8.21 netcdf/4.6.3
rm -f gmt.conf .gmtdefaults
gmt gmtset MAP_FRAME_TYPE plain

seisdir=seismograms_obs_syn
bt=-2.
sc=0.20


mkdir -p $seisdir

#IFS=$'\n'
#for evt in `cat ../DATA/inverse_problem/acqui_file.txt`;do
#cat  DATA/inverse_problem/acqui_file.txt |sed -n '13,14p' |awk 'NR % 2 != 0' |
#cat  DATA/inverse_problem/acqui_file.txt  |awk 'NR % 2 != 0' |
for iset in `seq 1 1 1`;do
cat ../src_rec/sources_set$iset.dat |
while read line;do
  ievt=`echo $line |awk '{printf $1}'`
  fwddir=../solver/M00.set$iset/$ievt/OUTPUT_FILES #solver/M00.set1/5/OUTPUT_FILES
  cmtf=../src_rec/CMTSOLUTION_$ievt
  misf=../misfits/M00.set${iset}_T005_T050_window_chi
  echo $ievt
  
  #=======================================================
  hpmag=`head -1 $cmtf |awk '{print $10}'`
  ctd_lat=`grep latitude $cmtf |awk '{print $2}'`
  ctd_lon=`grep longitude $cmtf |awk '{print $2}'`
  ctd_dep=`grep depth $cmtf |awk '{print $2}'`
  mrr=`grep Mrr $cmtf |awk '{print $2}'`
  mtt=`grep Mtt $cmtf |awk '{print $2}'`
  mpp=`grep Mpp $cmtf |awk '{print $2}'`
  mrt=`grep Mrt $cmtf |awk '{print $2}'`
  mrp=`grep Mrp $cmtf |awk '{print $2}'`
  mtp=`grep Mtp $cmtf |awk '{print $2}'`
  #=======================================================
  ###  plot 
  if true;then
  out=$seisdir/fsismo_${ievt}.abs.ps
  icomp=1
  for comp in Z R T;do
    saclst knetwk kstnm dist f $fwddir/*BX${comp}.obs.sac.T005_T050 >dist.dat 
    d_min=`cat dist.dat |minmax -C |awk '{print $7-10}'`
    d_max=`cat dist.dat |minmax -C |awk '{print $8+10}'`
    t_min=0.
    t_max=`echo $d_max |awk '{print $1/1.5}'`
    tx=`echo $t_min $t_max |awk '{print ($1+$2)/2}'`
    ty=`echo $d_min $d_max |awk '{print $2+($2-$1)*0.04}'`

    if [ $icomp -eq 1 ];then
       gmt psbasemap -JX5.5/18 -R$t_min/$t_max/$d_min/$d_max -Ba100f50:"Time (s)":/a50:"Distance (km)":WeSn -K -P -X3 -Y4 >$out
       gmt pstext -J -R -N -O -K >>$out <<eof
$tx $ty 16 0 0 CB Vertical
eof
    elif [ $icomp -eq 2 ];then
       gmt psbasemap -JX5.5/18 -R$t_min/$t_max/$d_min/$d_max -Ba100f50:"Time (s)":/a50weSn -K -O -X6.0 >>$out
       gmt pstext -J -R -N -O -K >>$out <<eof
$tx $ty 16 0 0 CB Radial
eof
    else
       gmt psbasemap -JX5.5/18 -R$t_min/$t_max/$d_min/$d_max -Ba100f50:"Time (s)":/a50weSn -K -O -X6.0 >>$out
       gmt pstext -J -R -N -O -K >>$out <<eof
$tx $ty 16 0 0 CB Transverse
eof
    fi
    # add windows
    cat dist.dat |while read line;do
       stnm=`echo $line |awk '{print $3}'`
       dist=`echo $line |awk '{print $4}'`
       ncmp=`grep $stnm $misf |wc -l`
       if [ $ncmp -gt 0 ];then
          for icmp in `seq 1 1 $ncmp`;do 
            cmp=`grep $stnm $misf |sed -n "${icmp}p" |awk '{print substr($4,3)}'`
            echo $cmop $stnm $ncmp $icmp $cmp
            if [ $cmp == "$comp" ];then
               grep $stnm $misf |awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $7,a-1*m,$7,a+m,$8,a+m,$8,a-1*m}' a=$dist m=3 | gmt psxy -J -R -W0.5 -G220/220/255 -O -K >>$out
            fi
          done
       fi
    done
    # get the absolute amplitude
    cat /dev/null >temp
    cat /dev/null >temp1
    saclst depmax f $fwddir/*BX${comp}.obs.sac.T005_T050 |awk '{print sqrt($2*$2)}' >>temp 
    saclst depmin f $fwddir/*BX${comp}.obs.sac.T005_T050 |awk '{print sqrt($2*$2)}' >>temp 
    saclst depmax f $fwddir/*BX${comp}.syn.sac.T005_T050 |awk '{print sqrt($2*$2)}' >>temp1 
    saclst depmin f $fwddir/*BX${comp}.syn.sac.T005_T050 |awk '{print sqrt($2*$2)}' >>temp1
    #sc1=`cat temp |sort -n |tail -1 |awk '{print 1/$1}'`
    sc1=`cat temp |awk 'BEGIN{sum=0.0;num=0}{sum=sum+$1;num=num+1}END{print num/sum*0.5}'`
    sc2=`cat temp1 |awk 'BEGIN{sum=0.0;num=0}{sum=sum+$1;num=num+1}END{print num/sum*0.5}'`
    rm temp temp1
    echo "$comp average amplitudes (obs, syn) =  " $sc1 $sc2
    gmt pssac $fwddir/*BX${comp}.obs.sac.T005_T050 -J -R  -Ektb -M50000/0 -W0.5p -K -O >>$out 
    gmt pssac $fwddir/*BX${comp}.syn.sac.T005_T050 -J -R  -Ektb -M50000/0 -W0.5p,red -K -O >>$out 
    cat dist.dat |awk '{printf"%f %f 8 0 0 CM %s.%s\n",-10,$4,$2,$3}' |gmt pstext -J -R -G255/255/255 -N -K -O >>$out 
    let icomp=icomp+1
  done
  #### tectonic map 
  gmt psbasemap -R-122./-117./34.5/37.7 -JM2.8i -Ba4f1/a4f1WesN -X-8. -Y20 -O -K >>$out
  gmt pscoast -J -R  -O -K -N1 -W1 -Dl >> $out
  ### add ray paths
  saclst evlo evla stlo stla f $fwddir/*BX${comp}.obs.sac.T005_T050 |awk '{printf">\n%f %f\n%f %f\n",$2,$3,$4,$5}' |gmt psxy -J -R -W0.5p -O -K >>$out
  ### add station symbols
  saclst stlo stla f $fwddir/*BX${comp}.obs.sac.T005_T050 |awk '{print $2,$3}' |gmt psxy -J -R -St0.05i -Gblue -O -K >>$out
  ### add Focal mechanism 
  gmt psmeca -J -R -Sd0.3i -Gred -L0.5p,black -W1p,red -N -O -K  >> $out <<END
$ctd_lon, $ctd_lat, $ctd_dep, $mrr, $mtt, $mpp, $mrt, $mrp, $mtp, $hpmag 0, 0, ${ievt}A
END
  cat /dev/null |gmt psxy -J -R -Sc0.1i -O >>$out
  ### finish plot
  fi

done # end source 
done # end iset
