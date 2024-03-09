#!/bin/bash

#module load gcc hdf5/1.8.20 netcdf/4.6.1
module load gcc/4.9.1 netcdf/4.4.1
gmtset ANNOT_FONT_SIZE_PRIMARY 12p
gmtset ANNOT_OFFSET_PRIMARY 0.1c
gmtset LABEL_FONT_SIZE 14p 
gmtset LABEL_OFFSET 0.15c
gmtset HEADER_FONT_SIZE 16p
gmtset HEADER_OFFSET -0.5c 
gmtset TICK_LENGTH -0.2c
gmtset BASEMAP_TYPE plain 

seisdir=seismograms_with_stf
prog=~/progs/teleseis_stf
bt=-10.
sc=0.20 ## plot normalized trace


mkdir -p $seisdir

#IFS=$'\n'
#for evt in `cat ../DATA/inverse_problem/acqui_file.txt`;do
#cat  DATA/inverse_problem/acqui_file.txt |sed -n '13,14p' |awk 'NR % 2 != 0' |
#cat  DATA/inverse_problem/acqui_file.txt  |awk 'NR % 2 != 0' |
for iset in `seq 1 1 13`;do
cat ../src_rec/sources_set$iset.dat |
while read line;do
  ievt=`echo $line |awk '{printf $1}'`
  fwddir=../solver/M00.set$iset/$ievt/OUTPUT_FILES #solver/M00.set1/5/OUTPUT_FILES
  #mkdir -p $seisdir/$ievt
  #saclst kstnm knetwk kcmpnm t0 f ../data/$ievt/*BXZ.sac |awk '{print $2,$3,$4,$5}' >$seisdir/$ievt/FKtimes.dat
  phasetab=../src_rec/FKtimes_${ievt}
  twb=5 #`head -1 $phasetab |awk '{print $4}'`
  twe=45 #`head -1 $phasetab |awk '{print $5}'`
  echo $ievt $twb $twe
  #=======================================================

  for comp in Z R ;do
    ls $fwddir/dat.*BX${comp}*.sac.T005_T050 >file1.lst
    cat file1.lst |sed -e "s/dat/syn/g" >file2.lst
    ls $fwddir/*BX.sac.T005_T050.dec >file3.lst
    ls $fwddir/PCs00?.sac        >file4.lst
    cat /dev/null >file5.lst
    for dfile in ` cat file3.lst`;do
       stnm=`echo $dfile |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'` 
       pdm=`grep -n $stnm ../src_rec/STATIONS_${ievt}_FILTERED |awk -F: '{printf"PDMs%03d",$1}'`
       echo $fwddir/${pdm}.sac    >>file5.lst
    done
    cat /dev/null >file6.lst
    for dfile in ` cat file1.lst`;do
       net=`echo $dfile |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'` 
       stnm=`echo $dfile |awk -F$fwddir/ '{print $2}' |awk -F. '{print $3}'` 
       echo $fwddir/$net.$stnm.BX.sac.T005_T050.rec${comp}  >>file6.lst
    done

    if [ $comp == 'R' ];then
       ls $fwddir/../SEM/*BXE*.sac >file7.lst
    else
       ls $fwddir/../SEM/*BX${comp}*.sac >file7.lst
    fi
    ls $fwddir/../SEM/*BX${comp}*.sac.T005_T050 >file8.lst

    nsta=`cat file2.lst |wc -l`
    len=20.
    npart=`echo $nsta $len |awk '{print $1/$2}' | awk '{printf("%d\n",$0+=$0<0?0:0.9999)}'`
    echo "nsta,len,npart=" $nsta $len $npart
    for i in `seq $npart`;do
      ib=`echo $i $len |awk '{print ($1-1)*$2+1}'`
      ie=`echo $i $len |awk '{print $1*$2}'`
      if [ ${ie} -gt $nsta ];then
          ie=$nsta
      fi
      ymax=`echo $len |awk '{print $1+1}'`
      echo $comp $ib $ie $ymax
    ##########################################
    if false;then
    out=$seisdir/fsismo_${ievt}_BX${comp}.seg${ib}-${ie}.ps
    psbasemap -R-10/120/-1/$ymax -JX8/11 -Ba20f10/a1:"Trace number":weSn  -P -K -X3 -Y15 > $out
    cat file1.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p -O -K >>$out
    cat file2.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p,255/0/0 -O -K >>$out
    ty=`echo $ymax |awk '{print $1*1.05}'`
    ty1=`echo $ymax |awk '{print $1*1.15}'`
    pstext -J -R -N -K -O >>$seisdir/fsismo_${ievt}_BX${comp}.seg${ib}-${ie}.ps <<EOF
45 $ty 14 0 0 CM (a) DATA and SYN without STF 
85 $ty1 16 0 0 CM EventID: $ievt
EOF
    cat /dev/null >temp
    cat /dev/null >temp1
    for fname in `cat file1.lst`;do
      netwk=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
      stnm=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $3}'`
      echo $netwk $stnm >>temp 
      grep $stnm $phasetab >>temp1 
    done
    cat temp |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",-10,NR,$1,$2}' |pstext -J -R -W255/0/0 -N -K -O >>$out 
    ## (1) P
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t-b,NR-0.25,0,0.2}' t=$bt b=$twb | psxy -J -R -SV0.05/0./0. -Gmagenta  -N -O -K >>$out
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t+e,NR-0.25,0,0.2}' t=$bt e=$twe | psxy -J -R -SV0.05/0./0. -Gmagenta  -N -O -K >>$out
    rm temp temp1  
    ### plot RF function
    psbasemap -R-10/120/-1/$ymax -JX8/11 -Ba20f10/a1weSn -O -K -X9 >> $out
    cat file3.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p -O -K >> $out
    cat file5.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p,red -O -K >> $out
    ty=`echo $ymax |awk '{print $1*1.05}'`
    pstext -J -R -N -O -K >>$out <<EOF
55 $ty 16 0 0 CM (b) Deconvolved vertical traces 
EOF
    cat /dev/null >temp
    for fname in `cat file3.lst`;do
      netwk=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $1}'`
      stnm=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
      echo $netwk $stnm >>temp
    done
    cat temp |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",-10,NR,$1,$2}' |pstext -J -R -W255/0/0 -N -K -O >>$out
    rm temp 
    ### plot PC analysis
    psbasemap -R-10/120/0.0/2 -JX8/2 -Ba20f10:"Time (s)":/a10:"Amplitude":WeSn -X-9 -Y-12.5 -K -O >>$out
    #pssac $fwddir/stf_pca001_Z.sac -J -R -Entb -M0.5 -W1p,red -O -K >>$out
    pssac $fwddir/PCs001.sac -J -R -Entb -M0.5 -W1p,red -O -K >>$out
    pstext -J -R -O -K >>$out <<EOF
105 0.25 12 0 0 CM STF
EOF

    xmax=`cat $fwddir/contribution.dat |wc -l |awk '{print $1+1}'`
    tx=`echo $xmax |awk '{print $1/2}'`
    psbasemap -R0/$xmax/-5/100 -JX8/3 -Ba5/a20f10:,-%:WeSn -Y3 -K -O >>$out
    cat $fwddir/contribution.dat |awk '{print NR,$1}' |
    psxy -J -R  -W1p -O -K >>$out
    cat $fwddir/contribution.dat |awk '{print NR,$1}' |
    psxy -J -R -Sd0.1i -W0.5p -Gblue -O -K >>$out
    pstext -J -R -O -K >>$out <<EOF
$tx 80 12 0 0 CM PC contribution
EOF


    psbasemap -R-10/120/-1/10 -JX8/4 -Ba20f10/a1:"PCs number":WeSn  -O -K -Y4 >>$out
    cat file4.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p -O -K >> $out 
    head -1 file4.lst | pssac -J -R -Entb -M$sc -W1p,red -O -K >> $out
    pstext -J -R -O -K -N >>$out <<EOF
45 0 12 0 0 CM PC time series
55 11.3 16 0 0 CM (c) PCA analysis
EOF

    ### plot reconstruted waveform
    cat /dev/null >temp
    saclst depmax f $fwddir/dat.*BX${comp}*.sac.T005_T050 |awk '{print sqrt($2*$2)}' >>temp 
    saclst depmin f $fwddir/dat.*BX${comp}*.sac.T005_T050 |awk '{print sqrt($2*$2)}' >>temp 
    #sc1=`cat temp |sort -n |tail -1 |awk '{print 1/$1}'`
    sc1=`cat temp |awk 'BEGIN{sum=0.0;num=0}{sum=sum+$1;num=num+1}END{print num/sum*0.5}'`
    echo "sc1: " $sc1
    psbasemap -R-10/120/0/$ymax -JX8/11 -Ba20f10:"Time (s)":/a1weSn  -O -K -X9 -Y-7 >>$out
    cat file1.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc1/0 -W1p -O -K >> $out 
    cat file6.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc1/0 -W1p,red -O -K >> $out 
    ty=`echo $ymax |awk '{print $1*1.05}'`
    pstext -J -R -N -O -K >>$out <<EOF
45 $ty 16 0 0 CM (d) DATA and SYN with STF 
EOF
    cat /dev/null >temp
    cat /dev/null >temp1
    for fname in `cat file1.lst`;do
      netwk=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
      stnm=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $3}'`
      echo $netwk $stnm >>temp 
      grep $stnm $phasetab >>temp1 
    done
    cat temp |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",-10,NR,$1,$2}' |pstext -J -R -W255/0/0 -N -K -O >>$out 
    ## (1) P
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t-b,NR-0.25,0,0.2}' t=$bt b=$twb | psxy -J -R -SV0.05/0./0. -Gmagenta  -N -O -K >>$out
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t+e,NR-0.25,0,0.2}' t=$bt e=$twe | psxy -J -R -SV0.05/0./0. -Gmagenta  -N -O -K >>$out
    rm temp temp1
    ## Get amplitude scale factors 
#    if true;then
#    cat /dev/null >ampl_fac_${ievt}_${comp}.dat
#
#    for dat in `ls $fwddir/*.num`;do
#       nt=`echo $dat |awk -F$fwddir/ '{print $2}' |awk -F. '{print $1}'`
#       stnm=`echo $dat |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
#       syn=$fwddir/$nt.$stnm.BX.sac.T005_T050.recp${comp}
#       pred_ttp=`grep $stnm $phasetab |awk '{print $3}'`
#       t0=`echo $pred_ttp |awk '{print $1+t-b}' t=$bt b=$twb`
#       t1=`echo $pred_ttp |awk '{print $1+t+e}' t=$bt e=$twe`
#       echo $dat $syn $pred_ttp $t0 $t1
#       sac <<!
#r $dat
#mulf $syn
#mul 1000000000000
#abs
#w s1.sac
#r $syn
#mulf $syn
#mul 1000000000000
#w s2.sac
#quit
#!
#       ~/progs/sactools_c/sac_ave $t0 $t1 s1.sac >temp
#       ave1=`tail -1 temp |awk '{print $2}'`
#       ~/progs/sactools_c/sac_ave $t0 $t1 s2.sac >temp
#       ave2=`tail -1 temp |awk '{print $2}'`
#       ampl_fac=`echo $ave1 $ave2 |awk '{print $1/$2}'`
#       echo $stnm $ampl_fac >>ampl_fac_${ievt}_${comp}.dat
#    done
#    rm temp
#    fi
#    cat ampl_fac_${ievt}_${comp}.dat |sed -n "${ib},${ie}p" |awk '{printf"110 %d 8 0 0 CB %f\n",NR,$2}' |pstext -J -R -Gred -N -O -K >>$out
    cat /dev/null >temp
    for fname in `cat file1.lst`;do
      netwk=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
      stnm=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $3}'`
      grep $stnm $fwddir/ampl_fac_${comp}.dat >>temp 
    done
    cat temp |sed -n "${ib},${ie}p" |awk '{printf"110 %f 8 0 0 CB %f\n",NR+0.2,$2}' |pstext -J -R -G255/0/0 -N -K -O >>$out 
    rm temp

    cat /dev/null |psxy -J -R -Sc0.1i -O >>$out
    ##
    fi ### end if plot stf
    #######################################
    if false;then
    out1=$seisdir/adjsource_${ievt}_BX${comp}.seg${ib}-${ie}.ps
    ### Adjoint sources
    psbasemap -R-10/120/0/$ymax -JX8/11 -Ba20f10:"Time (s)":/a1:"Trace number":weSn  -P -K -X4 -Y2 > $out1
    cat file7.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p,255/0/0 -O -K >>$out1
    ty=`echo $ymax |awk '{print $1*1.05}'`
    pstext -J -R -N -O -K >>$out1 <<EOF
45 $ty 16 0 0 CM (b) Adjoint sources
EOF
    for fname in `cat file1.lst`;do
      netwk=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
      stnm=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $3}'`
      echo $netwk $stnm >>temp 
    done
    cat temp |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",-10,NR,$1,$2}' |pstext -J -R -W255/0/0 -N -K -O >>$out1 
    rm temp

    ### Residuals
    psbasemap -R-10/120/0/$ymax -JX8/11 -Ba20f10:"Time (s)":/a1:"Trace number":wesn  -O -K  -Y13.5 >> $out1
    cat file8.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p,255/0/0 -O -K >>$out1
    ty=`echo $ymax |awk '{print $1*1.05}'`
    ty1=`echo $ymax |awk '{print $1*1.15}'`
    pstext -J -R -N -O -K >>$out1 <<EOF
45 $ty 16 0 0 CM (a) Residuals
45 $ty1 16 0 0 CM EventID: $ievt
EOF
    for fname in `cat file1.lst`;do
      netwk=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
      stnm=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $3}'`
      echo $netwk $stnm >>temp 
      grep $stnm $phasetab >>temp1 
    done
    cat temp |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",-10,NR,$1,$2}' |pstext -J -R -W255/0/0 -N -O -K >>$out1 
    ## (1) P
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t-b,NR-0.25,0,0.2}' t=$bt b=$twb | psxy -J -R -SV0.05/0./0. -Gmagenta  -N -O -K >>$out1
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t+e,NR-0.25,0,0.2}' t=$bt e=$twe | psxy -J -R -SV0.05/0./0. -Gmagenta  -N -O >>$out1
    rm temp temp1
    fi ### end if plot adj

    #########################################
    done # end ipart
    #######################################
    ### plot ampl_factor
    if true;then
    out2=$seisdir/sta_ampl_${ievt}_${comp}.ps 
    cat /dev/null >temp
    for fname in `cat file1.lst`;do
      netwk=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
      stnm=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $3}'`
      ampl_fac=`grep $stnm $fwddir/ampl_fac_${comp}.dat |awk '{print $2}'` 
      stla=`grep $stnm ../src_rec/STATIONS_${ievt} |awk '{print $3}'` 
      stlo=`grep $stnm ../src_rec/STATIONS_${ievt} |awk '{print $4}'` 
      echo $stlo $stla $ampl_fac >>temp 
    done
    makecpt -Csvel13.cpt -T0/0.2/0.01 -D > seis.cpt
    psbasemap -JM8 -R-101.5/-93.7/15.3/22.5 -Ba2/a2 -X8 -Y12 -K -P >$out2
    cat temp |psxy -J -R -Sc0.04i -Cseis.cpt -W0.5p,black -O -K >>$out2
    psscale -D4.0/-2/6/0.3h -Cseis.cpt -B0.1 -E -O -N >> $out2
    rm temp
    fi

    rm file?.lst
  done # end comp
done  # end source_set?.dat
done # end iset
