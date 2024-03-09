#!/bin/bash

module load gcc #hdf5/1.8.20 netcdf/4.6.1

seisdir=seismograms_telestf
#seisdir=seismograms_with_stf
#fktime=FKtimes.dat
prog=~/progs/teleseis_stf
is_sac2asc=true
bt=-10.


mkdir -p $seisdir

#IFS=$'\n'
#for evt in `cat ../DATA/inverse_problem/acqui_file.txt`;do
#cat  DATA/inverse_problem/acqui_file.txt |sed -n '13,14p' |awk 'NR % 2 != 0' |
#cat  DATA/inverse_problem/acqui_file.txt  |awk 'NR % 2 != 0' |
head -1 src_rec/sources_set1.dat |
while read line;do
  ievt=`echo $line |awk '{printf $1}'`
  fwddir=solver/M00.set1/$ievt/OUTPUT_FILES #solver/M00.set1/5/OUTPUT_FILES
  phasetab=$seisdir/$ievt/FKtimes.dat
  echo $ievt
  twb=5 #`grep "time_window true" DATA/DATA_SYNTHS/$ievt/FK$ievt |awk '{print $3}'`
  twe=45 #`grep "time_window true" DATA/DATA_SYNTHS/$ievt/FK$ievt |awk '{print $4}'`
  #=======================================================
  if $is_sac2asc;then
  ### Convert binary to ascii ###
  mkdir -p $seisdir/$ievt
  saclst kstnm knetwk kcmpnm t0 f data/$ievt/*BHZ.sac |awk '{print $2,$3,$4,$5}' >$seisdir/$ievt/FKtimes.dat
  cd $seisdir/$ievt
  ln -sf ../../$fwddir/*.T005_T050 . 
  #===================================================
  ### Apply time-domain iterative deconvolution
  module load intel/2018.2
  ls dat*BHZ*.T005_T050 >datZ.lst
  ls dat*BHZ*.T005_T050 |sed -e "s/BHZ/BHR/g" >datR.lst
  paste datZ.lst datR.lst >dat.lst
  ls syn*BHZ*.T005_T050 >synZ.lst
  ls syn*BHZ*.T005_T050 |sed -e "s/BHZ/BHR/g">synR.lst
  paste synZ.lst synR.lst >syn.lst
  cat FKtimes.dat |awk '{print $4}'>ttp.lst
  $prog/get_telestf dat.lst syn.lst ttp.lst
  cp PCs001.sac stf_pca001.sac
  cd ../..
  #rm dat*.lst syni*.lst
  exit
  ### Reconstruct syntheic using average STF 
  ista=1
  cat /dev/null >cost_valueR.dat
  cat /dev/null >cost_valueZ.dat
  cat filelist |
  while read line;do
    netwk=`echo $dat |awk -F. '{print $2}'`
    stnm=`echo $dat |awk -F. '{print $3}'`
    chan=`echo $dat |awk -F. '{print $4}'`
    syn=`echo $line |awk -F$seisdir/$ievt/dec '{printf"syn%s\n",$2}'`
    synR=`echo $syn |sed -e "s/BHZ/BHR/g"`
    dat=`echo $syn |awk -Fsyn '{printf"dat%s\n",$2}'` 
    datR=`echo $dat |sed -e "s/BHZ/BHR/g"`
    #synT=`echo $syn |sed -e "s/BHZ/BHT/g"`
    #rawsyn=`echo $syn |awk -F. '{printf"%s.%s.%s",$1,$2,$3}'`
    #rawsynR=`echo $rawsyn |sed -e "s/BHZ/BHR/g"`
    #rawdat=`echo $rawsyn |awk -Fsyn '{printf"dat%s\n",$2}'` 
    #rawdatR=`echo $rawdat |sed -e "s/BHZ/BHR/g"`
    rec=`echo $line |awk -F$seisdir/$ievt/dec '{printf"rec%s\n",$2}'`
    recR=`echo $rec |sed -e "s/BHZ/BHR/g"`
    recT=`echo $rec |sed -e "s/BHZ/BHT/g"`
    adj=`echo $line |awk -F$seisdir/$ievt/dec '{printf"adj%s\n",$2}'`
    adjR=`echo $adj |sed -e "s/BHZ/BHR/g"`
    res=`echo $line |awk -F$seisdir/$ievt/dec '{printf"res%s\n",$2}'`
    resR=`echo $res |sed -e "s/BHZ/BHR/g"`
    lnum=`echo $stalst $ista |awk '{print $1+$2-1}'`
    is_good=1 #`cat $evtdir/FK$ievt |sed -n "${lnum}p" |awk '{print $9}' `
    is_goodR=1 #`cat $evtdir/FK$ievt |sed -n "${lnum}p" |awk '{print $7}' `
    $prog/myconv $seisdir/$ievt/stf_pca001_Z.sac $seisdir/$ievt/$syn $seisdir/$ievt/$rec
    $prog/myconv $seisdir/$ievt/stf_pca001_R.sac $seisdir/$ievt/$synR $seisdir/$ievt/$recR
    #$prog/myconv $seisdir/$ievt/stf_pca001.sac $seisdir/$ievt/$synT $seisdir/$ievt/$recT
    tp=`cat $phasetab |sed -n "${ista}p" |awk '{printf"%f\n",$4}'` ##
    echo "is_good:" $syn $is_good $tp
    $prog/myadj $seisdir/$ievt/$dat $seisdir/$ievt/$syn $seisdir/$ievt/stf_pca001_Z.sac 0.03 0.2 $tp >>cost_valueZ.dat
    mv adjsource $seisdir/$ievt/$adj
    mv residual $seisdir/$ievt/$res
    $prog/myadj $seisdir/$ievt/$datR $seisdir/$ievt/$synR $seisdir/$ievt/stf_pca001_R.sac 0.03 0.2 $tp >>cost_valueR.dat
    mv adjsource $seisdir/$ievt/$adjR
    mv residual $seisdir/$ievt/$resR



    if [ $is_good -eq 0 ];then
      sac <<eof
r $seisdir/$ievt/$rec $seisdir/$ievt/$adj
mul 0
w over
quit
eof
    fi
    if [ $is_goodR -eq 0 ];then
      sac <<eof
r $seisdir/$ievt/$recR $seisdir/$ievt/$adjR
mul 0
w over
quit
eof
    fi
    let ista=ista+1
  done
  cat cost_valueZ.dat |awk 'BEGIN{sum=0;}{sum=sum+$2}END{printf"Total cost value: %f\n",sum}' >> cost_valueZ.dat
  cat cost_valueR.dat |awk 'BEGIN{sum=0;}{sum=sum+$2}END{printf"Total cost value: %f\n",sum}' >> cost_valueR.dat
  mv cost_value?.dat $seisdir/$ievt
  #====================================================  
  ### write STF to data dir 
  #~/progs/sactools_c/sac2col $seisdir/$ievt/stf_pca001_Z.sac |sed -n '2,$p' >$datadir/$ievt/STF 
  #~/progs/sactools_c/sac2col $seisdir/$ievt/stf_pca001_Z.sac |sed -n '2,$p' >>$datadir/$ievt/STF 
  #~/progs/sactools_c/sac2col $seisdir/$ievt/stf_pca001_Z.sac |sed -n '2,$p' >>$datadir/$ievt/STF 
  #is_stf=`grep source_time_function $datadir/$ievt/FK${ievt}`
  #if [ -z "$is_stf" ];then
  #  sed -i "13 a source_time_function '$datadir/$ievt/STF'" $datadir/$ievt/FK${ievt}
  #fi
  ###
  fi  # end of is_sac2asc
  #==================================================
  ### plot
  module load gcc hdf5/1.8.20 netcdf/4.6.1
gmtset ANNOT_FONT_SIZE_PRIMARY 12p
gmtset ANNOT_OFFSET_PRIMARY 0.1c
gmtset LABEL_FONT_SIZE 14p 
gmtset LABEL_OFFSET 0.15c
gmtset HEADER_FONT_SIZE 16p
gmtset HEADER_OFFSET -0.5c 
gmtset TICK_LENGTH -0.2c
sc=0.20

  for comp in BHR BHZ;do
    ls $seisdir/$ievt/dat.*${comp}*.sac.T005_T050 >file1.lst
    ls $seisdir/$ievt/syn.*${comp}*.sac.T005_T050 >file2.lst
    ls $seisdir/$ievt/dec.*BHZ*.sac.T005_T050 >file3.lst
    ls $seisdir/$ievt/PCs00?.sac        >file4.lst
    ls $seisdir/$ievt/PDMs*.sac    >file5.lst
    ls $seisdir/$ievt/rec.*${comp}*.sac.T005_T050 >file6.lst
    ls $seisdir/$ievt/adj.*${comp}*.sac.T005_T050 >file7.lst
    ls $seisdir/$ievt/res.*${comp}*.sac.T005_T050 >file8.lst
    nsta=`cat file2.lst |wc -l`
    len=40.
    npart=`echo $nsta $len |awk '{print $1/$2}' | awk '{printf("%d\n",$0+=$0<0?0:0.9)}'`
    echo "nsta,len,npart=" $nsta $len $npart
    for i in `seq $npart`;do
      ib=`echo $i $len |awk '{print ($1-1)*$2+1}'`
      ie=`echo $i $len |awk '{print $1*$2}'`
      if [ ${ie} -gt $nsta ];then
          ie=$nsta
      fi
      ymax=`echo $len |awk '{print $1+1}'`
      echo $ib $ie $ymax
    ##########################################
    out=$seisdir/fsismo_${ievt}_${comp}.seg${ib}-${ie}.ps
    psbasemap -R10/80/-1/$ymax -JX8/12 -Ba20f10/a1:"Trace number":weSn  -P -K -X3 -Y16 > $out
    cat file1.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p -O -K >>$out
    cat file2.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p,255/0/0 -O -K >>$out
    ty=`echo $ymax |awk '{print $1*1.05}'`
    pstext -J -R -N -K -O >>$seisdir/fsismo_${ievt}_${comp}.seg${ib}-${ie}.ps <<EOF
45 $ty 16 0 0 CM (a) DATA and SYN without STF 
EOF
    cat $seisdir/$ievt/FKtimes.dat |sed -n "${ib},${ie}p"  |awk '{printf"%f %f 8 0 0 CM %s\n",10,NR,$1}' | 
    pstext -J -R -W255/0/0 -N -K -O >>$out
    # (1) P
    cat $phasetab |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$4+t-b,NR-0.5,0,0.2}' t=$bt b=$twb | psxy -J -R -SV0.05/0./0. -Gmagenta  -N -O -K >>$out
    cat $phasetab |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$4+t+e,NR-0.5,0,0.2}' t=$bt e=$twe |
    psxy -J -R -SV0.05/0./0. -Gmagenta  -N -O -K >>$out
   
    ### plot RF function
    psbasemap -R-10/120/-1/$ymax -JX8/12 -Ba20f10/a1weSn -O -K -X9 >> $out
    cat file3.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p -O -K >> $out
    cat file5.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p,red -O -K >> $out
    ty=`echo $ymax |awk '{print $1*1.05}'`
    pstext -J -R -N -O -K >>$out <<EOF
55 $ty 16 0 0 CM (b) Deconvolved vertical traces 
EOF
    cat $seisdir/$ievt/FKtimes.dat |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s\n",-10,NR,$1}' | 
    pstext -J -R -W255/0/0 -N -K -O >>$out


    ### plot PC analysis
    psbasemap -R-10/120/0.0/2 -JX8/2 -Ba20f10:"Time (s)":/a10:"Amplitude":WeSn -X-9 -Y-13.5 -K -O >>$out
    pssac $seisdir/$ievt/stf_pca001.sac -J -R -Entb -M0.5 -W1p,red -O -K >>$out
    pstext -J -R -O -K >>$out <<EOF
105 0.25 12 0 0 CM STF
EOF

    xmax=`cat $seisdir/$ievt/contribution.dat |wc -l |awk '{print $1+1}'`
    tx=`echo $xmax |awk '{print $1/2}'`
    psbasemap -R0/$xmax/-5/100 -JX8/3 -Ba5/a20f10:,-%:WeSn -Y3 -K -O >>$out
    cat $seisdir/$ievt/contribution.dat |awk '{print NR,$1}' |
    psxy -J -R  -W1p -O -K >>$out
    cat $seisdir/$ievt/contribution.dat |awk '{print NR,$1}' |
    psxy -J -R -Sd0.1i -W0.5p -Gblue -O -K >>$out
    pstext -J -R -O -K >>$out <<EOF
$tx 80 12 0 0 CM PC contribution
EOF


    psbasemap -R-10/120/-1/10 -JX8/5 -Ba20f10/a1:"PCs number":WeSn  -O -K -Y4 >>$out
    cat file4.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p -O -K >> $out 
    head -1 file4.lst | pssac -J -R -Entb -M$sc -W1p,red -O -K >> $out
    pstext -J -R -O -K -N >>$out <<EOF
45 0 12 0 0 CM PC time series
55 11.3 16 0 0 CM (c) PCA analysis
EOF

    ### plot reconstruted waveform
    psbasemap -R10/80/0/$ymax -JX8/12 -Ba20f10:"Time (s)":/a1weSn  -O -K -X9 -Y-7 >>$out
    cat file1.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p -O -K >> $out 
    cat file6.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p,red -O -K >> $out 
    ty=`echo $ymax |awk '{print $1*1.05}'`
    pstext -J -R -N -O -K >>$out <<EOF
45 $ty 16 0 0 CM (d) DATA and SYN with STF 
EOF
    # (1) P
    cat $phasetab |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$4+t-b,NR-0.5,0,0.2}' t=$bt b=$twb | psxy -J -R -SV0.05/0./0. -Gmagenta  -N -O -K >>$out
    cat $phasetab |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$4+t+e,NR-0.5,0,0.2}' t=$bt e=$twe |
    psxy -J -R -SV0.05/0./0. -Gmagenta  -N -O -K >>$out
   

    cat $seisdir/$ievt/FKtimes.dat |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s\n",10,NR,$1}' | 
    pstext -J -R -W255/0/0 -N -O >>$out


    #######################################
    out1=$seisdir/adjsource_${ievt}_${comp}.seg${ib}-${ie}.ps

    psbasemap -R10/80/0/$ymax -JX8/12 -Ba20f10:"Time (s)":/a1:"Trace number":weSn  -P -K -X4 -Y4 > $out1
    cat file7.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p,255/0/0 -O -K >>$out1
    #cat $datadir/$ievt/CMTSOLUTION |sed -n '2p' |awk '{printf"50 34 16 0 1 CM %s: %s\n",$3,comp}' comp=$comp |
    #pstext -J -R -N -O -K >>$out1

    cat $seisdir/$ievt/FKtimes.dat |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s\n",10,NR,$1}' | 
    pstext -J -R -W255/0/0 -N -O -K >>$out1

    psbasemap -R10/80/0/$ymax -JX8/12 -Ba20f10:"Time (s)":/a1:"Trace number":wesn  -O -K  -Y12.5 >> $out1
    cat file8.lst |sed -n "${ib},${ie}p" | pssac -J -R -Entb -M$sc -W1p,255/0/0 -O -K >>$out1
    #cat $datadir/$ievt/CMTSOLUTION |sed -n '2p' |awk '{printf"%f 34 16 0 1 CM %s: %s\n",y+1,$3,comp}' y=$ymax comp=$comp |
    #pstext -J -R -N -O -K >>$out1

    cat $seisdir/$ievt/FKtimes.dat |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s\n",10,NR,$1}' | 
    pstext -J -R -W255/0/0 -N -O  >>$out1


    


    #rm file?.lst

    done
  done
done
