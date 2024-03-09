#!/bin/bash
# Author  : Kai Wang
# Inst    : University of Toronto
# Email   : wangkaim8@gmail.com
# history : last modified on Oct 09, 2019


indir=selectedSeismograms
outdir=data
osdir=src_rec # operation system directory for the tomography package (ANAT)
datadir=DATA/DATA_SYNTHS/
sactool=~/progs/sactools_c
acqui=acqui_file.txt
stalst=stations_pick.lst

tmin=-0.5
tmax=119.45
clat=35.975
clon=-119.95
khole=''
netwk=TO

rm -rf $outdir
mkdir -p $outdir

ievt=1
### loop over events
cat $acqui |awk '{if(NR%2!=0) print $0}' |
while read line;do
  evtdir=`echo $line |awk '{print $3}'`
  evid=`echo $evtdir |awk -F$datadir '{print $2}'`
  evtnm=`echo $evtdir |awk -F$datadir '{printf"%s.S%03d\n",n,$2}' n=$netwk`
  echo $evtdir $evid
  mkdir -p $outdir/$evtnm
  #######################################################
  ### calculate propogating t
  if [ $ievt -eq 1 ];then
     ref=645.
     tshift=-30.
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  elif [ $ievt -eq 2 ];then
     ref=715.
     tshift=0.
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  elif [ $ievt -eq 3 ];then
     ref=460.
     tshift=0.
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  elif [ $ievt -eq 4 ];then
     ref=625.
     tshift=0.
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  elif [ $ievt -eq 5 ];then
     ref=445.
     tshift=0.
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
   elif [ $ievt -eq 6 ];then
     ref=720.
     tshift=10.
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
   elif [ $ievt -eq 7 ];then
     ref=650.
     tshift=-15.
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  else
     ref=630.
     tshift=10.
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  fi
  
  ###############################################
  cd $datadir/$evid
  cat /dev/null >PIFfile
  cat /dev/null >PIFtimes.dat
  head -15 PIF$evid >>PIFfile
  cat ../../../PIFtimes.dat |
  while read line;do
    netwk=`echo $line |awk '{print $1}'`
    stnm=`echo $line |awk '{print $2}'`
    dt1=`echo $line |awk '{print $3}'`
    dt2=`echo $line |awk '{print $4}'`
    Ptime=`echo $line |awk '{print $3+t}' t=$tmin`
    grep $stnm PIF$evid |awk '{print $1,$2,$3,$4,$5,tp,$7,$8,$9}' tp=$Ptime >>PIFfile
    grep $stnm PIF$evid |awk '{printf"%s %s %f %f\n",net,st,t1,t2}' net=$netwk st=$stnm t1=$dt1 t2=$dt2 >>PIFtimes.dat
  done
  tail -1 PIF$evid >>PIFfile
  mv PIFfile PIF$evid
  cd -
  rm PIFtimes.dat 
  cp $datadir/$evid/PIFtimes.dat $osdir/PIFtimes_${evtnm}.dat 
  ###################################################

  ######################################################################
  ### loop over station
  ista=1
  cat $evtdir/$stalst |
  while read line1;do
    stnm=`echo $line1 |awk '{print $1}'`
    netwk=`echo $line1 |awk '{print $2}'`
    #stla=`echo $line1 |awk '{print $3}'`
    #stlo=`echo $line1 |awk '{print $4}'`
    #dist=`$sactool/distaz $stla $stlo $clat $clon |awk '{print $1*111.195}'`  ### dist: km
    #stbaz=`$sactool/distaz $stla $stlo $clat $clon |awk '{print $2}'`
    #staz=`$sactool/distaz $stla $stlo $clat $clon |awk '{print $3}'`
    

    saci_r=$indir/$evid/$netwk.$stnm.$khole.BHR.sac
    saco_r=$outdir/$evtnm/$netwk.$stnm.BXR.sac
    saci_t=$indir/$evid/$netwk.$stnm.$khole.BHT.sac
    saco_t=$outdir/$evtnm/$netwk.$stnm.BXT.sac
    saci_z=$indir/$evid/$netwk.$stnm.$khole.BHZ.sac
    saco_z=$outdir/$evtnm/$netwk.$stnm.BXZ.sac

    ### Shift the seismograms by the propogating time above ###
    tt=`cat $datadir/$evid/PIFtimes.dat |sed -n "${ista}p" |awk '{print -$4-a}' a=$tmin b=$t3` # we should shift the waveform by timeshift of GCMT
    #tt=`cat PIFtimes.dat |sed -n "${ista}p" |awk '{print -$3}'`  # In Wang Yi's paper, they did not shift the waveform
    echo $stnm $tt
    $sactool/sac_shift $tt $saci_r $saco_r
    $sactool/sac_shift $tt $saci_t $saco_t
    $sactool/sac_shift $tt $saci_z $saco_z

    ## for Z 
    stalistZ=`echo $stalst |awk -F. '{printf"%s_Z.lst",$1}'`
    is_data=`grep $stnm $datadir/$evid/$stalistZ` 
    if [ ! -z "$is_data" ];then
    $sactool/sac_shift $tt $saci_z $saco_z
    sac <<EOF
cut $tmin $tmax
r $saco_z
w over
quit
EOF
     fi
    ## for R
    stalistR=`echo $stalst |awk -F. '{printf"%s_R.lst",$1}'`
    is_data=`grep $stnm $datadir/$evid/$stalistR` 
    if [ ! -z "$is_data" ];then
    $sactool/sac_shift $tt $saci_r $saco_r
    sac <<EOF
cut $tmin $tmax
r $saco_r
w over
quit
EOF
     fi

  let ista=ista+1
  done
  ######################################################################
  let ievt=ievt+1
done



