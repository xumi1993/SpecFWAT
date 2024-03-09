#!/bin/bash

indir=selectedSeismograms
outdir=data
datadir=DATA/DATA_SYNTHS/
sactool=~/progs/sactools_c
acqui=acqui_file.txt
stalst=stations_pick.lst

tmin=-10
tmax=109.99
zbottom=-300000
clat=43.25
clon=-1.25
decim=5   # the original HH? data has a delta of 0.01, we resample it to 0.05.
khole=00
netwk=X7

rm -rf $outdir
mkdir -p $outdir

### loop over events
cat $acqui |awk '{if(NR%2!=0) print $0}' |
while read line;do
  evtdir=`echo $line |awk '{print $3}'`
  evid=`echo $evtdir |awk -F$datadir '{print $2}'`
  evtnm=`echo $evtdir |awk -F$datadir '{printf"%s.S%03d\n",n,$2}' n=$netwk`
  evbaz=`grep BACK_AZIMUTH $evtdir/FKmodel |awk '{print $2}'`
  inc_ang=`grep TAKE_OFF $evtdir/FKmodel |awk '{print $2}'`
  echo $evtdir $evid
  mkdir -p $outdir/$evtnm
  ### calculate propogating time from the origin wavefront to receiver
  #python3 correct_wavefront_time.py $evtdir/FKmodel $zbottom $evtdir/FK$evid > time.out
  #tt=`tail -1 time.out |awk '{print -1*$4}'`
  #tt=`echo $tt $tmin |awk '{print $1-$2}'`
  python3 calculate_FK_time.py $evtdir/FKmodel $evtdir/FK$evid > time.out
  cp FKtimes.dat $datadir/$evid
  cp FKtimes.dat src_rec/FKtimes_${evtnm}
   ###
  # find the time difference between first motion and peak of P
  #t3=`saclst t3 f $indir/$evid/fstack.sac |awk '{print $2}'`
  t3=`cat $datadir/$evid/CMTSOLUTION |sed -n '3p' |awk '{print $3}'`
  echo $t3
  ### loop over station
  ista=1
  cat $evtdir/$stalst |
  while read line1;do
    stnm=`echo $line1 |awk '{print $1}'`
    netwk=`echo $line1 |awk '{print $2}'`
    stla=`echo $line1 |awk '{print $3}'`
    stlo=`echo $line1 |awk '{print $4}'`
    dist=`$sactool/distaz $stla $stlo $clat $clon |awk '{print $1*111.195}'`  ### dist: km
    stbaz=`$sactool/distaz $stla $stlo $clat $clon |awk '{print $2}'`
    staz=`$sactool/distaz $stla $stlo $clat $clon |awk '{print $3}'`
    
    #echo $stnm $dist $evbaz $stbaz
    bazdiff=`echo $evbaz $stbaz |awk '{print $2-$1}'`

    saci_r=$indir/$evid/$netwk.$stnm.$khole.HHR.sac
    saco_r=$outdir/$evtnm/$netwk.$stnm.BXR.sac
    saci_t=$indir/$evid/$netwk.$stnm.$khole.HHT.sac
    saco_t=$outdir/$evtnm/$netwk.$stnm.BXT.sac
    saci_z=$indir/$evid/$netwk.$stnm.$khole.HHZ.sac
    saco_z=$outdir/$evtnm/$netwk.$stnm.BXZ.sac

    ### Shift the seismograms by the propogating time above ###
    #tt=`cat FKtimes.dat |sed -n "${ista}p" |awk '{print -$4-a+b}' a=$tmin b=$t3` # we should shift the waveform by timeshift of GCMT
    tt=`cat FKtimes.dat |sed -n "${ista}p" |awk '{print -$4}' a=$tmin b=$t3`  # In Wang Yi's paper, they did not shift the waveform
    echo $stnm $t3 $tt

    ## for Z 
    stalistZ=`echo $stalst |awk -F. '{printf"%s_Z.lst",$1}'`
    is_data=`grep $stnm $datadir/$evid/$stalistZ` 
    if [ ! -z "$is_data" ];then
    $sactool/sac_shift $tt $saci_z $saco_z
    sac <<EOF
cut $tmin $tmax
r $saco_z
decimate $decim
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
decimate $decim
w over
quit
EOF
     fi

  let ista=ista+1
  done

done



