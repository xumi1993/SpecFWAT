#!/bin/bash

indir=selectedSeismograms
outdir=preparedSeismograms
datadir=./DATA/DATA_SYNTHS/
sactool=~/progs/sactools_c
acqui=acqui_file.txt

tmin=-10
tmax=109.99
zbottom=-300000
clat=35.975
clon=-119.95
decim=5   # the original HH? data has a delta of 0.01, we resample it to 0.05.
dt=0.05
nt=2400

rm -rf $outdir
mkdir -p $outdir

### loop over events
cat $acqui |sed -n '3,4p' |awk '{if(NR%2!=0) print $0}' |
while read line;do
  evtdir=`echo $line |awk '{print $3}'`
  evid=`echo $evtdir |awk -F$datadir '{print $2}'`
  evbaz=`grep BACK_AZIMUTH $evtdir/FKmodel |awk '{print $2}'`
  inc_ang=`grep TAKE_OFF $evtdir/FKmodel |awk '{print $2}'`
  echo $evtdir $evid
  mkdir -p $outdir/$evid
  ### calculate propogating time from the origin wavefront to receiver
  #python3 correct_wavefront_time.py $evtdir/FKmodel $zbottom $evtdir/FK$evid > time.out
  #tt=`tail -1 time.out |awk '{print -1*$4}'`
  #tt=`echo $tt $tmin |awk '{print $1-$2}'`
  #**************remember to change the UTM zone in calculate_FK_time.py *************
  python calculate_FK_time.py $evtdir/FKmodel $evtdir/FK$evid > time.out
  cp FKtimes.dat $datadir/$evid
  ###
  # find the time difference between first motion and peak of P
  #t3=`cat $datadir/$evid/CMTSOLUTION |sed -n '3p' |awk '{print $3}'`
  t3=`grep "Mean_arrival_time:" $indir/$evid/HHZ.mcp |awk '{print $2}' `
  echo $t3
  ### loop over station
  ista=1
  cat $evtdir/stations_pick.lst |
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
    ### Assume that
    #STC=`echo $dist $inc_ang $$bazdiff |awk '{print $1*sin($2/180*3.14159265)}'`

    saci_r=$indir/$evid/$netwk.$stnm..HHR.sac
    saco_r=$outdir/$evid/$netwk.$stnm..HHR.sac
    saci_t=$indir/$evid/$netwk.$stnm..HHT.sac
    saco_t=$outdir/$evid/$netwk.$stnm..HHT.sac
    saci_z=$indir/$evid/$netwk.$stnm..HHZ.sac
    saco_z=$outdir/$evid/$netwk.$stnm..HHZ.sac

    ### Shift the seismograms by the propogating time above ###

    tt=`cat FKtimes.dat |sed -n "${ista}p" |awk '{print -$4-a+b}' a=$tmin b=$t3`
    echo $stnm $t3 $tt
    $sactool/sac_shift $tt $saci_r $saco_r
    $sactool/sac_shift $tt $saci_t $saco_t
    $sactool/sac_shift $tt $saci_z $saco_z
    sac <<EOF
cut $tmin $tmax
r $saco_z $saco_r $saco_t
decimate $decim
w over
quit
EOF
  let ista=ista+1
  done
  cat $evtdir/stations_pick.lst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s..HHR.sac\n",dir,$2,$1}' >$outdir/stations_evt$evid.HHR.lst
  cat $evtdir/stations_pick.lst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s..HHT.sac\n",dir,$2,$1}' >$outdir/stations_evt$evid.HHT.lst
  cat $evtdir/stations_pick.lst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s..HHZ.sac\n",dir,$2,$1}' >$outdir/stations_evt$evid.HHZ.lst
  nrec=`cat $outdir/stations_evt$evid.HHZ.lst |wc -l  `
  python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.HHR.lst $nrec $dt $nt HHR  $datadir/$evid/fsismo_dr.bin
  python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.HHT.lst $nrec $dt $nt HHT  $datadir/$evid/fsismo_dt.bin
  python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.HHZ.lst $nrec $dt $nt HHZ  $datadir/$evid/fsismo_dz.bin
done



