#!/bin/bash

indir=selectedSeismograms
outdir=preparedSeismograms
datadir=./DATA/DATA_SYNTHS/
sactool=~/progs/sactools_c
acqui=acqui_file.txt

tmin=-10
tmax=109.975
clat=36.5
clon=-119.5
decim=5   # the original BH? data has a delta of 0.01, we resample it to 0.05.
dt=0.025
nt=4800

rm -rf $outdir
mkdir -p $outdir

### loop over events
cat $acqui |awk '{if(NR%2!=0) print $0}' |
while read line;do
  evtdir=`echo $line |awk '{print $3}'`
  evid=`echo $evtdir |awk -F$datadir '{print $2}'`
  evbaz=`grep BACK_AZIMUTH $evtdir/FKmodel |awk '{print $2}'`
  inc_ang=`grep TAKE_OFF $evtdir/FKmodel |awk '{print $2}'`
  echo $evtdir $evid
  mkdir -p $outdir/$evid
  ### calculate propogating time from the origin wavefront to receiver
  #python3 correct_wavefront_time.py $evtdir/FKmodel $zbottom $evtdir/PIF$evid > time.out
  #tt=`tail -1 time.out |awk '{print -1*$4}'`
  #tt=`echo $tt $tmin |awk '{print $1-$2}'`
  #**************remember to change the UTM zone in calculate_FK_time.py *************
  ###
  # find the time difference between first motion and peak of P
  #tm=`cat $datadir/$evid/CMTSOLUTION |sed -n '3p' |awk '{print $3}'`
  tm=`grep "Mean_arrival_time:" $indir/$evid/*BHZ.mcp |awk '{print $2}' `
  tw=`grep "Window:" $indir/$evid/*BHZ.mcp |awk '{print $2}' `
  tb=2.0
  te=`echo $tw |awk '{printf"%.1f\n",$1-2.}' `
  echo $tm $tw $tb $te
  gsed -i "/time_window/c\time_window true $tb $te" ${evtdir}/PIF$evid 
  python calculate_FK_time.py $evtdir/FKmodel $evtdir/PIF$evid > time.out
  cp FKtimes.dat $datadir/$evid
  cat FKtimes.dat |awk '{print $1,$2,$4,a,b}' a=$tb b=$tw >src_rec/FKtimes_$evid
  cd $datadir/$evid
  cat /dev/null >FKfile
  head -15 PIF$evid >>FKfile
  cat FKtimes.dat |
  while read line;do
    stnm=`echo $line |awk '{print $2}'`
    #Ptime=`echo $line |awk '{print $4+t}' t=$tmin`
    Ptime=`echo $line |awk '{print $4}' `
    grep $stnm PIF$evid |awk '{print $1,$2,$3,$4,$5,tp,$7,$8,$9}' tp=$Ptime >>FKfile
  done
  tail -1 PIF$evid >>FKfile
  mv FKfile PIF$evid
  cd -
  ### loop over station
  if true;then
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

    saci_r=$indir/$evid/$netwk.$stnm..BHR.sac
    saco_r=$outdir/$evid/$netwk.$stnm.BXR.sac
    saci_t=$indir/$evid/$netwk.$stnm..BHT.sac
    saco_t=$outdir/$evid/$netwk.$stnm.BXT.sac
    saci_z=$indir/$evid/$netwk.$stnm..BHZ.sac
    saco_z=$outdir/$evid/$netwk.$stnm.BXZ.sac

    ### Shift the seismograms by the propogating time above ###

    tt=`cat FKtimes.dat |sed -n "${ista}p" |awk '{print -$4+$5-a+b}' a=$tmin b=$tm`
    ttp=`cat FKtimes.dat |sed -n "${ista}p" |awk '{print $4}'`
    echo $stnm $tm $tt
    $sactool/sac_shift $tt $saci_r $saco_r
    $sactool/sac_shift $tt $saci_t $saco_t
    $sactool/sac_shift $tt $saci_z $saco_z
    sac <<EOF
cut $tmin $tmax
r $saco_z $saco_r $saco_t
ch t0 $ttp
w over
quit
EOF
  let ista=ista+1
  done
  #cat $evtdir/stations_pick.lst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s..BHR.sac\n",dir,$2,$1}' >$outdir/stations_evt$evid.BHR.lst
  #cat $evtdir/stations_pick.lst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s..BHT.sac\n",dir,$2,$1}' >$outdir/stations_evt$evid.BHT.lst
  #cat $evtdir/stations_pick.lst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s..BHZ.sac\n",dir,$2,$1}' >$outdir/stations_evt$evid.BHZ.lst
  #nrec=`cat $outdir/stations_evt$evid.BHZ.lst |wc -l  `
  #python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.BHR.lst $nrec $dt $nt BHR  $datadir/$evid/fsismo_dr.bin
  #python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.BHT.lst $nrec $dt $nt BHT  $datadir/$evid/fsismo_dt.bin
  #python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.BHZ.lst $nrec $dt $nt BHZ  $datadir/$evid/fsismo_dz.bin
  fi
done



