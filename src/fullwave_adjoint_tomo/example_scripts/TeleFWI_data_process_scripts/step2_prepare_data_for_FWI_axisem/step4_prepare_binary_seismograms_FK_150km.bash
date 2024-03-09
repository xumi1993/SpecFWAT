#!/bin/bash
module load anaconda2
module load java
module load gcc


indir=selectedSeismograms
outdir=preparedSeismograms
datadir=DATA/DATA_SYNTHS/
sactool=~/progs/sactools_c
acqui=acqui_file.txt
stalst=stations_pick.lst

tmin=-10.
tmax=109.975
zbottom=-300000
clat=43.25
clon=-1.25
dt=0.025
nt=4800

rm -rf $outdir
mkdir -p $outdir

ievt=1
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
  python calculate_FK_time.py $evtdir/FKmodel $evtdir/FK$evid > time.out
  #python calculate_FK_time_with_topo_sph_correct.py $evtdir/FKmodel $evtdir/FK$evid $datadir/$evid/time.pred > time.out
  cp FKtimes.dat $datadir/$evid
  cd $datadir/$evid
  cat /dev/null >FKfile
  head -15 FK$evid >>FKfile
  cat FKtimes.dat |
  while read line;do
    stnm=`echo $line |awk '{print $2}'`
    Ptime=`echo $line |awk '{print $4+t}' t=$tmin`
    grep $stnm FK$evid |awk '{print $1,$2,$3,$4,$5,tp,$7,$8,$9}' tp=$Ptime >>FKfile
  done
  tail -1 FK$evid >>FKfile
  mv FKfile FK$evid
  cd -
  ###
  # find the time difference between first motion and peak of P
  #t3=`saclst t3 f $indir/$evid/fstack.sac |awk '{print $2}'`
  #t3=`cat $datadir/$evid/CMTSOLUTION |sed -n '3p' |awk '{print $3}'`
  #echo $t3
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
    ### Assume that
    #STC=`echo $dist $inc_ang $$bazdiff |awk '{print $1*sin($2/180*3.14159265)}'`

    khole=00
    saci_r=$indir/$evid/$netwk.$stnm.$khole.BHR.sac
    saco_r=$outdir/$evid/$netwk.$stnm.$khole.BHR.sac
    saci_t=$indir/$evid/$netwk.$stnm.$khole.BHT.sac
    saco_t=$outdir/$evid/$netwk.$stnm.$khole.BHT.sac
    saci_z=$indir/$evid/$netwk.$stnm.$khole.BHZ.sac
    saco_z=$outdir/$evid/$netwk.$stnm.$khole.BHZ.sac

    ### Shift the seismograms by the propogating time above ###

    #tt=`cat FKtimes.dat |sed -n "${ista}p" |awk '{print -$4-a+b}' a=$tmin b=$t3` # we should shift the waveform by timeshift of GCMT
    tt=`cat FKtimes.dat |sed -n "${ista}p" |awk '{print -$4-a}' a=$tmin`  # In Wang Yi's paper, they did not shift the waveform
    echo $stnm $tt
    $sactool/sac_shift $tt $saci_r $saco_r
    $sactool/sac_shift $tt $saci_t $saco_t
    $sactool/sac_shift $tt $saci_z $saco_z
    sac <<EOF
r $saco_z $saco_r $saco_t
interp delta $dt begin $tmin
w over
cut $tmin $tmax
r $saco_z $saco_r $saco_t
mul -1.
w over
quit
EOF

#    sac <<EOF
#cut $tmin $tmax
#r $saco_z $saco_r $saco_t
#mul -1.
#w over
#quit
#EOF

  let ista=ista+1
  done
  cat $evtdir/$stalst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s.00.BHR.sac\n",dir,$2,$1}' >$outdir/stations_evt$evid.R.lst
  cat $evtdir/$stalst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s.00.BHT.sac\n",dir,$2,$1}' >$outdir/stations_evt$evid.T.lst
  cat $evtdir/$stalst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s.00.BHZ.sac\n",dir,$2,$1}' >$outdir/stations_evt$evid.Z.lst
  nrec=`cat $outdir/stations_evt$evid.R.lst |wc -l  `
  python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.R.lst $nrec $dt $nt BHR  $datadir/$evid/fsismo_dr.bin
  nrec=`cat $outdir/stations_evt$evid.T.lst |wc -l  `
  python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.T.lst $nrec $dt $nt BHT  $datadir/$evid/fsismo_dt.bin
  nrec=`cat $outdir/stations_evt$evid.Z.lst |wc -l  `
  python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.Z.lst $nrec $dt $nt BHZ  $datadir/$evid/fsismo_dz.bin
  let ievt=ievt+1
done



