#!/bin/bash
#module load anaconda2
#module load java
#module load gcc


indir=selectedSeismograms
outdir=preparedSeismograms
datadir=DATA/DATA_SYNTHS/
sactool=~/progs/sactools_c
acqui=acqui_file.txt
stalst=stations_pick.lst

tmin=0.
tmax=119.95
clat=35.975
clon=-119.95
dt=0.05
nt=2400

rm -rf $outdir
mkdir -p $outdir

ievt=1
### loop over events
cat $acqui |awk '{if(NR%2!=0) print $0}' |
while read line;do
  evtdir=`echo $line |awk '{print $3}'`
  evid=`echo $evtdir |awk -F$datadir '{print $2}'`
  echo $evtdir $evid
  mkdir -p $outdir/$evid
  ###############################################
  ### calculate propogating t
  if [ $ievt -eq 1 ];then
     ref=645.
     tshift=-37.2085  ### time difference between reference and firt motion picked by AIMBAT
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  elif [ $ievt -eq 2 ];then
     ref=715.
     tshift=-2.5226
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  elif [ $ievt -eq 3 ];then
     ref=460.
     tshift=0.2979
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  elif [ $ievt -eq 4 ];then
     ref=625.
     tshift=-3.3693
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  elif [ $ievt -eq 5 ];then
     ref=445.
     tshift=-3.1886
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
   elif [ $ievt -eq 6 ];then
     ref=720.
     tshift=-0.9113
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
   elif [ $ievt -eq 7 ];then
     ref=650.
     tshift=-23.148
     cat $datadir/$evid/time.pred |awk '{printf"%s %s %f %f\n",$2,$1,$6-t,$6-t+t1}' t=$ref t1=$tshift >PIFtimes.dat
  else
     ref=630.
     tshift=-0.1946
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
  ###################################################
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

    khole=''
    saci_r=$indir/$evid/$netwk.$stnm.$khole.BHR.sac
    saco_r=$outdir/$evid/$netwk.$stnm.$khole.BHR.sac
    saci_t=$indir/$evid/$netwk.$stnm.$khole.BHT.sac
    saco_t=$outdir/$evid/$netwk.$stnm.$khole.BHT.sac
    saci_z=$indir/$evid/$netwk.$stnm.$khole.BHZ.sac
    saco_z=$outdir/$evid/$netwk.$stnm.$khole.BHZ.sac

    ### Shift the seismograms by the propogating time above ###

    tt=`cat $datadir/$evid/PIFtimes.dat |sed -n "${ista}p" |awk '{print -$4-a}' a=$tmin b=$t3` # we should shift the waveform by timeshift of GCMT
    #tt=`cat PIFtimes.dat |sed -n "${ista}p" |awk '{print -$3}'`  # In Wang Yi's paper, they did not shift the waveform
    echo $stnm $tt
    $sactool/sac_shift $tt $saci_r $saco_r
    $sactool/sac_shift $tt $saci_t $saco_t
    $sactool/sac_shift $tt $saci_z $saco_z
#    sac <<EOF
#r $saco_z $saco_r $saco_t
#interp delta $dt begin $tmin
#w over
#cut $tmin $tmax
#r $saco_z $saco_r $saco_t
#mul -1.
#w over
#quit
#EOF

    sac <<EOF
cuterr fillz
cut $tmin $tmax
r $saco_z $saco_r $saco_t
w over
quit
EOF

  let ista=ista+1
  done
  cat $evtdir/$stalst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s.%s.BHR.sac\n",dir,$2,$1,k}' k=$khole >$outdir/stations_evt$evid.R.lst
  cat $evtdir/$stalst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s.%s.BHT.sac\n",dir,$2,$1,k}' k=$khole >$outdir/stations_evt$evid.T.lst
  cat $evtdir/$stalst |awk -v dir=$outdir/$evid '{printf"%s/%s.%s.%s.BHZ.sac\n",dir,$2,$1,k}' k=$khole >$outdir/stations_evt$evid.Z.lst
  nrec=`cat $outdir/stations_evt$evid.R.lst |wc -l  `
  python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.R.lst $nrec $dt $nt BHR  $datadir/$evid/fsismo_dr.bin
  nrec=`cat $outdir/stations_evt$evid.T.lst |wc -l  `
  python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.T.lst $nrec $dt $nt BHT  $datadir/$evid/fsismo_dt.bin
  nrec=`cat $outdir/stations_evt$evid.Z.lst |wc -l  `
  python ./convert_sac_to_gather_binary.py $outdir/stations_evt$evid.Z.lst $nrec $dt $nt BHZ  $datadir/$evid/fsismo_dz.bin
  let ievt=ievt+1
done



