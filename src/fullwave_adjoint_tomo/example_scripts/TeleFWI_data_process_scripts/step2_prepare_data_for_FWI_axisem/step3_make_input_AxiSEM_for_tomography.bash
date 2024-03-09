#!/bin/bash
# Author  : Kai Wang
# Inst    : University of Toronto
# Email   : wangkaim8@gmail.com
# history : last modified on Oct 09, 2019

prog=~/progs/TauP-2.4.5/bin
sactool=~/progs/sactools_c
indir=../step1_download_data_sod/Picked_8events_20190801
outdir=selectedSeismograms
datadir=DATA/DATA_SYNTHS
evtlst=final_events.lst
cmtdir=$indir/GCMT_solutions
stalst=stations_pick.lst
osdir=src_rec # operation system directory for the tomography package (ANAT)
netwk=TO

rm -rf $osdir
mkdir -p $osdir

cat /dev/null >$osdir/sources.dat
ievt=1
IFS=$'\n'
#################### Main loop over events #####################
for evtdir in `cat $indir/final_events.lst`;do
  dirname=`echo $evtdir |awk '{print $1}'`
  evlo=`echo $evtdir |awk '{print $2}'`
  evla=`echo $evtdir |awk '{print $3}'`
  evdp=`echo $evtdir |awk '{print $4}'` # km
  evt=`echo $dirname |awk -FEvent_ '{print $2}'`
  evtnm=`echo $ievt |awk '{printf"S%03d",$1}'`
  evtdir=`echo $netwk $ievt |awk '{printf"%s.S%03d",$1,$2}'`
  echo $evt $evla $evlo $evdp
  #----------------------------------#
  #--------------------------------------------------#
  #------  2. event station list ----
  # STATIONS
  cp $datadir/$ievt/$stalst $osdir/STATIONS_${evtdir}
  for comp in R Z;do
    mylist=`echo $stalst |awk -F. '{printf"%s_%s.lst",$1,v}' v=$comp`
    cp $datadir/$ievt/$mylist $osdir/STATIONS_${evtdir}_$comp 
  done
  # CMTSOLUTION: Although the forward simulation do not use CMTSOLUTION
  # but we need the source lacation in CMTSOLUTION for rotation in the 
  # data processing of the tomography package
  cat >$osdir/CMTSOLUTION_${evtdir} <<eof
PDE  1999 01 01 00 00 00.00  $evlo $evla -25000 4.2 4.2 hom_explosion
event name:       hom_explosion
time shift:       0.0000
half duration:    0.0
latorUTM:         $evla
longorUTM:        $evlo
depth:            $evdp
Mrr:       1.000000e+23
Mtt:       1.000000e+23
Mpp:       1.000000e+23
Mrt:       0.000000
Mrp:       0.000000
Mtp:       0.000000
eof
  # sources.dat
  echo $evtnm $netwk $evla $evlo 0.0 0.0 >>$osdir/sources.dat 
  #----------------------------------------------------
  let ievt=ievt+1
done
########## END loop over events ###############################




