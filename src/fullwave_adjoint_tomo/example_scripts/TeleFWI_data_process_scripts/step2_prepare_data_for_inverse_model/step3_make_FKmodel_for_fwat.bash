#!/bin/bash
# Author  : Kai Wang
# Inst    : University of Toronto
# Email   : wangkaim8@gmail.com
# history : last modified on Jun 04, 2019

prog=~/progs/TauP-2.4.5/bin
sactool=~/progs/sactools_c
indir=../step1_download_data_sod
outdir=selectedSeismograms
datadir=DATA/DATA_SYNTHS
evtlst=final_events.lst
cmtdir=$indir/GCMT_solutions
stalst=stations_pick.lst
mod=ak135 ### prem or ak135
osdir=src_rec # operation system directory for the tomography package (ANAT)
phase=

#rm -rf $osdir
mkdir -p $osdir

cat /dev/null >$osdir/sources.dat
ievt=1 ### for TO
#ievt=44 ### for XJ
#ievt=59 ### for XE

IFS=$'\n'
#################### Main loop over events #####################
for evtdir in `cat $indir/final_events_P_TO2013-2015.lst`;do
#for evtdir in `cat $indir/final_events_P_XJ1997.lst`;do
#for evtdir in `cat $indir/final_events_P_XE2005-2007.lst`;do
  dirname=`echo $evtdir |awk '{print $1}'`
  evlo=`echo $evtdir |awk '{print $2}'`
  evla=`echo $evtdir |awk '{print $3}'`
  evdp=`echo $evtdir |awk '{print $4}'` # km
  evt=`echo $dirname |awk -FEvent_ '{print $2}'`
  evtnm=`echo $ievt |awk '{printf"S%03d",$1}'`
  echo $evt $evla $evlo $evdp
  #----------------------------------#
  #-----1. Obtain the original wavefront used in---------#
  #                   3DFK-SEM code                      #
  # Orignial wavefront in clat/clon/cdep. 3DFK-SEM code will use
  # the central of the mesh as xx0/yy0 in default. zz0 is set to
  # be the bottom z0 - 3*wavelength- ... (see the code for detail)
  #=========================================================================
  xmin=573223.498803403
  xmax=1042738.78269496
  ymin=3845599.93945242
  ymax=4251765.71778605
  xx0=`echo $xmin $xmax |awk '{printf"%.2f",($1+$2)/2}'`
  yy0=`echo $ymin $ymax |awk '{printf"%.2f",($1+$2)/2}'`
  zz0=-400000.0
  clon=-119.5
  clat=36.5                             # PARAMETERS TO BE SET BY USER
  cdep=400.0  ## depth: km
  tw=160.
  #=========================================================================
  dist=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $1}'`
  baz=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $2}'`
  az=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $3}'`
  $prog/taup_time  -mod prem -h $evdp  -ph P -deg $dist  >taup.out
  #$prog/taup_time  -mod prem -h $evdp  -ph P -deg $dist --stadepth $cdep >taup.out
  inc_ang=`cat taup.out |sed -n '6p' |awk '{print $7}'`  ## Incident angle
  #-----------------------------------------#
 cat >$osdir/FKmodel_${phase}${ievt} <<EOF
#  input file for embedded FK modeiling
#
#  for each layer we give :
#  LAYER ilayer rho vp vs ztop
#  the last layer is the homogeneous half space
#
#
# model description  ---------------------
EOF
  #******* Use PREM layers ***************
  nlay=`python make_${mod}_fkmodel.py |wc -l`
  echo "NLAYER             $nlay"  >>./${osdir}/FKmodel_${phase}${ievt}
  python make_${mod}_fkmodel.py >>./$osdir/FKmodel_${phase}${ievt}
  #******* Use simple crust/mantle model from Ping et al., GRL, 2014
#  cat >>./$datadir/$ievt/FKmodel <<EOF
#NLAYER             2
#LAYER 1 2600.000 5800.000 3198.000   0000.000
#LAYER 2 3380.000 8080.000 4485.000 -30000.000
#EOF
  #****************************************
  cat >>./$osdir/FKmodel_${phase}${ievt} <<EOF
#----------------------------------------
# incident wave p or sv
INCIDENT_WAVE      p
#----------------------------------------
# anlges of incomming wave
BACK_AZIMUTH               $baz
TAKE_OFF               $inc_ang
#----------------------------------------
FREQUENCY_MAX      4
#----------------------------------------
TIME_WINDOW       $tw
#----------------------------------------
# optionnal
 ORIGIN_WAVEFRONT  $xx0       $yy0      $zz0
# ORIGIN_TIME 
EOF

  #--------------------------------------------------#
  #------  2. event station list ----
  # STATIONS
  #for comp in R Z;do
    #mylist=`echo $stalst |awk -F. '{printf"%s_%s.lst",$1,v}' v=$comp`
    #cp $datadir/$ievt/$mylist $osdir/STATIONS_${ievt}_$comp 
    cp $datadir/$ievt/$stalst $osdir/STATIONS_${phase}${ievt} 
  #done
  # CMTSOLUTION: Although the forward simulation do not use CMTSOLUTION
  # but we need the source lacation in CMTSOLUTION for rotation in the 
  # data processing of the tomography package
  cat >$osdir/CMTSOLUTION_${phase}${ievt} <<eof
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
  echo $phase$ievt $evla $evlo 0.0 0.0 >>$osdir/sources.dat 
  #----------------------------------------------------
  let ievt=ievt+1
done
########## END loop over events ###############################




