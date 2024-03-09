#!/bin/bash

prog=~/progs/TauP-2.4.5/bin
sactool=~/progs/sac_msc
indir=../step1_download_data_sod
outdir=selectedSeismograms
datadir=DATA/DATA_SYNTHS
evtlst=final_events.lst
cmtdir=$indir/GCMT_solutions
stlst=station_TO.lst


ievt=1
IFS=$'\n'
#################### Main loop over events #####################
for evtdir in `cat $indir/final_events.lst`;do
  dirname=`echo $evtdir |awk '{print $1}'`
  evlo=`echo $evtdir |awk '{print $2}'`
  evla=`echo $evtdir |awk '{print $3}'`
  evdp=`echo $evtdir |awk '{print $4}'` # km
  evt=`echo $dirname |awk -FEvent_ '{print $2}'`
  echo $evt $evla $evlo $evdp
  #----------------------------------#
  #-----1. Obtain the original wavefront used in---------#
  #                   3DFK-SEM code                      #
  # Orignial wavefront in lat/lon/depth. 3DFK-SEM code will use
  # the central of the mesh as xx0/yy0 in default. zz0 is set to
  # be the bottom z0 - 3*wavelength- ... (see the code for detail)
  #=========================================================================
  xx0=-119.95
  yy0=35.975                             # PARAMETERS TO BE SET BY USER
  zz0=300.0  ## depth: km
  tw=160.
  #=========================================================================
  dist=`$sactool/distaz $yy0 $xx0 $evla $evlo |awk '{print $1}'`
  baz=`$sactool/distaz $yy0 $xx0 $evla $evlo |awk '{print $2}'`
  az=`$sactool/distaz $yy0 $xx0 $evla $evlo |awk '{print $3}'`
  $prog/taup_time  -mod prem -h $evdp  -ph P -deg $dist  >taup.out
  #$prog/taup_time  -mod prem -h $evdp  -ph P -deg $dist --stadepth $zz0 >taup.out
  inc_ang=`cat taup.out |sed -n '6p' |awk '{print $7}'`  ## Incident angle
  #-----------------------------------------#
 cat >./$datadir/$ievt/FKmodel <<EOF
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
  nlay=`python make_prem_fkmodel.py |wc -l`
  echo "NLAYER             $nlay"  >>./$datadir/$ievt/FKmodel
  python make_prem_fkmodel.py >>./$datadir/$ievt/FKmodel
  #******* Use simple crust/mantle model from Ping et al., GRL, 2014
#  cat >>./$datadir/$ievt/FKmodel <<EOF
#NLAYER             2
#LAYER 1 2600.000 5800.000 3198.000   0000.000
#LAYER 2 3380.000 8080.000 4485.000 -30000.000
#EOF
  #****************************************
  cat >>./$datadir/$ievt/FKmodel <<EOF
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
 ORIGIN_WAVEFRONT  768145.1       3990028.      -300000.0
# ORIGIN_TIME $t0
EOF

  nsta=`cat $datadir/$ievt/stations_pick.lst |wc -l`
  #========================================================
  # set it to the same as used in SEM
  dt=0.05                    # PARAMETERS TO BE SET BY USER
  nt=2400
  #========================================================
  cat >$datadir/$ievt/FK$ievt <<EOF
event_name FK$ievt
source_type moment
source_components './DATA/DATA_SYNTHS/$ievt/CMTSOLUTION'
modeling_tool 'fk' './DATA/DATA_SYNTHS/$ievt/FKmodel'
data_components d rtz
cartloc_mesh_origin 0. 0. 0.
data_origin_time 0
number_of_station $nsta
data_time_step $dt
data_sample_number $nt
is_time_pick false
time_window false 0 0
station_coord_system geo
###################################################
STARTLIST
EOF

#source_time_function './DATA/DATA_SYNTHS/$ievt/STF'

  cat $datadir/$ievt/stations_pick.lst >>$datadir/$ievt/FK$ievt
  echo ENDLIST >>$datadir/$ievt/FK$ievt
  #---------------------------------------------------------------#
  #-----------------------2. get Source Time Function ---------------#
  #python write_ascii_source_signature.py $outdir/$ievt/fstack.sac.tp.stf $datadir/$ievt/STF
  #-----------------------3 .Set CMTSOLUTION -------------------------#
  evtnm=`echo $dirname |awk '{print $1}' |awk -F/ '{print $3}'`
  cp $cmtdir/CMTSOLUTION_${evtnm} $datadir/$ievt/CMTSOLUTION
  #----------------------------------------------------#
  let ievt=ievt+1
done
########## END loop over events ###############################
