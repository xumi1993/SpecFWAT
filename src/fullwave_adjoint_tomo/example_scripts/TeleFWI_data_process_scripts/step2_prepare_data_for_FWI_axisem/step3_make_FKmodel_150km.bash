#!/bin/bash
# Author  : Kai Wang
# Inst    : University of Toronto
# Email   : wangkaim8@gmail.com
# history : last modified on Jun 04, 2019

module load anaconda2
module load java



#module load anaconda3
prog=~/progs/TauP-2.4.5/bin
sactool=~/progs/sactools_c
indir=../step1_download_data_sod
outdir=selectedSeismograms
datadir=DATA/DATA_SYNTHS
evtlst=final_events.lst
cmtdir=$indir/GCMT_solutions
stalst=stations_pick.lst
mod=ak135 ### prem or ak135



ievt=1
IFS=$'\n'
#################### Main loop over events #####################
for evtdir in `cat $indir/final_events_wang2016.lst`;do
  dirname=`echo $evtdir |awk '{print $1}'`
  evlo=`echo $evtdir |awk '{print $2}'`
  evla=`echo $evtdir |awk '{print $3}'`
  evdp=`echo $evtdir |awk '{print $4}'` # km
  evt=`echo $dirname |awk -FEvent_ '{print $2}'`
  echo $evt $evla $evlo $evdp
  #----------------------------------#
  #-----1. Obtain the original wavefront used in---------#
  #                   3DFK-SEM code                      #
  # Orignial wavefront in clat/clon/cdep. 3DFK-SEM code will use
  # the central of the mesh as xx0/yy0 in default. zz0 is set to
  # be the bottom z0 - 3*wavelength- ... (see the code for detail)
  #=========================================================================
  xx0=640778.881
  yy0=4791747.05
  zz0=-300000.0
  clon=-1.2250
  clat=43.25                             # PARAMETERS TO BE SET BY USER
  cdep=300.0  ## depth: km
  tw=160.
  #=========================================================================
  dist=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $1}'`
  baz=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $2}'`
  az=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $3}'`
  $prog/taup_time  -mod prem -h $evdp  -ph P -deg $dist  >taup.out
  #$prog/taup_time  -mod prem -h $evdp  -ph P -deg $dist --stadepth $cdep >taup.out
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
  nlay=`python make_${mod}_fkmodel_150km.py |wc -l`
  echo "NLAYER             $nlay"  >>./$datadir/$ievt/FKmodel
  python make_${mod}_fkmodel_150km.py >>./$datadir/$ievt/FKmodel
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
 ORIGIN_WAVEFRONT  $xx0       $yy0      $zz0
# ORIGIN_TIME 
EOF

  nsta=`cat $datadir/$ievt/$stalst |wc -l`
  #========================================================
  # set it to the same as used in SEM
  dt=0.025                    # PARAMETERS TO BE SET BY USER
  nt=4800
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
is_time_pick true
time_window true 5. 45.
station_coord_system geo
###################################################
STARTLIST
EOF

#source_time_function './DATA/DATA_SYNTHS/$ievt/STF'

  cat $datadir/$ievt/$stalst |
  while read line;do
    stnm=`echo $line |awk '{print $1}'`  
    i=0
    for comp in R T Z;do
      mylist=`echo $stalst |awk -F. '{printf"%s_%s.lst",$1,v}' v=$comp`
      is_data=`grep $stnm $datadir/$ievt/$mylist`
      if [ -z "$is_data" ];then
         dc[$i]=0
      else
         dc[$i]=1
      fi
      let i=i+1
    done  
    echo $stnm ${dc[0]} ${dc[1]} ${dc[2]}
    grep $stnm $datadir/$ievt/$stalst |awk '{print $1,$2,$3,$4,$5,$6,a,b,c}' a=${dc[0]} b=${dc[1]} c=${dc[2]} >>$datadir/$ievt/FK$ievt
  done
  echo ENDLIST >>$datadir/$ievt/FK$ievt
  #---------------------------------------------------------------#
  #-----------------------2. get Source Time Function ---------------#
  #python write_ascii_source_signature.py $outdir/$ievt/fstack.sac.tp.stf $datadir/$ievt/STF
  #-----------------------3 .Set CMTSOLUTION -------------------------#
  evtnm=`echo $dirname |awk '{print $1}' |awk -F/ '{print $3}'`
  #cp $cmtdir/CMTSOLUTION_${evtnm} $datadir/$ievt/CMTSOLUTION
  #----------------------------------------------------#
  let ievt=ievt+1
done
########## END loop over events ###############################
