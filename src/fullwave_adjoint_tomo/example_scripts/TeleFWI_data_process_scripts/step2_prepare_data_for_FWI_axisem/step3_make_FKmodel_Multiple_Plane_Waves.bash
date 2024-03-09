#!/bin/bash
# Author  : Kai Wang
# Inst    : University of Toronto
# Email   : wangkaim8@gmail.com
# history : last modified on Jun 04, 2019


module load anaconda2
module load java
module load gnuplot

prog=~/progs/TauP-2.4.5/bin
zrayamp=~/progs/teleseis_stf/zrayamp/CODE
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
  #---------- Calculate P, pP, PcP and pPcP using zrayamp ------------------
  dist=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $1}'`
  distkm=`echo $dist |awk '{print 6371.0*$1/180.0*3.14159265}'`
  baz=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $2}'`
  az=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $3}'`
  # P
  $prog/taup_time  -mod ak135 -h $evdp  -ph P -deg $dist  >taup.out
  inc_ang1=`cat taup.out |sed -n '6p' |awk '{print $7}'`  ## Incident angle
  pred_t1=`cat taup.out |sed -n '6p' |awk '{print $4}'`  
  # pP
  $prog/taup_time  -mod ak135 -h $evdp  -ph pP -deg $dist  >taup.out
  inc_ang2=`cat taup.out |sed -n '6p' |awk '{print $7}'`  ## Incident angle
  pred_t2=`cat taup.out |sed -n '6p' |awk '{print $4}'`  
  # PcP
  $prog/taup_time  -mod ak135 -h $evdp  -ph PcP -deg $dist  >taup.out
  inc_ang3=`cat taup.out |sed -n '6p' |awk '{print $7}'`  ## Incident angle
  pred_t3=`cat taup.out |sed -n '6p' |awk '{print $4}'`  
  # pPcP
  #$prog/taup_time  -mod prem -h $evdp  -ph pPcP -deg $dist  >taup.out
  #inc_ang4=`cat taup.out |sed -n '6p' |awk '{print $7}'`  ## Incident angle
  echo $dist $distkm $evdp
  ### zrayamp
  zraydir=zrayamp$ievt
  inp=zrayamp.dat
  mkdir -p $zraydir
  ### make data for EVENT 01
  echo "MODEL AK135" >$inp
  echo " 2 0 1" >>$inp   ### K2, MPRINT, KX(1=Two-point ray tracing, 2=initial-value ray tracing)
  python ./make_ak135_fkmodel_for_zrayamp.py >>$inp
  cat >>$inp <<eof
1. .0002
1 -1 1 1 0 1 0 1 0 99 1
$distkm
0. $evdp 0. .1  3.32   .03125
.001 .001 1.5707 -.001 -.001 -1.5707 
1 10  3  4  5  6  6  5  4  3  2  1
1 12 1 2 3 4 5 6 6 5 4 3 2 1
-1 13 1 1 2 3 4 5 6 6 5 4 3 2 1
-1 14 2 1 1 2 3 4 5 6 6 5 4 3 2 1
-1 15 3 2 1 1 2 3 4 5 6 6 5 4 3 2 1
-1 17 5 4 3 2 1 1 2 3 4 5 6 6 5 4 3 2 1
-1 18 6 5 4 3 2 1 1 2 3 4 5 6 6 5 4 3 2 1
0 0
eof
  rm -f zrayamp.out zrayamp.res rays.dat timeampl.dat
  $zrayamp/zrayamp
  $zrayamp/posi
  cp $zrayamp/plot.all.g .
  gnuplot plot.all.g
  mv plot.all.* $zraydir
  mv zrayamp.dat zrayamp.out zrayamp.res rays.dat timeampl.dat $zraydir
  ### Get amplitude
  ampl1=`cat $zraydir/timeampl.dat |awk '{if(sqrt(($2-t)*($2-t))<0.5) print $3}' t=$pred_t1 |uniq` 
  ampl2=`cat $zraydir/timeampl.dat |awk '{if(sqrt(($2-t)*($2-t))<0.5) print $3}' t=$pred_t2 |uniq` 
  ampl3=`cat $zraydir/timeampl.dat |awk '{if(sqrt(($2-t)*($2-t))<0.5) print $3}' t=$pred_t3 |uniq` 
  dt2=`echo $pred_t1 $pred_t2 |awk '{if(($2-$1)<=50.) print 1;else print 0}'`
  dt3=`echo $pred_t1 $pred_t3 |awk '{if(($2-$1)<=50.) print 1;else print 0}'`
  if [ ! -z "$ampl2" ] && [ $dt2 -eq 1 ];then
     is_wave2=1
  else
     is_wave2=0
  fi
  if [ ! -z "$ampl3" ] && [ $dt3 -eq 1 ];then
     is_wave3=1
  else
     is_wave3=0
  fi
  echo $ampl1 $ampl2 $ampl3 $is_wave2 $is_wave3
  if [ -z "$ampl1" ] && [ -z "$ampl2" ] && [ -z "$ampl3" ];then
     echo "No valid amplitude for P, pP, or PcP"
     echo "Check traveltime tables:"
     cat $zraydir/timeampl.dat
     echo TAUP $pred_t1　$pred_t2　$pred_t3
     exit
  fi
  echo "P   ==> time, ampl: " $pred_t1,$ampl1
  echo "pP  ==> time, ampl: " $pred_t2,$ampl2
  echo "PcP ==> time, ampl: " $pred_t3,$ampl3
  if [ -z "$ampl2" ] && [ -z "$ampl3" ];then
     num_wave=1
  elif [ $is_wave2 -eq 1 ] && [ $is_wave3 -eq 1 ];then
     num_wave=3
     relt1=0.0
     relt2=`echo $pred_t1 $pred_t2 |awk '{print $2-$1}'`
     relt3=`echo $pred_t1 $pred_t3 |awk '{print $2-$1}'`
  elif [ $is_wave2 -eq 1 ] && [ $is_wave3 -eq 0 ];then
     num_wave=2
     am1=$ampl1
     am2=$ampl2
     relt1=0.0
     relt2=`echo $pred_t1 $pred_t2 |awk '{print $2-$1}'`
     toff1=$inc_ang1
     toff2=$inc_ang2
  elif [ $is_wave2 -eq 0 ] && [ $is_wave3 -eq 1 ];then
    num_wave=2
    am1=$ampl1
    am2=$ampl3
    relt1=0.0
    relt2=`echo $pred_t1 $pred_t3 |awk '{print $2-$1}'`
    toff1=$inc_ang1
    toff2=$inc_ang3
  fi
  echo "num_wave: " $num_wave
  if [ $num_wave -eq 2 ];then
     echo $am1 $am2 $relt1 $relt2
  elif [ $num_wave -eq 3 ];then
     echo $ampl1 $ampl2 $ampl3 $relt1 $relt2 $relt3
  fi
  #-----------------------------------------#
  if true;then
  #-------------------------------------------------------------------------
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
  nlay=`python make_${mod}_fkmodel.py |wc -l`
  echo "NLAYER             $nlay"  >>./$datadir/$ievt/FKmodel
  python make_${mod}_fkmodel.py >>./$datadir/$ievt/FKmodel
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
EOF
 if [ $num_wave -eq 2 ];then
  cat >>./$datadir/$ievt/FKmodel <<EOF
# for multiple plane wave injection
NUMBER_OF_WAVES     2
PLANE 1  $relt1 $am1 $toff1
PLANE 2  $relt2 $am2 $toff2
EOF
 elif [ $num_wave -eq 3 ];then
  cat >>./$datadir/$ievt/FKmodel <<EOF
# for multiple plane wave injection
NUMBER_OF_WAVES     3
PLANE 1  $relt1 $ampl1 $inc_ang1
PLANE 2  $relt2 $ampl2 $inc_ang2
PLANE 3  $relt3 $ampl3 $inc_ang3
EOF
 fi
  cat >>./$datadir/$ievt/FKmodel <<EOF
#----------------------------------------
# anlges of incomming wave
BACK_AZIMUTH               $baz
TAKE_OFF               $inc_ang1
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
  cp $cmtdir/CMTSOLUTION_${evtnm} $datadir/$ievt/CMTSOLUTION
  fi
  #----------------------------------------------------#

  let ievt=ievt+1
done
########## END loop over events ###############################
