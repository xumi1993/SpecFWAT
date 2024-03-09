#!/bin/bash
# Author  : Kai Wang
# Inst    : University of Toronto
# Email   : wangkaim8@gmail.com
# history : last modified on Jun 04, 2019

#module load anaconda2
#module load java
#module load anaconda3



prog=~/progs/TauP-2.4.5/bin
sactool=~/progs/sactools_c
indir=../step1_download_data_sod/Picked_8events_20190801
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
  #=========================================================================
  clon=-119.95
  clat=35.975                             # PARAMETERS TO BE SET BY USER
  tw=160.
  #=========================================================================
  dist=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $1}'`
  baz=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $2}'`
  az=`$sactool/distaz $clat $clon $evla $evlo |awk '{print $3}'`
 
  #----------------------------------#
  nsta=`cat $datadir/$ievt/$stalst |wc -l`
  #========================================================
  # set it to the same as used in SEM
  dt=0.05                    # PARAMETERS TO BE SET BY USER
  nt=2400
  #========================================================
  cat >$datadir/$ievt/PIF$ievt <<EOF
event_name PIF$ievt
source_type moment
source_components './DATA/DATA_SYNTHS/$ievt/CMTSOLUTION'
modeling_tool 'axisem' './DATA/AxiSEM_tractions/$ievt'
data_components d rtz
cartloc_mesh_origin $clat $clon 0.
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
    grep $stnm $datadir/$ievt/$stalst |awk '{print $1,$2,$3,$4,$5,$6,a,b,c}' a=${dc[0]} b=${dc[1]} c=${dc[2]} >>$datadir/$ievt/PIF$ievt
  done
  echo ENDLIST >>$datadir/$ievt/PIF$ievt
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
