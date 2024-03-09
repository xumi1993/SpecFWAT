#!/bin/bash

prog=~/progs/TauP-2.4.5/bin
sactool=~/progs/sac_msc
astack=~/progs/astack  # adapted stacking from Rawlinson to align waveforms


indir=../step1_download_data_sod
outdir=selectedSeismograms
datadir=src_rec
phase=S
evtlst=final_events_S_TO2013-2015.lst

ievt=44 ### for TO
#ievt=44 ### for XJ
#ievt=59 ### for XE
IFS=$'\n'
#################### Main loop over events #####################
#for evtdir in `cat $indir/$evtlst`;do
for evtdir in `cat $indir/$evtlst |sed -n '44,$p'`;do
  dirname=`echo $evtdir |awk '{print $1}'`
  evlo=`echo $evtdir |awk '{print $2}'`
  evla=`echo $evtdir |awk '{print $3}'`
  evdp=`echo $evtdir |awk '{print $4}'` # km
  evt=`echo $dirname |awk -FEvent_ '{print $2}'`
  echo $ievt: $evt $evla $evlo $evdp

  #----------------------------------------------------#
  #-------------2. Pick good data and their array stacking seismogram by AIBAT-----------------------------#
  if true;then
  # I did AIMBAT manually in a seperate terminal.
  # As in the FK-SEM code, the direct P arrival of synthetics is at the peak of wave, however the predict P 
  # from reference Earth model is the time of the first motion of P. So, there are time delays between FK time
  # and the predict P. To correct this, my idea is to ...  
  cd $outdir/$phase$ievt
  aimbat-sac2pkl *BHZ.sac -s -o sac.BHZ.pkl
  aimbat-sac2pkl *BHT.sac -s -o sac.BHT.pkl
  aimbat-sac2pkl *BHR.sac -s -o sac.BHR.pkl
  #ttpick.py -p P sac.BHZ.pkl -t -5 45 -x -10 70
  #mv 2*.mcp BHZ.mcp
  #ttpick.py -p P sac.BHR.pkl -t -5 45 -x -10 70
  #mv 2*.mcp BHR.mcp
  #sac2pkl.py -p sac.pkl
  cd ../..
  # How to pick the direct P arrival?
  # 1. Delete bad waveforms by click them (CC<0.90)
  # 2. Press ICCS-A/Align
  # 3. After hitting the Align button, place the cursor on the array stack where the peak of direct P,
  #    either up or down, occurs. Press t and 2 simultaneously on the keyboard to select the arrival time.
  #    Now press Sync. Use the mouse to drag and select the desired time window on the seismogram
  #    on the array stack. Press w key to save the window.
  # 4. Now press refine.
  # 5. Finalize and save header.
  fi
  ########
  if false;then
  #for comp in Z R;do
  for comp in Z T;do
  cat /dev/null >$datadir/$ievt/stations_pick_${comp}.lst
  cat $datadir/$ievt/stations.lst |
  while read line;do
    stnm=`echo $line |awk '{print $1}'`
    ispick=`grep $stnm $outdir/$ievt/2*.BH${comp}.mcp`
    if [ -n "$ispick" ];then
      echo $line >>$datadir/$ievt/stations_pick_${comp}.lst
    else
      echo $stnm is not picked
    fi
  done
  done
  fi
  ### Or choose all data and pick it after forward simulation with STF ###
  #cp $datadir/$ievt/stations.lst $datadir/$ievt/stations_pick_Z.lst
  #cp $datadir/$ievt/stations.lst $datadir/$ievt/stations_pick_R.lst
  ###
  #touch $datadir/$ievt/stations_pick_T.lst ### NOT use T component
  #touch $datadir/$ievt/stations_pick_R.lst ### NOT use R component
  #sort -u $datadir/$ievt/stations_pick_R.lst $datadir/$ievt/stations_pick_Z.lst >$datadir/$ievt/stations_pick.lst 
  ###


  let ievt=ievt+1
done
########## END loop over events ###############################
