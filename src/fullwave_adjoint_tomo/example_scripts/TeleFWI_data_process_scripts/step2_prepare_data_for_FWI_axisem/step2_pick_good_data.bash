#!/bin/bash

prog=~/progs/TauP-2.4.5/bin
sactool=~/progs/sac_msc
astack=~/progs/astack  # adapted stacking from Rawlinson to align waveforms


indir=../step1_download_data_sod
outdir=selectedSeismograms
datadir=DATA/DATA_SYNTHS
evtlst=final_events.lst


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

  #----------------------------------------------------#
  #-------------2. Pick good data and their array stacking seismogram by AIBAT-----------------------------#
  if false;then
  # I did AIMBAT manually in a seperate terminal.
  # As in the FK-SEM code, the direct P arrival of synthetics is at the peak of wave, however the predict P 
  # from reference Earth model is the time of the first motion of P. So, there are time delays between FK time
  # and the predict P. To correct this, my idea is to ...  
  cd $outdir/$ievt
  sac2pkl.py *HHZ.sac -s -o sac.HHZ.sac
  ttpick.py -p P sac.HHZ.pkl -t -5 45 -x -10 70
  sac2pkl.py -p sac.pkl
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
  ########
  for comp in R Z;do
  cat /dev/null >$datadir/$ievt/stations_pick_${comp}.lst
  cat $datadir/$ievt/stations.lst |
  while read line;do
    stnm=`echo $line |awk '{print $1}'`
    ispick=`grep $stnm $outdir/$ievt/*mcp`
    if [ -n "$ispick" ];then
      echo $line >>$datadir/$ievt/stations_pick_${comp}.lst
    else
      echo $stnm is not picked
    fi
  done
  done
  # merge station list of different components without duplicates 
  sort -u stations_pick_R.lst stations_pick_Z.lst >stations_pick.lst 
  ###
  fi


  let ievt=ievt+1
done
########## END loop over events ###############################
