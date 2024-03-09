#!/bin/bash

#module load anaconda2
#module load java

prog=~/progs/TauP-2.4.5/bin
sactool=~/progs/sactools_c


indir=../step1_download_data_sod
outdir=selectedSeismograms
datadir=DATA/DATA_SYNTHS
evtlst=final_events.lst

#rm -rf $outdir $datadir
mkdir -p $outdir
mkdir -p $datadir
cat /dev/null >acqui_file.txt

ievt=1
IFS=$'\n'
#################### Main loop over events #####################
for evtdir in `cat ${indir}/Picked_8events_20190801/$evtlst`;do
  dirname=`echo $evtdir |awk '{print $1}'`
  evlo=`echo $evtdir |awk '{print $2}'`
  evla=`echo $evtdir |awk '{print $3}'`
  evdp=`echo $evtdir |awk '{print $4}'` # km
  #evt=`echo $dirname |awk -F/ '{print $2}' |awk -F_ '{print $1}'`
  echo $ievt $evla $evlo $evdp
  mkdir -p $outdir/$ievt
  #------ 1. Make acqui_file.txt -------#
  echo "event_rep : ./$datadir/$ievt" >>acqui_file.txt
  echo "event_name : FK$ievt" >>acqui_file.txt
  #----------------------------------#
  ################## Loop over stations #################################
  #--------  2. copy selected seismograms --------------------#
  cat /dev/null >$datadir/$ievt/time.pred
  mkdir $datadir/$ievt
  cat /dev/null >$datadir/$ievt/stations.lst
  for file in `ls $indir/$dirname/*HHZ.sac`;do
     echo $file
    fname=`echo $file |awk -F$indir/$dirname/ '{print $2}'`

    b=`$sactool/saclst b f $file |awk '{print $2}'`  ## event origin time - 120 s
    dt=`$sactool/saclst delta f $file |awk '{print $2}'`
    stla=`saclst stla f $file |awk '{print $2}'`
    stlo=`saclst stlo f $file |awk '{print $2}'`
    stnm=`saclst kstnm f $file |awk '{print $2}'`
    netwk=`saclst knetwk f $file |awk '{print $2}'`
    stel=`saclst stel f $file |awk '{print $2}'`
    dist=`saclst dist f $file |awk '{printf"%.5f",$2}'`
    az=`saclst az f $file |awk '{print $2}'`
    baz=`saclst baz f $file |awk '{print $2}'`
    khole=`saclst khole f $file |awk '{print $2}'`
    khole=''

    sacf=$outdir/$ievt/$netwk.$stnm.$khole.BHZ.sac


    #-----------------------------------------------------------------------
    $prog/taup_time -mod ak135 -h $evdp  -ph P -km $dist >taup.out
    t_ref_taup=`cat taup.out |sed -n '6p' |awk '{print $4}'`  ## referece P arrival
    tmin=`echo $t_ref_taup |awk '{print $1-120.}'`  ## referece P arrival - 120. s
    tmax=`echo $t_ref_taup |awk '{print $1+179.95}'`  ## referece P arrival + 180. s
    echo $stnm $netwk $stla $stlo $stel $t_ref_taup >>$datadir/$ievt/time.pred
    #----------------------------------------------------------------------
    file_e=$indir/$dirname/$netwk.$stnm.$khole.HHE.sac
    file_n=$indir/$dirname/$netwk.$stnm.$khole.HHN.sac
    file_r=$indir/$dirname/$netwk.$stnm.$khole.HHR.sac
    file_t=$indir/$dirname/$netwk.$stnm.$khole.HHT.sac
    sacf_e=$outdir/$ievt/$netwk.$stnm.$khole.BHE.sac
    sacf_n=$outdir/$ievt/$netwk.$stnm.$khole.BHN.sac
    sacf_r=$outdir/$ievt/$netwk.$stnm.$khole.BHR.sac
    sacf_t=$outdir/$ievt/$netwk.$stnm.$khole.BHT.sac
    ### Check data here ###
    if [ -f $file ] && [ -f $file_e ] && [ -f $file_n ];then
       sac <<EOF
r $file $file_e $file_n
ch LCALDA true
bp n 4 p 2 c 0.01 0.2
w $sacf $sacf_e $sacf_n
r $sacf_n
ch CMPAZ 0 CMPINC 90
wh
r $sacf_e
ch CMPAZ 90 CMPINC 90
wh
r $sacf_n $sacf_e
rot to gcp
bp n 4 p 2 c 0.01 0.2
w $sacf_r $sacf_t
quit
EOF
       sacz_zero=`$sactool/sac_zero $sacf |awk '{print $2}'`
       sace_zero=`$sactool/sac_zero $sacf_e |awk '{print $2}'`
       sacn_zero=`$sactool/sac_zero $sacf_n |awk '{print $2}'`
       echo $ievt: $stnm $sacz_zero $sace_zero $sacn_zero
       if [ $sacz_zero -eq 0 ] && [ $sace_zero -eq 0 ] && [ $sacn_zero -eq 0 ] ;then
          echo  $stnm $netwk $stla $stlo $stel 0.00000 >>$datadir/$ievt/stations.lst
       fi
    fi
    if [ -f $file ] && [ -f $file_r ] && [ -f $file_t ];then
       sac <<EOF
cuterr fillz
cut $tmin $tmax
r $file $file_r $file_t
ch b -120.
ch t0 0.
decimate 5
bp n 4 p 2 c 0.01 0.2
w $sacf $sacf_r $sacf_t
quit
EOF

       sacz_zero=`$sactool/sac_zero $sacf |awk '{print $2}'`
       sacr_zero=`$sactool/sac_zero $sacf_r |awk '{print $2}'`
       sact_zero=`$sactool/sac_zero $sacf_t |awk '{print $2}'`
       echo $ievt: $stnm $sacz_zero $sacr_zero $sact_zero
       if [ $sacz_zero -eq 0 ] && [ $sacr_zero -eq 0 ] && [ $sact_zero -eq 0 ] ;then
          echo  $stnm $netwk $stla $stlo $stel 0.00000 >>$datadir/$ievt/stations.lst
       fi
    fi
  done
  ################### END loop over stations ######################
  #----------------------------------------------------#
  let ievt=ievt+1
done
########## END loop over events ###############################
