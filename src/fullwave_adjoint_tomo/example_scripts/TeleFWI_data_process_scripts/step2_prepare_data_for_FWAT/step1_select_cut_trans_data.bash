#!/bin/bash

prog=~/progs/TauP-2.4.5/bin
sactool=~/progs/sactools_c


indir=../step1_download_data_sod
outdir=selectedSeismograms
#datadir=DATA/DATA_SYNTHS
datadir=src_rec
evtlst=final_events_S_TO2013-2015.lst
phase=S

#rm -rf $outdir $datadir
mkdir -p $outdir
mkdir -p $datadir
cat /dev/null >$datadir/sources.dat

ievt=1
IFS=$'\n'
#################### Main loop over events #####################
for evtdir in `cat $indir/$evtlst`;do
  echo $evtdir
  dirname=`echo $evtdir |awk '{print $1}'`
  evlo=`echo $evtdir |awk '{print $2}'`
  evla=`echo $evtdir |awk '{print $3}'`
  evdp=`echo $evtdir |awk '{print $4}'` # km
  echo $ievt $evla $evlo $evdp
  mkdir -p $outdir/$phase$ievt
  #------ 1. Make acqui_file.txt -------#
  echo "$phase$ievt" >>$datadir/sources.dat
  #----------------------------------#
  ################## Loop over stations #################################
  #--------  2. copy selected seismograms --------------------#
  cat /dev/null >time.pred
  mkdir $datadir/$phase$ievt
  cat /dev/null >$datadir/$phase$ievt/stations.lst
  for file in `ls $indir/$dirname/*.?HZ.sac`;do
    #fname=`echo $file |awk -F$indir/$dirname/ '{print $2}'`
    sac <<EOF
r $file
wh
quit
EOF
    b=`saclst b f $file |awk '{print $2}'`  ## event origin time - 120 s
    dt=`saclst delta f $file |awk '{print $2}'`
    stla=`saclst stla f $file |awk '{print $2}'`
    stlo=`saclst stlo f $file |awk '{print $2}'`
    stnm=`saclst kstnm f $file |awk '{print $2}'`
    netwk=`saclst knetwk f $file |awk '{print $2}'`
    stel=`saclst stel f $file |awk '{print $2}'`
    dist=`saclst dist f $file |awk '{printf"%.5f",$2}'`
    az=`saclst az f $file |awk '{print $2}'`
    baz=`saclst baz f $file |awk '{print $2}'`
    cmp=`saclst kcmpnm f $file |awk '{print substr($2,1,2)}'`
    khole=`saclst khole f $file |awk '{print $2}'`
   

    #-----------------------------------------------------------------------
    echo "$prog/taup_time -mod prem -h $evdp  -ph $phase -km $dist"
    $prog/taup_time -mod prem -h $evdp  -ph $phase -km $dist >taup.out
    t_ref_taup=`cat taup.out |sed -n '6p' |awk '{print $4}'`  ## referece P/S arrival
    tmin=`echo $t_ref_taup |awk '{print $1-120.}'`  ## referece P/S arrival - 120. s
    tmax=`echo $t_ref_taup |awk '{print $1+179.95}'`  ## referece P/S arrival + 180. s
    inc_ang=`cat taup.out |sed -n '6p' |awk '{print $7}'`  ## Incident angle
    tdiff=`echo $b $t_ref_taup |awk '{print sqrt(($1+120-$2)*($1+120-$2))}'`
    echo $stnm $netwk $stla $stlo $stel $b $t_ref_taup $tdiff >>time.pred
    #is_small=`echo $tdiff |awk '{if($1< 0.5) print 1;else print 0}'`
    #if [ $is_small -eq 0 ];then
    #   echo "Attention!!! reference time difference is too large"
    #fi
    #echo "==> $t_ref_taup $tdiff $inc_ang"
    #-----------------------------------------------------------------------

    file_e=$indir/$dirname/$netwk.$stnm..${cmp}E.sac
    file_n=$indir/$dirname/$netwk.$stnm..${cmp}N.sac
    file_r=$indir/$dirname/$netwk.$stnm..${cmp}R.sac
    file_t=$indir/$dirname/$netwk.$stnm..${cmp}T.sac

    sacf=$outdir/$phase$ievt/$netwk.$stnm..BHZ.sac
    sacf_e=$outdir/$phase$ievt/$netwk.$stnm..BHE.sac
    sacf_n=$outdir/$phase$ievt/$netwk.$stnm..BHN.sac
    sacf_r=$outdir/$phase$ievt/$netwk.$stnm..BHR.sac
    sacf_t=$outdir/$phase$ievt/$netwk.$stnm..BHT.sac
    ### Check data here ###
    if [ -f $file ] && [ -f $file_e ] && [ -f $file_n ];then
       echo "Roatation N, E to R, T"
       sac <<EOF
r $file $file_e $file_n
ch b -120.
ch t0 0.
w $sacf $sacf_e $sacf_n
r $sacf_n $sacf_e
rot to gcp
ch t0 0.
bp n 4 p 2 c 0.01 0.2
w $sacf_r $sacf_t
quit
EOF
       sacz_zero=`$sactool/sac_zero $sacf |awk '{print $2}'`
       sace_zero=`$sactool/sac_zero $sacf_e |awk '{print $2}'`
       sacn_zero=`$sactool/sac_zero $sacf_n |awk '{print $2}'`
       echo $ievt: $stnm $sacz_zero $sace_zero $sacn_zero
       if [ $sacz_zero -eq 0 ] && [ $sace_zero -eq 0 ] && [ $sacn_zero -eq 0 ] ;then
          echo  $stnm $netwk $stla $stlo $stel 0.00000 >>$datadir/$phase$ievt/stations.lst
       fi
    fi
    
    if [ -f $file ] && [ -f $file_r ] && [ -f $file_t ];then
#============ For resampling the SAC files ================
       delta=`saclst delta f $file |awk '{printf"%f\n",$2}'`
       is_dec=`echo $delta |awk '{if(0.05/$1>1) {print 10} else if(0.05/$1==1) {print 1} else {print 0}}'`       
       dec_mod=`echo $delta |awk '{if(0.05%$1>0) {print 1} else {print 0}}'`       
       dec_num=`echo $delta |awk '{printf"%d",0.05/$1+0.5}'`       
       echo "delta is_dec dec_mod dec_num:" $delta $is_dec $dec_mod $dec_num
       if [ $is_dec -gt 1 ] && [ $dec_mod -ne 0 ];then  ### 
          echo "Decimate by $dec_num and interp 0.05 " 
          sac <<EOF
cuterr fillz
cut $tmin $tmax
r $file $file_r $file_t
ch b -120.
ch t0 0.
decimate $dec_num
interp delta 0.05
bp n 4 p 2 c 0.02 0.2
w $sacf $sacf_r $sacf_t
quit
EOF
       elif [ $is_dec -gt 1 ] && [ $dec_mod -eq 0 ];then  ###
          echo "Decimate by $dec_num" 
          sac <<EOF
cuterr fillz
cut $tmin $tmax
r $file $file_r $file_t
ch b -120.
ch t0 0.
decimate $dec_num
bp n 4 p 2 c 0.02 0.2
w $sacf $sacf_r $sacf_t
quit
EOF
       elif [ $is_dec -lt 1 ];then
          echo "Interp 0.05"
          sac <<EOF
cuterr fillz
cut $tmin $tmax
r $file $file_r $file_t
ch b -120.
ch t0 0.
interp delta 0.05
bp n 4 p 2 c 0.02 0.2
w $sacf $sacf_r $sacf_t
quit
EOF
       elif [ $is_dec -eq 1 ];then
	  echo "No decimate or interp"
          sac <<EOF
cuterr fillz
cut $tmin $tmax
r $file $file_r $file_t
ch b -120.
ch t0 0.
bp n 4 p 2 c 0.02 0.2
w $sacf $sacf_r $sacf_t
quit
EOF
      fi
#============================================================

       sacz_zero=`$sactool/sac_zero $sacf |awk '{print $2}'`
       sacr_zero=`$sactool/sac_zero $sacf_r |awk '{print $2}'`
       sact_zero=`$sactool/sac_zero $sacf_t |awk '{print $2}'`
       echo $ievt: $stnm $cmp $sacz_zero $sacr_zero $sact_zero
       if [ $sacz_zero -eq 0 ] && [ $sacr_zero -eq 0 ] && [ $sact_zero -eq 0 ] ;then
          echo  $stnm $netwk $stla $stlo $stel 0.00000 >>$datadir/$phase$ievt/stations.lst
       else
	  echo "zero component found !!!"
          rm $sacf $sacf_r $sacf_t
       fi
   fi

  done
  ################### END loop over stations ######################
  #----------------------------------------------------#
  let ievt=ievt+1
done
########## END loop over events ###############################
