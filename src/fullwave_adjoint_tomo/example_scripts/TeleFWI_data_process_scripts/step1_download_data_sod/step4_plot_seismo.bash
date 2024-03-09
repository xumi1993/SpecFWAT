#!/bin/bash
#cd evt_catalog_d0-30/processedSeismograms
#cd evt_catalog_d200-10000/processedSeismograms
#cd evt_catalog_d30-200/processedSeismograms
#cd evt_catalog_P_XJ1997/processedSeismograms
cd evt_catalog_P_XE2005-2007/processedSeismograms

for evtdir in `ls -d Event_*`;do
  out=seismo_${evtdir}.ps
  icomp=1
  for comp in BHZ BHR BHT;do
    echo $icomp $comp
    sac <<eof
r $evtdir/*.$comp.sac
bp n 4 p 2 c 0.02 0.2
w append .lp
quit
eof
    saclst b e dist f $evtdir/*.$comp.sac.lp >dist.dat 
    t_min=`cat dist.dat |minmax -C |awk '{print $3}'`
    t_max=`cat dist.dat |minmax -C |awk '{print $6}'`
    d_min=`cat dist.dat |minmax -C |awk '{print $7-10}'`
    d_max=`cat dist.dat |minmax -C |awk '{print $8+10}'`
    tx=`echo $t_min $t_max |awk '{print ($1+$2)/2}'`
    ty=`echo $d_min $d_max |awk '{print $2+($2-$1)*0.04}'`
    echo $tx $ty 

    if [ $icomp -eq 1 ];then
       psbasemap -JX16/8 -R$t_min/$t_max/$d_min/$d_max -Ba50/a10:"Distance (km)":WeSn -K -P -X4 -Y20 >$out
    pstext -J -R -N -O -K >>$out <<eof
$tx $ty 16 0 0 CB $evtdir
eof
    elif [ $icomp -eq 2 ];then
       psbasemap -JX16/8 -R$t_min/$t_max/$d_min/$d_max -Ba50/a10:"Distance (km)":Wesn -K -O -Y-8.75 >>$out
    else
       psbasemap -JX16/8 -R$t_min/$t_max/$d_min/$d_max -Ba50:"Time (s)":/a10:"Distance (km)":WeSn -K -O -Y-8.75 >>$out
    fi

    pssac $evtdir/*.$comp.sac.lp -J -R -Ektb -M1. -W0.5p -K -O >>$out 
    rm $evtdir/*.$comp.sac.lp
    let icomp=icomp+1
  done
  cat /dev/null |psxy -J -R -Sc0.1i -O >>$out
done
cd ../..
