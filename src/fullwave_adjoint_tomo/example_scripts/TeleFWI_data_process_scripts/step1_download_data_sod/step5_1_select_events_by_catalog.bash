#!/bin/bash
#cl=d0-30  # catalog <30km or >200km
#cl=d200-10000  # catalog <30km or >200km
#cl=d30-200  # catalog 30-200km for additional direct P
#cl=P_XJ1997
#cl=P_XE2005-2007
cl=S_TO2013-2015
out=final_events_${cl}.lst

cat /dev/null >$out

#cat event_mag5.8_dep200kmup_TO_dist30-90.lst |
#while read line;do 
i=1
IFS=$'\n'
for line in `cat event_mag5.5_${cl}.lst`;do
#for line in `cat event_mag5.8_${cl}.lst`;do
evlo=`echo $line |awk '{print $1}'`
evla=`echo $line |awk '{print $2}'`
evdp=`echo $line |awk '{print $3}'`
evnm=`echo $line |awk '{print $4}'`
evdir=`echo $evnm |awk -F_ '{printf"%04d_%02d_%02d_%02d_%02d_%02d",$1,$2,$3,$4,$5,$6}'`
echo "$i: evlo evla evdp evdir  $evlo $evla $evdp $evdir"

 echo evt_catalog_${cl}/processedSeismograms/Event_${evdir} $evlo $evla $evdp >>$out
#if [ -f evt_catalog_${cl}/processedSeismograms/seismo_Event_${evdir}.ps ];then
#  gv -scale=0.6  evt_catalog_${cl}/processedSeismograms/seismo_Event_${evdir}.ps &
#  echo "Save this event or Not? Yes(y)No(n)"
#  read issave
#  if [ $issave == "y" ];then
#    echo evt_catalog_${cl}/processedSeismograms/Event_${evdir} $evlo $evla $evdp >>$out
#  fi
#fi
#
#let i=i+1
#if [ $i -eq 20 ];then
#   killall gv
#fi
done

killall gv
