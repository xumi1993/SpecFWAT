#!/bin/bash
# This script is to download specific events when there is a event list.

evtlst=events_tt.lst

cat /dev/null >events_tele_S_tt.lst
for line in `cat $evtlst`;do
   echo $line
   year=`echo $line |awk '{printf"%04d\n",substr($1,1,4)}'`
   month=`echo $line |awk '{printf"%02d\n",substr($1,5,2)}'`
   day=`echo $line |awk '{printf"%02d\n",substr($1,7,2)}'`
   hr=`echo $line |awk '{printf"%02d\n",substr($1,9,2)}'`
   day1=`echo $line |awk '{printf"%02d\n",substr($1,7,2)+1}'`
   hr1=`echo $line |awk '{printf"%02d\n",substr($1,9,2)+1}'`
   echo $year $month $day $hr == $day1 $hr1

   find_events -d36/-120/30/180 -b ${year}-${month}-${day}-${hr} -e ${year}-${month}-${day}-${hr1} -m 5.5  -D0-1000 >> events_tele_S_tt.lst 
   
done


