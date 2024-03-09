#!/bin/bash

cat /dev/null >events_zz.lst
IFS=$'\n'
for line in `cat final_events_directP.lst`;do
    evtpath=`echo $line |awk '{print $1}'`
    evnm=`echo $evtpath |awk -F/ '{print $3}'`
    year=`echo $evnm |awk -F_ '{print $2}'`
    month=`echo $evnm |awk -F_ '{print $3}'`
    day=`echo $evnm |awk -F_ '{print $4}'`
    hr=`echo $evnm |awk -F_ '{print $5}'`
    sec=`echo $evnm |awk -F_ '{print $6}'`
    msec=`echo $evnm |awk -F_ '{print $7}'`
    echo $evnm
    echo "${year}${month}${day}${hr}${sec}${msec}" >>events_zz.lst 
done

