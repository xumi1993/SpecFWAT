#!/bin/bash

cd DATA/DATA_SYNTHS
cat /dev/null >station_all.lst
for evid in `ls -d *`;do
   echo $evid
   cat $evid/stations_pick.lst >> ../station_all.lst
done
cd ../../

