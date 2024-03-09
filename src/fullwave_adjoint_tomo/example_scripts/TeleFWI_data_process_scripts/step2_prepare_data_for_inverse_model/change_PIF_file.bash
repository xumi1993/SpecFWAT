#!/bin/bash


cd DATA/DATA_SYNTHS
for evid in `ls -d *`;do
   echo $evid
   sed -i "/data_time_step/c\data_time_step 0.025" $evid/PIF$evid
   sed -i "/data_sample_number/c\data_sample_number 4800" $evid/PIF$evid
done
cd ../../
