#!/bin/bash

### Step1: find events ###
#find_events -d36/-120/30/90 -b 2013-12-13 -e 2015-10-12 -m 5.8  -D0-30 >event_mag5.8_dep30kmless_TO_dist30-90.lst
#find_events -d36/-120/30/90 -b 2013-12-13 -e 2015-10-12 -m 5.8  -D0-30 -r > recipes/Cecal/event1.xml
#find_events -d36/-120/30/90 -b 2013-12-13 -e 2015-10-12 -m 5.8  -D200-10000 >event_mag5.8_dep200kmup_TO_dist30-90.lst
#find_events -d36/-120/30/90 -b 2013-12-13 -e 2015-10-12 -m 5.8  -D200-10000 -r > recipes/Cecal/event2.xml
#find_events -d36/-120/30/90 -b 2013-12-13 -e 2015-10-12 -m 5.8  -D30-200 >event_mag5.8_d30-200_TO_dist30-90.lst
#find_events -d36/-120/30/90 -b 2013-12-13 -e 2015-10-12 -m 5.8  -D30-200 -r > recipes/Cecal/event3.xml
### Step2: find stations ###
#find_stations -R-122/-118/35/37 -n TO >station_TO.lst
find_stations -R-122.2/-117.5/35.0/37.25 -n TO,BK,CE,CI >station_Cecal.lst
### Step3: download CMTSOLUTION ###
# go to https://www.globalcmt.org/CMTsearch.html
