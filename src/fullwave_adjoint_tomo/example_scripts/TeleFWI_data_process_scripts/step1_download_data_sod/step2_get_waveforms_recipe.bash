#!/bin/bash

### for Teleseismic FWI by hybrid methods (Pick good waveform of P at both R and Z components)
#  Mag >5.8
#  Dep < 30 km or > 200 km
#  Dist >30 deg < 90 deg

#find_events -d36/-120/30/90 -b 2013-12-13 -e 2015-10-12 -m 5.8  -D200-10000 -r |
#find_stations -R-122/-118/35/37 -n TO -r |
#find_seismograms -B -2ttp -E 3ttp -c HH* -r >recipes/Cecal/waveform_d200+.xml 
#
#find_events -d36/-120/30/90 -b 2013-12-13 -e 2015-10-12 -m 5.8  -D0-30 -r |
#find_stations -R-122/-118/35/37 -n TO -r |
#find_seismograms -B -2ttp -E 3ttp -c HH* -r >recipes/Cecal/waveform_d0-30.xml 

### for Teleseismic travetime inversion by hybrid methods (pick direct P at Z comp and direct S at T comp)
#   Mag >5.5
#   Dist >30 deg


### addtional data for XJ,CI,BK for year 1997
#recipe=recipes/Cecal/waveform_XJ_CI_BK_1997.xml
##find_events -d36/-120/30/90 -b 1997-06-01 -e 1997-12-31 -m 5.5 >event_mag5.5_P_XJ1997.lst
#find_events -d36/-120/30/90 -b 1997-06-01 -e 1997-12-31 -m 5.5  -D0-1000 -r |
#find_stations -R-122.2/-117.5/35/38 -n CI,NN,NC -r >$recipe        # SCEDC
##find_stations -R-122.2/-117.5/35/38 -n XJ,CI,BK -r >$recipe  # IRIS
#sed -i '$ d' $recipe
#sed -i '$ d' $recipe
#cat >> $recipe <<!
#    <waveformVectorArm>
#        <phaseRequest>
#            <beginPhase>ttp</beginPhase>
#            <beginOffset>
#                <unit>MINUTE</unit>
#                <value>-2</value>
#            </beginOffset>
#            <endPhase>ttp</endPhase>
#            <endOffset>
#                <unit>MINUTE</unit>
#                <value>3</value>
#            </endOffset>
#        </phaseRequest>
#        <fdsnDataSelect>
#            	<host>service.scedc.caltech.edu</host>
#        </fdsnDataSelect>
#	<bestChannelAtStation/>
#        <fullCoverage/>
#        <printlineSeismogramProcess/>
#        <sacWriter/>
#        <responseGain/>
#        <rMean/>
#        <rTrend/>
#        <integrate/>
#                <sampleSyncronize/>
#                <vectorTrim/>
#                <rotateGCP/>
#        <sacWriter>
#            <workingDir>processedSeismograms</workingDir>
#        </sacWriter>
#        <legacyExecute>
#            <command>echo Sod saved this file</command>
#        </legacyExecute>
#
#    </waveformVectorArm>
#</sod>
#!
# rm -rf Sod*
# sod -f $recipe

######## Additional data for XE, TA, CI for year 2005-2007
recipe=recipes/Cecal/waveform_XE_TA_CI_2005-2007.xml
#find_events -d36/-120/30/90 -b 2005-06-01 -e 2007-12-31 -m 5.8 >event_mag5.8_P_XE2005-2007.lst
find_events -d36/-120/30/90 -b 2005-06-01 -e 2007-12-31 -m 5.8  -D0-1000 -r |
find_stations -R-122.2/-117.5/35/38 -n CI,NC,NN -r >$recipe  ### SCEDC
#find_stations -R-122.2/-117.5/35/38 -n XE,CI,TA,BK -r >$recipe  ### IRIS
sed -i '$ d' $recipe
sed -i '$ d' $recipe
cat >> $recipe <<!
    <waveformVectorArm>
        <phaseRequest>
            <beginPhase>ttp</beginPhase>
            <beginOffset>
                <unit>MINUTE</unit>
                <value>-2</value>
            </beginOffset>
            <endPhase>ttp</endPhase>
            <endOffset>
                <unit>MINUTE</unit>
                <value>3</value>
            </endOffset>
        </phaseRequest>
        <fdsnDataSelect>
        	<host>service.scedc.caltech.edu</host>
        </fdsnDataSelect>
	<bestChannelAtStation/>
        <fullCoverage/>
        <printlineSeismogramProcess/>
        <sacWriter/>
        <responseGain/>
        <rMean/>
        <rTrend/>
        <integrate/>
                <sampleSyncronize/>
                <vectorTrim/>
                <rotateGCP/>
        <sacWriter>
            <workingDir>processedSeismograms</workingDir>
        </sacWriter>
        <legacyExecute>
            <command>echo Sod saved this file</command>
        </legacyExecute>

    </waveformVectorArm>
</sod>
!
 rm -rf Sod*
 sod -f $recipe


