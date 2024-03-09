#!/bin/bash
#cmtlst=GCMT_event_mag5.8_all.lst
cmtlst=GCMT.ndk
type=S

#=====================================================
if true;then
#rm -rf GCMT_solutions_${type}
mkdir -p GCMT_solutions_${type}
###@ get CMTSOLUTION for final events @###
#cat final_events_P_XJ1997.lst |
#cat final_events_P_XE2005-2007.lst |
cat final_events_S_TO2013-2015.lst |
while read line;do
evtnm=`echo $line |awk '{print $1}' |awk -F/ '{print $3}'`
y=`echo $evtnm |awk -F_ '{print $2}'`
m=`echo $evtnm |awk -F_ '{print $3}'`
d=`echo $evtnm |awk -F_ '{print $4}'`
min=`echo $evtnm |awk -F_ '{print $5}'`
sec=`echo $evtnm |awk -F_ '{print $6}'`
msec=`echo $evtnm |awk -F_ '{print $7}'`
echo $evtnm

### for CMT SOLUTION format
#key=`echo $y $m $d $min $sec $msec |awk '{printf" PDEW%4d%3d%3d%3d%3d%3d\n",$1,$2,$3,$4,$5,$6}'`
#echo $key
#is_cmt=`grep -n "$key" $cmtlst`
#if [ -z "$is_cmt" ];then
#  echo "NO CMT FOUND, exit !"
#  exit
#fi
#ln1=`grep -n "$key" $cmtlst |awk -F: '{print $1}'`
#ln11=`grep -n "$key" $cmtlst |awk -F: '{print $1+1}'`
#ln2=`echo $ln1 |awk '{print $1+12}'`
#cat $cmtlst |sed -n "${ln1},${ln2}p" >GCMT_solutions_${type}/CMTSOLUTION_${evtnm}
### for ndk format
# PDE  2005/01/01 01:20:05.4  13.78  -88.78 193.1 5.0 0.0 EL SALVADOR
# C200501010120A   B:  4    4  40 S: 27   33  50 M:  0    0   0 CMT: 1 TRIHD:  0.6
# CENTROID:     -0.3 0.9  13.76 0.06  -89.08 0.09 162.8 12.5 FREE S-20050322125201
# 23  0.838 0.201 -0.005 0.231 -0.833 0.270  1.050 0.121 -0.369 0.161  0.044 0.240
# V10   1.581 56  12  -0.537 23 140  -1.044 24 241   1.312   9 29  142 133 72   66
key=`echo $y $m $d $min $sec $msec |awk '{printf"%04d/%02d/%02d %02d:%02d\n",$1,$2,$3,$4,$5}'`
echo $key
is_cmt=`grep -n "$key" $cmtlst`
if [ -z "$is_cmt" ];then
  echo "NO CMT FOUND, exit !"
else
ln1=`grep -n "$key" $cmtlst |awk -F: '{print $1}'`
ln2=`grep -n "$key" $cmtlst |awk -F: '{print $1+1}'`
ln3=`grep -n "$key" $cmtlst |awk -F: '{print $1+2}'`
ln4=`grep -n "$key" $cmtlst |awk -F: '{print $1+3}'`
cat $cmtlst |sed -n "${ln1}p"
cat $cmtlst |sed -n "${ln2}p"
cat $cmtlst |sed -n "${ln3}p"
cat $cmtlst |sed -n "${ln4}p"
hypo=`cat $cmtlst |sed -n "${ln1}p" |awk '{print substr($0,28,53)}'`
hpcl=`cat $cmtlst |sed -n "${ln1}p" |awk '{print $1}'`
hpdate=`cat $cmtlst |sed -n "${ln1}p" |awk '{print $2}' |awk -F/ '{printf"%4d%3d%3d",$1,$2,$3}'`
hptime=`cat $cmtlst |sed -n "${ln1}p" |awk '{print $3}' |awk -F: '{printf"%3d%3d%4.1f",$1,$2,$3}'`
hplat=`cat $cmtlst |sed -n "${ln1}p" |awk '{print $4}'`
hplon=`cat $cmtlst |sed -n "${ln1}p" |awk '{print $5}'`
hpdep=`cat $cmtlst |sed -n "${ln1}p" |awk '{print $6}'`
eventname=`cat $cmtlst |sed -n "${ln2}p" |awk '{print $1}' |awk -FC '{print $2}'`
halfdur=`cat $cmtlst |sed -n "${ln2}p" |awk '{print substr($0,70,11)}' |awk -F: '{printf"%.4f",$2}'`
ctd_t=`cat $cmtlst |sed -n "${ln3}p" |awk '{printf"%.4f",$2}'`
ctd_lat=`cat $cmtlst |sed -n "${ln3}p" |awk '{printf"%.4f",$4}'`
ctd_lon=`cat $cmtlst |sed -n "${ln3}p" |awk '{printf"%.4f",$6}'`
ctd_dep=`cat $cmtlst |sed -n "${ln3}p" |awk '{printf"%.4f",$8}'`
echo $ctd_t $ctd_lat $ctd_lon $ctd_dep
exp=`cat $cmtlst |sed -n "${ln4}p" |awk '{print $1}'`
mrr=`cat $cmtlst |sed -n "${ln4}p" |awk '{printf"%.6e",$2*10**s}' s=$exp`
mtt=`cat $cmtlst |sed -n "${ln4}p" |awk '{printf"%.6e",$4*10**s}' s=$exp`
mpp=`cat $cmtlst |sed -n "${ln4}p" |awk '{printf"%.6e",$6*10**s}' s=$exp`
mrt=`cat $cmtlst |sed -n "${ln4}p" |awk '{printf"%.6e",$8*10**s}' s=$exp`
mrp=`cat $cmtlst |sed -n "${ln4}p" |awk '{printf"%.6e",$10*10**s}' s=$exp`
mtp=`cat $cmtlst |sed -n "${ln4}p" |awk '{printf"%.6e",$12*10**s}' s=$exp`
cat /dev/null >GCMT_solutions_${type}/CMTSOLUTION_${evtnm}
cat >>GCMT_solutions_${type}/CMTSOLUTION_${evtnm} <<!
 ${hpcl}${hpdate}${hptime}  $hypo
event name:     $eventname
time shift:     $ctd_t
half duration:  $halfdur
latitude:      $ctd_lat
longitude:     $ctd_lon
depth:         $ctd_dep
Mrr:       ${mrr}
Mtt:       ${mtt}
Mpp:       ${mpp}
Mrt:       ${mrt}
Mrp:       ${mrp}
Mtp:       ${mtp}
!
fi
done
fi
exit
#==============================================================
echo "******* Now plot *********"
###@ plot event beach ball @###
ps=final_events_directP.ps

psbasemap -JX15/8 -Rg -Ba45/a30WESN -X2 -Y16 -P -K >$ps
psxy -J -R -Wred -O -K >>$ps <<EOF
>
-122. 35.
-118. 35.
>
-118. 35.
-118. 37.0
>
-118. 37.0
-122. 37.0
>
-122. 37.0
-122. 35.0
EOF
pscoast -J -R -Dc -A10000 -W0.5p,gray -O -K >>$ps 
for cmtf in `ls GCMT_solutions_${type}/CMT*`;do
echo $cmtf
evnm=`cat $cmtf |sed -n '2p' |awk '{printf"%s\n",$3}'`
tshift=`cat $cmtf |sed -n '3p' |awk '{printf"%f\n",$3}'`
hdur=`cat $cmtf |sed -n '4p' |awk '{printf"%f\n",$3}'`
lat=`cat $cmtf |sed -n '5p' |awk '{printf"%f\n",$2}'`
lon=`cat $cmtf |sed -n '6p' |awk '{printf"%f\n",$2}'`
dep=`cat $cmtf |sed -n '7p' |awk '{printf"%f\n",$2}'`
scale=`cat $cmtf |sed -n '8,$p' |awk '{print $2}' |gmtinfo -I1 -C |awk '{if($1*$1>$2*$2) {print $1} else print $2}' |awk -Fe '{print $2}'`
mrr=`cat $cmtf |sed -n '8p' |awk '{printf"%f\n",$2/10**s}' s=$scale`
mtt=`cat $cmtf |sed -n '9p' |awk '{printf"%f\n",$2/10**s}' s=$scale`
mpp=`cat $cmtf |sed -n '10p' |awk '{printf"%f\n",$2/10**s}' s=$scale`
mrt=`cat $cmtf |sed -n '11p' |awk '{printf"%f\n",$2/10**s}' s=$scale`
mrp=`cat $cmtf |sed -n '12p' |awk '{printf"%f\n",$2/10**s}' s=$scale`
mtp=`cat $cmtf |sed -n '13p' |awk '{printf"%f\n",$2/10**s}' s=$scale`
~/progs/sactools_c/distaz 36.0 -120. $lat $lon
psmeca -J -R -Sm0.5 -Gblue -K -O >>$ps <<EOF
$lon $lat $dep $mrr $mtt $mpp $mrt $mrp $mtp $scale X Y
EOF
psxy -J -R -W1p,red -K -O >>$ps <<EOF
-120. 36.0
$lon $lat
EOF
done
# plot line of az=60
project -C-120./36.0 -A60. -G90 -L0/90 |psxy -J -R -W3p,darkgreen,dashed -K -O >>$ps
project -C-120./36.0 -A120. -G90 -L0/90 |psxy -J -R -W3p,darkgreen,dashed -K -O >>$ps
project -C-120./36.0 -A240. -G120 -L0/120 |psxy -J -R -W3p,darkgreen,dashed -K -O >>$ps
project -C-120./36.0 -A300. -G120 -L0/120 |psxy -J -R -W3p,darkgreen,dashed  -O >>$ps

