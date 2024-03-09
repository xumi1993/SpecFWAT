#!/bin/bash
cmtlst=GCMT_event_mag5.8_all.lst
rm -rf GCMT_solutions
mkdir -p GCMT_solutions

if true;then
###@ get CMTSOLUTION for final events @###
cat final_events.lst |
while read line;do
evtnm=`echo $line |awk '{print $1}' |awk -F/ '{print $3}'`
y=`echo $evtnm |awk -F_ '{print $2}'`
m=`echo $evtnm |awk -F_ '{print $3}'`
d=`echo $evtnm |awk -F_ '{print $4}'`
min=`echo $evtnm |awk -F_ '{print $5}'`
sec=`echo $evtnm |awk -F_ '{print $6}'`
msec=`echo $evtnm |awk -F_ '{print $7}'`
echo $evtnm
key=`echo $y $m $d $min $sec $msec |awk '{printf" PDEW%4d%3d%3d%3d%3d%3d\n",$1,$2,$3,$4,$5,$6}'`
grep -n "$key" $cmtlst
ln1=`grep -n "$key" $cmtlst |awk -F: '{print $1}'`
ln11=`grep -n "$key" $cmtlst |awk -F: '{print $1+1}'`
ln2=`echo $ln1 |awk '{print $1+12}'`
#evnm=`cat $cmtlst |sed -n "${ln11}p" |awk '{printf"%s\n",$3}'`
cat $cmtlst |sed -n "${ln1},${ln2}p" >GCMT_solutions/CMTSOLUTION_${evtnm}

done
fi

###@ plot event beach ball @###
ps=final_events.ps

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
pscoast -J -R -Dc -A10000 -W0.5p,black -O -K >>$ps 
for cmtf in `ls GCMT_solutions/CMT*`;do
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

