#!/bin/bash
module load gcc/7.3.0 hdf5/1.8.20 netcdf/4.6.1

gmt gmtset FONT_TITLE 24p
gmt gmtset FONT_ANNOT_PRIMARY 14p
gmt gmtset FONT_LABEL 16p
gmt gmtset MAP_FRAME_TYPE plain

indir=.
infile=model.vs_smooth.M33.regridded.xyz

PS=$indir/${infile}.ps
width=2.2
PROJ=M${width}i
tick=a2f1/a1
lonmin=-122.0
lonmax=-117.0
latmin=34.5
latmax=37.5
RMAP=$lonmin/$lonmax/$latmin/$latmax

# create color table for dvp
CPT=dvp.cpt
cmax=13
makecpt -Csvel13.cpt -T-$cmax/$cmax/2 -D >$CPT


if true;then
./convert_utm2lonlat.pl $infile 10 >model.regridded.lonlat
cat $infile |awk '{print $4}' > temp
cat model.regridded.lonlat |awk '{print $1,$2,$3}' >temp1
paste temp1 temp >model.3D
echo "finish convert_utm2lonlat"
fi

i=1
#for dep in 10 40 70 100 130 160 190 210;do
for dep in 15 25 45 60 75 90 110 130 ;do
#for dep in 20 40 60 80 100 120 140 160;do
    echo $i
    z=`echo $dep |awk '{print -$1*1000}'`
    cat $indir/model.3D |awk '{if($3==a&&$4!=-1000) print $1,$2,$4}' a=$z >$indir/dep$dep.xyz

    ref_vs=`cat $indir/dep$dep.xyz |awk 'BEGIN{sum=0.0;num=0} {{sum=sum+$3;num=num+1}} END{print sum/num}' `
    echo $dep $ref_vs
    # for vs 
    cat $indir/dep$dep.xyz |awk '{print $1,$2,($3-a)/a*100}' a=$ref_vs >slow.slice
    # for RA of M28
    #cat $indir/dep$dep.xyz |awk '{print $1,$2,$3}' a=$ref_vs >slow.slice
    cp slow.slice slow.slice.$dep
    blockmedian slow.slice -R$RMAP -I0.02 >slow.slice.m
    surface slow.slice.m -R$RMAP -I0.02 -Gslice.grd0
    #xyz2grd slow.slice -R$RMAP -I0.02 -Gslice.grd
    #triangulate slow.slice -R$RMAP -I0.25/0.25 -Gslice.grd -E > tria.out
    grdfilter slice.grd0 -D1 -Fg100 -Gslice.grd

    if [ $i -le 4 ];then
    if [ $i -eq 1 ];then
    gmt psbasemap -R$RMAP -J$PROJ -B"$tick"::WenN -X4 -Y22 -K -P > $PS
    elif [ $i -eq 4 ];then
    gmt psbasemap -R$RMAP -J$PROJ -B"$tick"::WSen -Y-5.0 -O -K >> $PS
    else
    gmt psbasemap -R$RMAP -J$PROJ -B"$tick"::Wesn -Y-5.0 -O -K >> $PS
    fi
    else
    if [ $i -eq 5 ];then
    gmt psbasemap -R$RMAP -J$PROJ -B"$tick"::wEsN -X6.5 -Y15.0 -O -K >> $PS
    elif [ $i -eq 8 ];then
    gmt psbasemap -R$RMAP -J$PROJ -B"$tick"::wEnS -Y-5.0 -O -K >> $PS
    else
    gmt psbasemap -R$RMAP -J$PROJ -B"$tick"::wEns -Y-5.0 -O -K >> $PS
    fi
    fi
    gmt grdimage slice.grd -J$PROJ -R$RMAP -C$CPT -B"$tick"::wens -O -K >> $PS

    vstick=`echo $cmax |awk '{print $1/4}'`
    gmt pscoast -J -R  -O -K -W1p,black -Na/1p,black -Dh >> $PS
    gmt psxy jennings.xy -J -R -W0.4p,0/0/0 -O -K >> $PS
gmt pstext -R -J -W0.5p -G255/255/255 -O -K -N <<EOF >> $PS
-121.7 37.2 12 0 1 ML $dep km
EOF
gmt psxy -J -R  -W1p -O -K <<EOF >>$PS
-121.56 35.93 
-117.5 36.63
EOF
gmt psxy Cal_geo_provins/geoprov.dat -J -R -W1p,255/0/0,10_10:- -O -K >> $PS
cat station.lst |awk '{print $4,$3}' |gmt psxy -J -R -St0.04i -W0.1p,black -O -K >>$PS
gmt pstext -J -R -F+f10,1,black+a+j  -Gblack -To -t60 -O -K >>$PS <<eof
-118.55 36.0  0  6 SNB
-119.8 35.25  -45 6 SAF
-116.3 33.8 -45 6 SAF
-119.57 35.75 -60 6 GV
-120.75 35.75 -60 6 SCR
-117.3 36.3 -60 6 WL
eof
gmt pstext -J -R -F+f10,1,white+a+j  -O -K >>$PS <<eof
-118.55 36.0  0  6 SNB
-119.8 35.25  -45 6 SAF
-116.3 33.8 -45 6 SAF
-119.57 35.75 -60 6 GV
-120.75 35.75 -60 6 SCR
-117.3 36.3 -60 6 WL
eof

    let i=i+1
done 
gmt psscale -Dx-6.5/-2.+w4i/0.3+h+e -C$CPT -B4::/:"dln vs (%) ": -O -N >> $PS
#psscale -D-0.25/-2/4i/0.3h -C$CPT -B"$vstick"::/:"vs RA (%) ": -E -O -N >> $PS
exit

