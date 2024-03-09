#!/bin/bash
gmtset MEASURE_UNIT cm
gmtset HEADER_FONT_SIZE 24p
gmtset BASEMAP_TYPE plain 
indir=.
infile=model.vs.M22.regridded.xyz

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
makecpt -Csvel13m.cpt -T-$cmax/$cmax/2 -D >$CPT


if false;then
./convert_utm2lonlat.pl $infile 10 >model.regridded.lonlat
cat $infile |awk '{print $4}' > temp
cat model.regridded.lonlat |awk '{print $1,$2,$3}' >temp1
paste temp1 temp >model.3D
echo "finish convert_utm2lonlat"
fi

i=1
#for dep in 10 40 70 100 130 160 190 210;do
#for dep in 10 30 50 70 90 110 130 150;do
for dep in 20 40 60 80 100 120 140 160;do
    echo $i
    z=`echo $dep |awk '{print -$1*1000}'`
    cat $indir/model.3D |awk '{if($3==a&&$4!=-1000) print $1,$2,$4}' a=$z >$indir/dep$dep.xyz

    ref_vs=`cat $indir/dep$dep.xyz |awk 'BEGIN{sum=0.0;num=0} {{sum=sum+$3;num=num+1}} END{print sum/num}' `
    echo $dep $ref_vs
    # for vs 
    cat $indir/dep$dep.xyz |awk '{print $1,$2,($3-a)/a*100}' a=$ref_vs >slow.slice
    # for RA of M28
    #cat $indir/dep$dep.xyz |awk '{print $1,$2,$3*100}' a=$ref_vs >slow.slice
    cp slow.slice slow.slice.$dep
    blockmedian slow.slice -R$RMAP -I0.02 >slow.slice.m
    surface slow.slice.m -R$RMAP -I0.02 -Gslice.grd
    #xyz2grd slow.slice -R$RMAP -I0.02 -Gslice.grd
    #triangulate slow.slice -R$RMAP -I0.25/0.25 -Gslice.grd -E > tria.out

    if [ $i -le 4 ];then
    if [ $i -eq 1 ];then
    psbasemap -R$RMAP -J$PROJ -B"$tick"::WenN -X4 -Y22 -K -P > $PS
    elif [ $i -eq 4 ];then
    psbasemap -R$RMAP -J$PROJ -B"$tick"::WSen -Y-5.0 -O -K >> $PS
    else
    psbasemap -R$RMAP -J$PROJ -B"$tick"::Wesn -Y-5.0 -O -K >> $PS
    fi
    else
    if [ $i -eq 5 ];then
    psbasemap -R$RMAP -J$PROJ -B"$tick"::wEsN -X6.5 -Y15.0 -O -K >> $PS
    elif [ $i -eq 8 ];then
    psbasemap -R$RMAP -J$PROJ -B"$tick"::wEnS -Y-5.0 -O -K >> $PS
    else
    psbasemap -R$RMAP -J$PROJ -B"$tick"::wEns -Y-5.0 -O -K >> $PS
    fi
    fi
    grdimage slice.grd -J$PROJ -R$RMAP -C$CPT -B"$tick"::wens -O -K >> $PS

    vstick=`echo $cmax |awk '{print $1/4}'`
    pscoast -J -R  -O -K -W1p,black -Na/1p,black -Dh >> $PS
pstext -R -J -W255/255/255,o,0.5p -O -K -N <<EOF >> $PS
-121.7 37.2 12 0 1 ML $dep km
EOF
psxy -J -R  -W1p -O -K <<EOF >>$PS
-121.56 35.93 
-117.5 36.63
EOF
psxy Cal_geo_provins/geoprov.dat -J -R -W1p,255/0/0,10_10:- -O -K >> $PS

    let i=i+1
done 
psscale -D-0.25/-2/4i/0.3h -C$CPT -B4::/:"dln vs (%) ": -E -O -N >> $PS
#psscale -D-0.25/-2/4i/0.3h -C$CPT -B"$vstick"::/:"vs RA (%) ": -E -O -N >> $PS
exit

