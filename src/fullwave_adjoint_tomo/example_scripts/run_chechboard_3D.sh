#!/bin/bash
v0=4.0
dv=12.
### SEM simulation domain ####
latmin=3845964.  # Ymin
latmax=4166167.  # Ymax
lonmin=609836.8  # Xmin
lonmax=1012890.4  # Xmax
echo "lonmin,lonmax,latmin,latmax=" $lonmin,$lonmax,$latmin,$latmax
###### center of a profile #####
lonA=`head -1 plots/loc_A_Ap.xyz |awk '{print $1}'` # UTMx
latA=`head -1 plots/loc_A_Ap.xyz |awk '{print $2}'` # UTMy
lonB=`tail -1 plots/loc_A_Ap.xyz |awk '{print $1}'` # UTMx
latB=`tail -1 plots/loc_A_Ap.xyz |awk '{print $2}'` # UTMy
wlon_p=`echo $lonA $lonB |awk '{print $2-$1}'`
wlat_p=`echo $latA $latB |awk '{print $2-$1}'`
clon=`echo $lonA $lonB |awk '{printf"%d",($1+$2)/2 }'`
clat=`echo $latA $latB |awk '{printf"%d",($1+$2)/2 }'`
wlon_l=`echo $lonmin $clon |awk '{print $2-$1 }'`
wlon_r=`echo $lonmax $clon |awk '{print $1-$2 }'`
wlat_l=`echo $latmin $clat |awk '{print $2-$1 }'`
wlat_r=`echo $latmax $clat |awk '{print $1-$2 }'`
echo "profile center:" $clon $clat 
echo "wlon_l,wlon_r,wlat_l,wlat_r=" $wlon_l $wlon_r $wlat_l $wlat_r 
#####
z1=-210000.           # Zmin
z2=-00000.      # Zmax
dlon=20000.
dlat=20000.
dz=20000.

dvlon=20000.
dvlat=20000.
dvz=20000.
####
nlon_l=`echo $wlon_l $dvlon |awk '{print $1/$2}' | awk '{printf("%d\n",$0+=$0<0?0:0.9)}'`
nlon_r=`echo $wlon_r $dvlon |awk '{print $1/$2}' | awk '{printf("%d\n",$0+=$0<0?0:0.9)}'`
nlat_l=`echo $wlat_l $dvlat |awk '{print $1/$2}' | awk '{printf("%d\n",$0+=$0<0?0:0.9)}'`
nlat_r=`echo $wlat_r $dvlat |awk '{print $1/$2}' | awk '{printf("%d\n",$0+=$0<0?0:0.9)}'`
echo "number of grid along x and y:"
echo "nlon_l,nlon_r,nlat_l,nlat_r=" $nlon_l $nlon_r $nlat_l $nlat_r 
lon0=`echo $clon $dvlon $nlon_l |awk '{print $1-$2*($3-1)}'` # away from clon0 by multiple dlon (odd)
lon1=`echo $clon $dvlon $nlon_r |awk '{print $1+$2*($3-1)}'` # away from clon0 by multiple dlon (odd)
lat0=`echo $clat $dvlat $nlat_l |awk '{print $1-$2*($3-1)}'` # away from clat0 by multiple dat (odd)
lat1=`echo $clat $dvlat $nlat_r |awk '{print $1+$2*($3-1)}'` # away from clat0 by multiple dat (odd)
echo "lon0,lon1,lat0,lat1=" $lon0,$lon1,$lat0,$lat1

grid=check_"$dvlon"X"$dvlat".dat
~/progs/FWAT_tools/check_board/check_board_3Dx2<<END
$v0
$dv
$lon0 $lon1 $lat0 $lat1 $z1 $z2
$dlon $dlat $dz
$dvlon $dvlat $dvz
$grid
END
#### Do rotation of anomalies #####
npts=`cat $grid |wc -l`

echo $wlon_p $wlat_p |awk '{printf"scale=3;a(%f)\n",$2/$1}'
rad_angle=`echo $wlon_p $wlat_p |awk '{printf"scale=3;a(%f)\n",$2/$1}' |bc -l` ### atan(wlat/wlon)
angle=`echo $rad_angle |awk '{print $1*180/3.1415926}'`
echo $rad_angle $angle
./rotation_2d $grid $npts $clon $clat $angle >$grid.rt
##################################
cat $grid.rt |awk '{print $1,$2,$3,($4-a)/a}' a=$v0 > checkboard_grid.dat
##############################3
PS=Map_chb_test_3D.ps
xwidth=6
tick=a100000/a100000
CPT=pert.cpt
makecpt -Cpolar -I -T-$dv/$dv/0.5 -D >$CPT 
###
depth=-70000.
echo "plot model slice at depth of " $depth
#RMAP=$lon0/$lon1/$lat0/$lat1
RMAP=$lonmin/$lonmax/$latmin/$latmax
ywidth=`echo $xwidth |awk '{print (y2-y1)/(x2-x1)*$1}' x1=$lonmin x2=$lonmax y1=$latmin y2=$latmax `
PROJ=X${xwidth}/${ywidth}
psbasemap -R$RMAP -J$PROJ -B"$tick"WenN -X4 -Y18 -K -P > $PS
cat $grid.rt |awk -v d=$depth '{if($3==d) print $1,$2,($4-a)/a*100}' a=$v0 > slice.xyz
#blockmean slice.xyz -R$RMAP -I$dlon/$dlat >slice.xyz.m 
#surface slice.xyz.m -R$RMAP -Gslice.grd  -I$dlon/$dlat
xyz2grd slice.xyz -R$RMAP -Gslice.grd  -I$dlon=/$dlat=
grdimage slice.grd -J$PROJ -R$RMAP -C$CPT -B"$tick"::wens -O -K >> $PS
cat slice.xyz |psxy -J -R -Sc0.05i -G255/255/255 -O -K >>$PS
psxy -J -R -W1p -O -K >>$PS <<EOF
$lon0 $lat0
$lon1 $lat0
$lon1 $lat1
$lon0 $lat1
$lon0 $lat0
EOF
cat plots/loc_A_Ap.xyz |psxy -J -R -Sa0.1i -Gyellow -O -K >>$PS
cat plots/loc_A_Ap.xyz |psxy -J -R -W1p -O -K >>$PS
###
y=$clat
echo "plot vertical cross-section at y= " $y
RMAP=$lonmin/$lonmax/-220000/0
ywidth=`echo $xwidth |awk '{print (y2-y1)/(x2-x1)*$1}' x1=$lonmin x2=$lonmax y1=-220000 y2=0. `
PROJ=X${xwidth}/${ywidth}

psbasemap -R$RMAP -J$PROJ -B"$tick"WenN -O -K -X10 >> $PS
cat $grid |awk -v a=$y '{if($2==a) print $1,$3,($4-b)/b*100}' b=$v0 > slice.xyz
#blockmean slice.xyz -R$RMAP -I$dlon/$dlat >slice.xyz.m 
#surface slice.xyz.m -R$RMAP -Gslice.grd  -I$dlon/$dlat
xyz2grd slice.xyz -R$RMAP -Gslice.grd  -I$dlon=/$dz=
grdimage slice.grd -J$PROJ -R$RMAP -C$CPT -B"$tick"::wens -O -K >> $PS
cat slice.xyz |psxy -J -R -Sc0.05i -G255/255/255 -O -K >>$PS
psscale -D4./7.5/6.0/0.5h -C$CPT -Ba5 -E -O -N >> $PS
############
