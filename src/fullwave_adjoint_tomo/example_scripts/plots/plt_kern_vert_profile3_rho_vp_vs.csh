#!/bin/csh
module load gcc
gmtset BASEMAP_TYPE plain
gmtset MEASURE_UNIT cm
gmtset ANOT_FONT_SIZE 12p LABEL_FONT_SIZE 12p
gmtset LABEL_FONT Helvetica ANNOT_FONT_PRIMARY Helvetica ANNOT_FONT_SECONDARY Helvetica HEADER_FONT Helvetica
gmtset LABEL_OFFSET 0.1c
gmtset COLOR_BACKGROUND 28/48/104
gmtset COLOR_FOREGROUND white
set only_plot = 0
set ie = 1
while ( $ie <= 1 ) # loop over all virtual evts
  set evtfile=`cat profiles.dat |sed -n "${ie}p"` # loop over all virtual evts
  set eid = `echo $evtfile |awk '{printf"%s",$5}'`
  set elon = `echo $evtfile |awk '{print $1}'`
  set elat = `echo $evtfile |awk '{print $2}'`
  set stnm = `echo $evtfile |awk '{printf"%s'\''",$5}'`
  set stnmp = `echo $evtfile |awk '{printf"%sp",$5}'`
  set slon = `echo $evtfile |awk '{print $3}'`
  set slat = `echo $evtfile |awk '{print $4}'`
  echo "$eid $elon/$elat $stnm $slon/$slat "
## 
#set dir = model_input_smooth
#set iter = input
#set md = $dir/rho.xyz
#set md1 = $dir/vp.xyz
#set md2 = $dir/vs.xyz

set mod = M00.set1.5
set dir = .
set md = $dir/grad.rhop_kernel.${mod}.regridded.xyz
set md1 = $dir/grad.alpha_kernel.${mod}.regridded.xyz
set md2 = $dir/grad.beta_kernel.${mod}.regridded.xyz

#set dir = model_output
#set iter = output
#set md = $dir/000001_model_rh_output.xyz
#set md1 = $dir/000001_model_vp_output.xyz
#set md2 = $dir/000001_model_vs_output.xyz

set cmax = 4.2
set cmin = 2.8
set swathwidth = 1.0
set slicewidth = 1.6 ## y slice width
set xtick=a100f50
set ytick=a100f50
#-----------------------
set iflagdeg = 2 # iflagdeg = 1 plot degree vs. depth else plot distance vs. depth
set PS = fig_kern_profile_${mod}_${eid}_${stnmp}.ps
set PROG1 = preplot_ray_map100_nan
# create color table for dvp
set CPT = rho.cpt
set CPT1 = vp.cpt
set CPT2 = vs.cpt
#$PDIR/svcpt <<!>> junk
#$cmin $cmax
#!
#mv svel13.cpt $CPT
makecpt -Csvel13.cpt -T-2/2.0/0.02 -D >$CPT
makecpt -Csvel13.cpt -T-2/2.0/0.02 -D >$CPT1
makecpt -Csvel13.cpt -T-2/2.0/0.02 -D >$CPT2
#makecpt -Cdrywet -T2.6/3.6/0.02 -D >$CPT
#makecpt -Cdrywet -T5.0/8.9/0.02 -D >$CPT1
#makecpt -Cdrywet -T3.0/5.0/0.02 -D >$CPT2
#makecpt -Cseis -T2.6/3.5/0.02 -D >$CPT
#makecpt -Cseis -T5.6/8.6/0.02 -D >$CPT1
#makecpt -Cseis -T3.3/4.8/0.02 -D >$CPT2


#makecpt -Cseis -T-0.5/0.5/0.01 -D >$CPT
#--------------------------------------------
set RMAP =  -122.2/-116.8/34.75/38.25
#set GCP1 = ( -120/32.2 -119/32.2 -118/32.2 -117/32.2    )# -125/39.5  -125/44.1
#set GCP2 = ( -120/36.8 -119/36.8 -118/36.8 -117/36.8     )#-100/39.5  -100/44.1
set GCP1 = ( $elon/$elat    )# -125/39.5  -125/44.1
set GCP2 = ( $slon/$slat     )#-100/39.5  -100/44.1
set plane = ( yzpqrs xzpqrs xzpqrs xzpqrs yzpqrs )
set t1 = ( ${eid}  "B"  "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" )
set t2 = ( ${stnm} "B'" "C'" "D'" "E'" "F'" "G'" "H'" "I'" "J'" "K'" "L'" )
###------- make image.3D with fourth cloumns: x, y, z, radius
if ( $only_plot == 0 ) then
./convert_utm2lonlat.pl $md 10 >model.regridded.lonlat
# M28.vsh
rm -f image.3D
cat $md |awk '{print $4}' > temp
cat model.regridded.lonlat |awk '{print $1,$2,6371+$3/1000}' >temp1
paste temp1 temp >image.3D
#M28.vsv
rm -f image.3D.1
cat $md1 |awk '{print $4}' > temp
paste temp1 temp >image.3D.1
# M28.voigt_vs
rm -f image.3D.2
cat $md2 |awk '{print $4}' > temp
paste temp1 temp >image.3D.2
date

date
echo "extract 3D data finished!"
endif
###----------------------
set FV = image.proj
set kmdeg = `echo 6371 | awk '{print 2*3.1416*$1/360.}'`
set DLON2km = `echo 70  | awk '{print $1*'$kmdeg' }'`
@ i = 0
foreach SL ( $GCP1 )
@ i ++
set GC1 = $SL
set GC2 = $GCP2[$i]

###-----------
\rm -f $FV $FV.1 $FV.2 #$FV.3
date
### method 1: shung-huei's method for regular mesh 
#project image.3D -C$GC1 -E$GC2 -W-$DLON2km/$DLON2km -F$plane[$i] -Q -Lw > $FV
#project image.3D -C$GC1 -E$GC2 -F$plane[$oi] -Q -Lw > $FV
###  method 2: seisman
if ( $only_plot == 0 ) then
set scale = 11
project image.3D -C$GC1 -E$GC2 -G0.02 -Q -Lw >track.dat
foreach dep  ( `seq 0 2 200` )
  echo $dep
  set rad = `echo $dep |awk '{print 6371-$1}'`;
  # M16
  cat image.3D |awk -v a=$rad '{if(sqrt(($3-a)**2)<1) print $1,$2,$4*10**s}' s=$scale >dep.xyd
  #xyz2grd dep.xyd -R$RMAP -I0.02 -Gdep.grd 
  surface dep.xyd -R$RMAP -I0.02 -Gdep.grd 
  grdtrack track.dat -Gdep.grd |awk -v d=$dep '{print $2,6371-d,$4,$3}'>>$FV
  # M21
  cat image.3D.1 |awk -v a=$rad '{if(sqrt(($3-a)**2)<1) print $1,$2,$4*10**s}' s=$scale >dep.xyd
  #xyz2grd dep.xyd -R$RMAP -I0.02 -Gdep.grd 
  surface dep.xyd -R$RMAP -I0.02 -Gdep.grd 
  grdtrack track.dat -Gdep.grd |awk -v d=$dep '{print $2,6371-d,$4,$3}'>>$FV.1
  # Barak2015
  cat image.3D.2 |awk -v a=$rad '{if(sqrt(($3-a)**2)<1) print $1,$2,$4*10**s}' s=$scale >dep.xyd
  #xyz2grd dep.xyd -R$RMAP -I0.02 -Gdep.grd 
  surface dep.xyd -R$RMAP -I0.02 -Gdep.grd 
  grdtrack track.dat -Gdep.grd |awk -v d=$dep '{print $2,6371-d,$4,$3}'>>$FV.2
  # CVM-S.26
  #cat image.3D.3 |awk -v a=$dep '{if($3==(6371-a)) print $1,$2,$4}' >dep.xyd
  #xyz2grd dep.xyd -R$RMAP -I0.02 -Gdep.grd 
  #surface dep.xyd -R$RMAP -I0.02 -Gdep.grd 
  #grdtrack track.dat -Gdep.grd |awk -v d=$dep '{print $2,6371-d,$4,$3}'>>$FV.3
end
echo "extract 2D profiles finished"
endif
date
### moho
 #project moho_lupei_zhu.dat -C$GC1 -E$GC2 -G0.02 -Q -Lw >track.dat
 #xyz2grd moho_lupei_zhu.dat -R$RMAP -I0.1 -Gmoho.grd
 #grdtrack track.dat -Gmoho.grd |awk  '{print $3,$4}' >moho.dat
 # cp moho.dat moho.clip
 # tail -1 moho.dat |awk '{print $1,60.0}' >>moho.clip
 # head -1 moho.dat |awk '{print $1,60.0}' >>moho.clip
 # head -1 moho.dat  >>moho.clip
### rename and store 3D data
if ( $only_plot == 0 ) then
# M16
echo $FV
cat $FV | awk '{print $1, 6371-$2, $3 }' > $FV.deg
cat $FV | awk '{print $4, 6371-$2, $3 }' > $FV.km
mv $FV $FV.pqrs

if ( $iflagdeg == 1 ) then
   cp $FV.deg $FV
else
   cp $FV.km ${eid}_${stnm}.$FV
endif
# M21
echo $FV.1
cat $FV.1 | awk '{print $1, 6371-$2, $3 }' > $FV.deg.1
cat $FV.1 | awk '{print $4, 6371-$2, $3 }' > $FV.km.1
mv $FV.1 $FV.pqrs.1

if ( $iflagdeg == 1 ) then
   cp $FV.deg.1 $FV.1
else
   cp $FV.km.1 ${eid}_${stnm}.$FV.1
endif
# Barak2015
echo $FV.2
cat $FV.2 | awk '{print $1, 6371-$2, $3 }' > $FV.deg.2
cat $FV.2 | awk '{print $4, 6371-$2, $3 }' > $FV.km.2
mv $FV.2 $FV.pqrs.2

if ( $iflagdeg == 1 ) then
   cp $FV.deg.2 $FV.2
else
   cp $FV.km.2 ${eid}_${stnm}.$FV.2
endif
# CVM-S.26
#echo $FV.3
#cat $FV.3 | awk '{print $1, 6371-$2, $3 }' > $FV.deg.3
#cat $FV.3 | awk '{print $4, 6371-$2, $3 }' > $FV.km.3
#mv $FV.3 $FV.pqrs.3
#
#if ( $iflagdeg == 1 ) then
#   cp $FV.deg.3 $FV.3
#else
#   cp $FV.km.3 ${eid}_${stnm}.$FV.3
#endif
endif
###------
set FV = ${eid}_${stnm}.$FV
#set range = `minmax -C $FV | awk ' { print $1 "/" $2 "/" $3 "/" $4 }'`
set range = `minmax -C $FV | awk ' { print $1 "/" $2 "/" "0/200" }'`
#set range = `minmax -C $FV | awk ' { print "0/220/0/100" }'`
#set range = `minmax -C $FV | awk ' {if($2<=200) {print $1 "/" $2 "/0/50"} else {print "0/400/-4/60"} }'`
minmax -C $FV | awk ' { print $1 "/" $2 "/" $3 "/" $4 }'
#if ( $slicewidth != 0 ) then
echo "range=" $range
#if ( $i == 1 ) then
set scaleunit = `echo $range | awk -F/ '{print '$slicewidth'/($4-$3) } '`
echo $scaleunit
#set xwidth = `echo $slicewidth | awk '{print $1-1}'`
set ywidth = `echo $range | awk -F/ '{print ($4-$3)*'$scaleunit' } '`
#else
set xwidth = `echo $range | awk -F/ '{print ($2-$1)*'$scaleunit' }'` 
#endif
#else
#set xwidth = 3.0
#set ywidth = 1.5
#endif
echo $i $xwidth $ywidth
### interpolate
#triangulate $FV -R$range -I$DD/$DDEP -Gslice.grd -E > tria.out
#triangulate $FV -R$range -I10/10 -Gslice.grd -E > tria.out
set xint = 2.0
set zint = 1.0
blockmean $FV -R$range -I$xint/$zint > $FV.bm
surface $FV.bm -R$range -I$xint/$zint -Gslice.grd 
#xyz2grd $FV -R$range  -I2/5.12 -Gslice.grd

blockmean $FV.1 -R$range  -I$xint/$zint > $FV.bm.1
surface $FV.bm.1 -R$range -I$xint/$zint -Gslice.grd.1 
#xyz2grd $FV.1 -R$range  -I2/5.12 -Gslice.grd.1

blockmean $FV.2 -R$range -I$xint/$zint > $FV.bm.2
surface $FV.bm.2 -R$range -I$xint/$zint -Gslice.grd.2
#xyz2grd $FV.2 -R$range  -I2/5.12 -Gslice.grd.2

#blockmean $FV.3 -R$range  -I$xint/$zint > $FV.bm.3
#surface $FV.bm.3 -R$range -I$xint/$zint -Gslice.grd.3
#xyz2grd $FV.3 -R$range  -I2/1 -Gslice.grd.3
#set dist2 = `echo $range | awk -F/ ' { print ($1+$2)*0.5 } '`
echo "plot"
set sizearrow = 0.2i
set scalarrow = 0.03i/0.08i/0.04i
set SCALE = X${xwidth}i/-${ywidth}i
#M16
psbasemap -R$range -J$SCALE -B"$xtick"/"$ytick":"Depth (km)":Wesn -K -X8 -Y20 -P > $PS
#psclip moho.clip -J -R -N -O -K >>$PS
grdimage slice.grd -R -J -B"$xtick"/"$ytick"::wens -C$CPT -O -K -P >> $PS
cat RF_west_bound.csv |sed -n '7,$p' |awk -F, '{print $1,-$2}' |psxy -J -R -N -W1p -O -K >>$PS
cat RF_east_bound1.csv |sed -n '7,$p' |awk -F, '{print $1,-$2}' |psxy -J -R -N -W1p -O -K >>$PS
cat RF_east_bound2.csv |sed -n '7,$p' |awk -F, '{print $1,-$2}' |psxy -J -R -N -W1p -O -K >>$PS
#grdcontour slice.grd -W0.75p,black -C0.05 -L3.5/4.7 -A+g+o+p1+s8 -Gn1+r0.75c -S -R -J -O -K >> $PS
#grdcontour slice.grd.3 -W1p,black -C0.05 -L3.0/4.7 -A+g+o+f1+p1+s8+N -S -R -J -O -K >> $PS
#psclip -J -R -C -O -K >>$PS
psscale -D3.2i/0.8i/1.3i/0.5 -C$CPT -Ba0.4f0.2/:"Rho (g/cm@+3@+)": -V -E -O -K -N >> $PS
#grdcontour slice.grd -J -R -O -K -Ccont.d -W1/255 >> $PS 
### plot moho
#cat moho.dat |psxy -J -R -W1p -O -K >> $PS

set tx1 = `echo $range | awk -F/ ' { print $1 } '`
set tx2 = `echo $range | awk -F/ ' { print $2 } '`
set tx = `echo $range | awk -F/ ' { print ($1+$2)/2 } '`
pstext -R -J -O -K -N -P <<EOF>>$PS
$tx1 -20 12 0 2 MC $t1[$i]
$tx2 -20 12 0 2 MC $t2[$i]
EOF
#echo $tx1 $tx2 $i $t1[$i] $t2[$i]
#M21
psbasemap -R$range -J$SCALE -B"$xtick"/"$ytick":"Depth (km)":Wens -K -Y-1.9i -O  >> $PS  
#psclip moho.clip -J -R -N -O -K >>$PS
grdimage slice.grd.1 -R -J -B"$xtick"/"$ytick"::wens -C$CPT1 -O -K -P >> $PS
cat RF_west_bound.csv |sed -n '7,$p' |awk -F, '{print $1,-$2}' |psxy -J -R -N -W1p -O -K >>$PS
cat RF_east_bound1.csv |sed -n '7,$p' |awk -F, '{print $1,-$2}' |psxy -J -R -N -W1p -O -K >>$PS
cat RF_east_bound2.csv |sed -n '7,$p' |awk -F, '{print $1,-$2}' |psxy -J -R -N -W1p -O -K >>$PS

#grdcontour slice.grd.1 -W0.75p,black -C0.05 -L3.5/4.7 -A+g+o+p1+s8 -Gn1+r0.75c -S -R -J -O -K >> $PS
#grdcontour slice.grd.3 -W1p,black -C0.05 -L3.0/4.7 -A+g+o+f1+p1+s8+N -S -R -J -O -K >> $PS
#psclip -J -R -C -O -K >>$PS
psscale -D3.2i/0.8i/1.3i/0.5 -C$CPT1 -Ba0.5/:"Vp (km/s)": -V -E -O -K -N >> $PS
#grdcontour slice.grd -J -R -O -K -Ccont.d -W1/255 >> $PS 
### plot moho
#cat moho.dat |psxy -J -R -W1p -O -K >> $PS
#pstext -R -J -O -K -N -P <<EOF>>$PS
#$tx1 -10 12 0 2 MC $t1[$i]
#$tx2 -10 12 0 2 MC $t2[$i]
#EOF

#Barak2015
psbasemap -R$range -J$SCALE -B"$xtick":"Distance (km)":/"$ytick":"Depth (km)":WenS -K -Y-1.9i -O  >> $PS
#psclip moho.clip -J -R -N -O -K >>$PS
grdimage slice.grd.2 -R -J -B -C$CPT2 -O -K -P >> $PS
cat RF_west_bound.csv |sed -n '7,$p' |awk -F, '{print $1,-$2}' |psxy -J -R -N -W1p -O -K >>$PS
cat RF_east_bound1.csv |sed -n '7,$p' |awk -F, '{print $1,-$2}' |psxy -J -R -N -W1p -O -K >>$PS
cat RF_east_bound2.csv |sed -n '7,$p' |awk -F, '{print $1,-$2}' |psxy -J -R -N -W1p -O -K >>$PS

#psclip -J -R -C -O -K >>$PS
psscale -D3.2i/0.8i/1.3i/0.5 -C$CPT2 -Ba0.4f0.2/:"Vs (km/s)": -V -E -O -K -N >> $PS
#grdcontour slice.grd -J -R -O -K -Ccont.d -W1/255 >> $PS 
### plot moho
#cat moho.dat |psxy -J -R -W1p -O -K >> $PS
#pstext -R -J -O -K -N -P <<EOF>>$PS
#$tx1 -10 12 0 2 MC $t1[$i]
#$tx2 -10 12 0 2 MC $t2[$i]
#EOF



end # end of foreach slice

cat /dev/null |psxy -J -R -Sc0.1i -O >>$PS
  @ ie++
end # end loop of sources



