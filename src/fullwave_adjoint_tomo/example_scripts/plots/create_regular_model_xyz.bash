#!/bin/sh
# AK135 form TAUP_2.4.5
# creates a simple, formatted tomography model with a constant velocity gradient
# for a block model with dimensions 134000 x 134000 x 60000

# origin points
ORIG_X=`grep "UTM X min" ../OUTPUT_FILES/output_meshfem3D.txt | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 |awk '{printf"%f\n",$1}'`
ORIG_Y=`grep "UTM Y min" ../OUTPUT_FILES/output_meshfem3D.txt | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 |awk '{printf"%f\n",$1}'`
ORIG_Z=0.

# end points
END_X=`grep "UTM X max" ../OUTPUT_FILES/output_meshfem3D.txt | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 |awk '{printf"%f\n",$1}'`
END_Y=`grep "UTM Y max" ../OUTPUT_FILES/output_meshfem3D.txt | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 |awk '{printf"%f\n",$1}'`
END_Z=-210000. #`grep ^DEPTH_BLOCK_KM ../DATA/meshfem3D_files/Mesh_Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 |awk '{printf"%f\n",-$1*1000.}'`

# spacing of given tomography points
SPACING_X=2500.
SPACING_Y=2500.
SPACING_Z=-2500.

# number of cell increments
NX=`echo $ORIG_X $END_X  $SPACING_X |awk '{printf"%d\n",($2-$1)/$3+1}'`
NY=`echo $ORIG_Y $END_Y  $SPACING_Y |awk '{printf"%d\n",($2-$1)/$3+1}'`
NZ=`echo $ORIG_Z $END_Z  $SPACING_Z |awk '{printf"%d\n",($2-$1)/$3+1}'`

# min/max values
VP_MIN=5800.
VP_MAX=8500.
VS_MIN=3460.
VS_MAX=4609.
RHO_MIN=2720.
RHO_MAX=3457.


# header info
echo "creating header info..."

echo "$ORIG_X $ORIG_Y $ORIG_Z $END_X $END_Y $END_Z  " #> tmp.xyz
echo "$SPACING_X $SPACING_Y $SPACING_Z " #>> tmp.xyz
echo "$NX $NY $NZ " #>> tmp.xyz
echo "$VP_MIN $VP_MAX $VS_MIN $VS_MAX $RHO_MIN $RHO_MAX" #>> tmp.xyz


# adds point location and velocity model values
echo "adding model values..."

# format: lists first all x, then y, then z
echo "1" | awk '{ for(k=0;k<NZ;k++){ for(j=0;j<NY;j++){for(i=0;i<NX;i++){ \
x=ORIG_X+i*SPACING_X;y=ORIG_Y+j*SPACING_Y;z=ORIG_Z+k*SPACING_Z; \
if(z<-210000)      {vp=8300+(z+210000)*(8482.5-8300)/(-260000+210000);vs=4523+(z+210000)*(4609-4523)/(-260000+210000);rho=3425.8+(z+210000)*(3456.1-3425.8)/(-260000+210000);} \
else if(z<-165000) {vp=8175+(z+165000)*(8300-8175)/(-210000+165000);vs=4509+(z+165000)*(4518-4509)/(-210000+165000);rho=3398.5+(z+165000)*(3425.8-3398.5)/(-210000+165000);} \
else if(z<-120000) {vp=8050+(z+120000)*(8175-8050)/(-165000+120000);vs=4500+(z+120000)*(4509-4500)/(-165000+120000);rho=3371.3+(z+120000)*(3398.5-3371.3)/(-165000+120000);} \
else if(z<-77500)  {vp=8045+(z+77500)*(8050-8045)/(-120000+77500);vs=4490+(z+77500)*(4500-4490)/(-120000+77500);rho=3345.5+(z+77500)*(3371.3-3345.5)/(-120000+77500);} \
else if(z<-35000)  {vp=8040+(z+35000)*(8045-8040)/(-77500+35000);vs=4480+(z+35000)*(4490-4480)/(-77500+35000);rho=3319.8+(z+35000)*(3345.5-3319.8)/(-77500+35000);} \
else if(z<-20000)  {vp=6500;vs=3850;rho=2920;} \
else               {vp=5800;vs=3460;rho=2720;} \
print x,y,z,vp,vs,rho }}} }' \
NX=$NX NY=$NY NZ=$NZ SPACING_X=$SPACING_X SPACING_Y=$SPACING_Y SPACING_Z=$SPACING_Z ORIG_X=$ORIG_X ORIG_Y=$ORIG_Y ORIG_Z=$ORIG_Z > model_regrid.xyz

echo "created file: tomography_model.xyz"
