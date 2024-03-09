#!/usr/bin/env python3

import sys
import os
import numpy as np
from io import StringIO
from pyproj import Proj

if len(sys.argv) != 4:
  print("Usage: correct_wavefront_time.py FKmodel zbottom FKfile")
  exit()
else:
  print(str(sys.argv))

fkmodel=sys.argv[1]
zbottom=float(sys.argv[2])
fkfile=sys.argv[3]

###### read FKmodel ##########
f=open(fkmodel,'r')
print(fkmodel)
for line in f:
  if len(line)<0 or line[0]=='#':
    continue
  data=line.split()
  print(data)
  if data[0]=='NLAYER':
    nlayer=int(data[1])
    ilayer_fk_input=np.zeros(nlayer,dtype=int)
    rho_fk_input=np.zeros(nlayer)
    vp_fk_input=np.zeros(nlayer)
    vs_fk_input=np.zeros(nlayer)
    ztop_fk_input=np.zeros(nlayer)
  if data[0]=='LAYER':
    ilayer=int(data[1])
    ilayer_fk_input[ilayer-1]=int(data[1])
    rho_fk_input[ilayer-1]=float(data[2])
    vp_fk_input[ilayer-1]=float(data[3])
    vs_fk_input[ilayer-1]=float(data[4])
    ztop_fk_input[ilayer-1]=float(data[5])
  if data[0]=='INCIDENT_WAVE':
    if data[1]=='p':
      kpsv=1
    elif data[1]=='s':
      kpsv=2
  if data[0]=='AZIMUTH':
    phi_FK=float(data[1])+180.
  if data[0]=='BACK_AZIMUTH':
    phi_FK=float(data[1])
  if data[0]=='TAKE_OFF':
    theta_FK=float(data[1])
  if data[0]=='ORIGIN_WAVEFRONT':
    x0=float(data[1])
    y0=float(data[2])
    z0=float(data[3])
  if data[0]=='ORIGIN_TIME':
    t0=float(data[1])

#######################################################################################
#     Calculate the traveltime from the  starting point (O=(x0,y0,z0)) of wavefront   #
#     to the central of the study area at the surface (Q=(x0,y0,0))                   #
#######################################################################################
tt=0.
h_FK=np.zeros(nlayer)
hsum=0.
for ilayer in ilayer_fk_input:
  if ilayer < nlayer:
      h_FK[ilayer-1]=ztop_fk_input[ilayer-1]-ztop_fk_input[ilayer]
      hsum=hsum+h_FK[ilayer-1]
      if kpsv==1:
        tt=tt+h_FK[ilayer-1]/np.cos(theta_FK*np.pi/180.)/vp_fk_input[ilayer-1]
      elif kpsv==2:
        tt=tt+h_FK[ilayer-1]/np.cos(theta_FK*np.pi/180.)/vs_fk_input[ilayer-1]
  elif ilayer==nlayer:
      h_FK[ilayer-1]=ztop_fk_input[ilayer-1]-zbottom
      if kpsv==1:
        tt=tt+(np.abs(zbottom)*np.cos(theta_FK*np.pi/180.)-hsum/np.cos(theta_FK*np.pi/180.))/vp_fk_input[ilayer-1]
      elif kpsv==2:
        tt=tt+(np.abs(zbottom)*np.cos(theta_FK*np.pi/180.)-hsum/np.cos(theta_FK*np.pi/180.))/vs_fk_input[ilayer-1]

  print(h_FK[ilayer-1],vp_fk_input[ilayer-1],vs_fk_input[ilayer-1],tt)

#######################################################################################
#     Calculate the station correction time for each receiver   #
#######################################################################################
eta_p=np.cos(theta_FK/180*np.pi)/vp_fk_input[nlayer-1]
tdelay=eta_p*(0-zbottom)
print(tdelay)
xp=x0-zbottom*np.sin(theta_FK*np.pi/180.)*np.cos(theta_FK*np.pi/180.)*np.sin(phi_FK*np.pi/180.)
yp=y0-zbottom*np.sin(theta_FK*np.pi/180.)*np.cos(theta_FK*np.pi/180.)*np.cos(phi_FK*np.pi/180.)
zp=zbottom*(1-np.sin(theta_FK*np.pi/180.)**2)
A=(x0-xp)
B=(y0-yp)
C=-zp
PQ=np.sqrt((x0-xp)**2+(y0-yp)**2+(zp)**2)

print('phi=%f' %phi_FK)
print('xp,yp,zp=%f %f %f' %(xp,yp,zp) )
print('A,B,C,PQ=%f %f %f %f' %(A,B,C,PQ))
p = Proj(proj='utm', ellps='WGS84', zone='10')

fp=open(fkfile,'r')
fp1=open('correct_sta.dat','w')
for i, line in enumerate(fp):
  data=line.split()
  if i==7:
    nsta=int(data[1])
    print('nsta:%d' % nsta)
  if i>=15 and i<nsta+15:
    print(data)
    stnm=data[0]
    netwk=data[1]
    stla=data[2]  # Lat
    stlo=data[3]  # Lon
    stel=data[4]  # Elevation: m
    xi,yi=p(stlo,stla)
    zi=float(stel)
    SR=(A*(xi-xp)+B*(yi-yp)+C*(zi-zp))/PQ
    RRp=SR-PQ
    if kpsv==1:
      stc=RRp/vp_fk_input[0]
    elif kpsv==2:
      stc=RRp/vs_fk_input[0]
    print('net.stnm xi yi SR RRp tt stc %s %s %f %f %f %f %f %f\n' %(netwk,stnm,xi,yi,SR,RRp,tt,stc))
    fp1.write('%s %s %f %f %f\n' %(netwk,stnm,RRp,tt,tt+stc))
fp.close()
fp1.close()
