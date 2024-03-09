#!/usr/bin/env python

import sys
import os
import numpy as np
from io import StringIO
from pyproj import Proj

###
#compute travetime of code phases for FK method in layered halfspace model
def traveltime(H,p,al,be,nlayer,x0,y0,z0,x,y,z,phi_FK):
  t=np.zeros(7)
  t[:]=p*(xi-x0)*np.cos(phi_FK/180.*np.pi)+p*(yi-y0)*np.sin(phi_FK/180*np.pi)
  # figure out the lcoation of z with respect to layer stack
  if z<z0:
    print('z should be greater that z0')
    exit()
  if z<=0: # in lower halfspace
    print('z in the lower half space')
    exit()
  else:
    h=np.zeros(nlayer)
    ilayer=nlayer-1
    for j in range(nlayer-2,-1,-1):
      if z <= np.sum(H[j:nlayer-1]):
        ilayer = j
        break
    h[ilayer+1:nlayer-1]=H[ilayer+1:nlayer-1]
    h[ilayer] = z - np.sum(H[ilayer+1:nlayer-1])
    h[nlayer-1] = 0 - z0
    if h[ilayer] < 0:
      print('Error setting layer thickness')
      exit()
    print('z, ilayer = %s %s' %(str(z),str(ilayer)))
  ## Calculate the predict times of coda phases
  # (1) P
  for j in range(nlayer-1,ilayer-1,-1):
    eta_al = np.sqrt(1 / al[j]**2 - p**2)
    t[0] = t[0] + eta_al * h[j]
  # (2) Ps
  for j in range(nlayer-1,ilayer-1,-1):
    eta_al = np.sqrt(1 / al[j]**2 - p**2)
    eta_be = np.sqrt(1 / be[j]**2 - p**2)
    if j > ilayer:
      t[1] = t[1] + eta_al * h[j]
    else:
      t[1] = t[1] + eta_be * h[j]
  eta_al = np.sqrt(1 / al[ilayer]**2 - p**2)
  eta_be = np.sqrt(1 / be[ilayer]**2 - p**2)
  # (3) PpPmp
  t[2] = t[0] + 2 * eta_al * h[ilayer]
  # (4) PpPms, PsPmP, PpSmp
  t[3] = t[0] + eta_al * h[ilayer] + eta_be * h[ilayer]
  # (5) PsSmp, PsPms, PpSms
  t[4] = t[1] + eta_al * h[ilayer] + eta_be * h[ilayer]
  # (6) PpPmpPmp
  t[5] = t[2] + 2 * eta_al * h[ilayer]
  # (7) PsSms
  t[6] = t[1] + 2 * eta_be * h[ilayer]
  print(t)
  return t;
########################################################

if len(sys.argv) != 3:
  print("Usage: correct_wavefront_time.py FKmodel FKfile")
  exit()
else:
  print(str(sys.argv))

fkmodel=sys.argv[1]
fkfile=sys.argv[2]

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
    phi_FK=90.0-float(data[1])
  if data[0]=='BACK_AZIMUTH':
    phi_FK=-float(data[1])-90.
  if data[0]=='TAKE_OFF':
    theta_FK=float(data[1])
  if data[0]=='ORIGIN_WAVEFRONT':
    x0=float(data[1])
    y0=float(data[2])
    z0=float(data[3])
  if data[0]=='ORIGIN_TIME':
    t0=float(data[1])
  zbot_fk_input=np.zeros(nlayer)
  for i in range(nlayer-1):
    zbot_fk_input[i]=ztop_fk_input[i+1]
  zbot_fk_input[nlayer-1]=ztop_fk_input[nlayer-1]
  H=ztop_fk_input-zbot_fk_input
  z_ref_fk=zbot_fk_input[nlayer-1]
#######################################################################################
#     Calculate the traveltime from the  starting point (O=(x0,y0,z0)) of wavefront   #
#     to the central of the study area at the surface (Q=(x0,y0,0))                   #
#######################################################################################
tdelay=0.
is_utm=1
#     Calculate the station correction time for each receiver   #
if kpsv==1:
  p=np.sin(theta_FK/180*np.pi)/vp_fk_input[nlayer-1]
  eta_p=np.cos(theta_FK/180*np.pi)/vp_fk_input[nlayer-1]
elif kpsv==2:
  p=np.sin(theta_FK/180*np.pi)/vs_fk_input[nlayer-1]
  eta_p=np.cos(theta_FK/180*np.pi)/vs_fk_input[nlayer-1]

utm2latlon = Proj(proj='utm', ellps='WGS84', zone='10')
fp=open(fkfile,'r')
fp1=open('FKtimes.dat','w')
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
    if is_utm==1:
      xi,yi=utm2latlon(stlo,stla)
      zi=float(stel)
    else:
      xi=float(stla)
      yi=float(stlo)
    zi=0 #float(stel)
    dist=np.sqrt((xi-x0)**2+(yi-y0)**2)
    if kpsv==1:
      #tdelay=p*(xi-x0)*np.cos(phi_FK/180.*np.pi)+p*(yi-y0)*np.sin(phi_FK/180*np.pi)+eta_p*(0-z0)
      tdelay=traveltime(H, p, vp_fk_input, vs_fk_input, nlayer, x0, y0, z0-z_ref_fk, xi, yi, zi-z_ref_fk, phi_FK)
    elif kpsv==2:
      tdelay=p*(xi-x0)*np.cos(phi_FK/180.*np.pi)+p*(yi-y0)*np.sin(phi_FK/180*np.pi)+eta_s*(0-z0)
    print('net.stnm xi yi tdelay %s %s %f %f %f\n' %(netwk,stnm,xi,yi,tdelay[0]))
    fp1.write('%s %s %f %f %f %f %f %f %f %f\n' %(netwk,stnm,dist,tdelay[0],tdelay[1],tdelay[2],tdelay[3],tdelay[4],tdelay[5],tdelay[6]))
fp.close()
fp1.close()
