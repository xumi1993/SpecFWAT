#!/usr/bin/python

import sys
import os
import numpy as np
from io import StringIO

f=open("ak135.tvel",'r')
# read lines 1-2
lay1=f.readline()
data=lay1.split()
#print('LAYER 1 '+str(float(data[3])*1000)+' '+str(float(data[1])*1000)+' '+str(float(data[2])*1000)+' 0.\n')
print("LAYER %d %.3f %.3f %.3f %.3f" % (1,float(data[3])*1000,float(data[1])*1000,float(data[2])*1000,-float(data[0])*1000+000.)) # add 5km topography
f.readline()
# read lines 3-4
lay2=f.readline()
data=lay2.split()
#print('LAYER 2 '+str(float(data[3])*1000)+' '+str(float(data[1])*1000)+' '+str(float(data[2])*1000)+' 0.\n')
print("LAYER %d %.3f %.3f %.3f %.3f" % (2,float(data[3])*1000,float(data[1])*1000,float(data[2])*1000,-float(data[0])*1000))
f.readline()
# read lines 4-8
data=np.genfromtxt(f, dtype="f8,f8,f8,f8", names=['dep', 'vp', 'vs', 'rho'],max_rows=6)
f.close()

# We use PREM as the layered model. As the velocity increases with depth in the PREM,
# we need to use the average value of PREM for each layer.
ilay=3
dep=data['dep']
vp=data['vp']
vs=data['vs']
rho=data['rho']
for i in range(len(dep)-1):
    vpi=(vp[i]+vp[i+1])/2.*1000
    vsi=(vs[i]+vs[i+1])/2.*1000
    rhoi=(rho[i]+rho[i+1])/2.*1000
    depi=-dep[i]*1000
    print("LAYER %d %.3f %.3f %.3f %.3f" %(ilay,rhoi,vpi,vsi,depi))
    ilay=ilay+1
