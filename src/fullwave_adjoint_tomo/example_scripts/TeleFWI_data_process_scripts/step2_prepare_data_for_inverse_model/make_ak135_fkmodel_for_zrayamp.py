#!/usr/bin/python

import sys
import os
import numpy as np
from io import StringIO

f=open("StdModels/ak135.tvel",'r')
#print("MODEL AK135")
#print(" 2 0 1")
# read lines 1-2
f.readline()
f.readline()
# read line 3-70
data=np.genfromtxt(f, dtype="f8,f8,f8,f8", names=['dep', 'vp', 'vs', 'rho'],max_rows=69)
f.close()

dep=data['dep']
vp=data['vp']
vs=data['vs']
rho=data['rho']
dep_new=np.zeros(len(dep))
vp_min=np.zeros(len(dep))
vp_max=np.zeros(len(dep))
rho_min=np.zeros(len(dep))
rho_max=np.zeros(len(dep))
vs2vp_min=np.zeros(len(dep))
vs2vp_max=np.zeros(len(dep))
j=0
for i in range(len(dep)-1):
    if i == 0:
       dep_new[j]=dep[i]
       vp_min[j]=0.
       vp_max[j]=vp[i]
       rho_min[j]=0.
       rho_max[j]=rho[i]
       vs2vp_min[j]=0.
       vs2vp_max[j]=vs[i]/vp[i]
       j=j+1
    if i > 0:
       if dep[i] == dep[i-1]:
          vp_min[j-1]=vp[i-1] 
          vp_max[j-1]=vp[i] 
          rho_min[j-1]=rho[i-1]
          rho_max[j-1]=rho[i]
          vs2vp_min[j-1]=vs[i-1]/vp[i-1]
          vs2vp_max[j-1]=vs[i]/vp[i]
          #print("Interface at %f") %(dep[i],)
       else:
          dep_new[j]=dep[i]
          vp_min[j]=0.
          vp_max[j]=vp[i]
          rho_min[j]=0.
          rho_max[j]=0.
          vs2vp_min[j]=0.
          vs2vp_max[j]=vs[i]/vp[i]
          j=j+1
         
for i in range(j):
    print(" %6.1f %.4f %.4f %.4f %.4f %.4f %.4f %.1f") %(dep_new[i],vp_min[i],vp_max[i],rho_min[i],rho_max[i],vs2vp_min[i],vs2vp_max[i],0.)
    #ilay=ilay+1
print(" %6.1f %.4f %.4f %.4f %.4f %.4f %.4f %.1f") %(-1,0.,0.,0.,0.,0.,0.,0.)
