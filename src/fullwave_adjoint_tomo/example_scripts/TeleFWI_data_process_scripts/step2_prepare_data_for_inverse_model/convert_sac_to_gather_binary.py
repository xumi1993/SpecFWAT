#!/usr/bin/env python

import sys
import os
import numpy as np
from io import StringIO
from obspy import read

#===================================================
# Write to binary data
def write_binary_data(filename,nrec,nt,data):
  print('Write binary data:')
  print('filename, nrec, nt = ',filename,nrec,nt)

  f=open(filename,'wb')
  for i in range(nt):
    data[:,i].tofile(f)
  f.close()
  return
#===================================================

if len(sys.argv) != 7:
  print("Usage: convert_sac_to_gather_binary.py station.lst nrec dt nt component outfile ")
  exit()
else:
  print(str(sys.argv))

stlst=sys.argv[1]
nrec=int(sys.argv[2])
dt=float(sys.argv[3])
nt=int(sys.argv[4])
comp=sys.argv[5]
outfile=sys.argv[6]
#print('nrec,dt,nt,comp=',nrec,dt,nt,comp)
trace=np.zeros((nrec,nt),dtype='float32')
irec=0
f=open(stlst,'r')
for line in f:
  data=line.split()
  sacf=data[0]
  st=read(sacf)
  tr=st[0]
  nt0=tr.stats.npts
  dt0=tr.stats.delta
  if nt0!=nt:
    print("Error! NT not match !!!")
  if np.abs(dt0-dt)>1.e10-10:
    print("Error! DT not match !!!")
  trace[irec,:]=tr.data
  irec=irec+1

write_binary_data(outfile,nrec,nt,trace)

