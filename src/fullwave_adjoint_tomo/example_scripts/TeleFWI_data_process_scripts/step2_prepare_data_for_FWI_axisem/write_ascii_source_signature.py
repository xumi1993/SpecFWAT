#!/usr/bin/env python
import sys
import os
import numpy as np
from io import StringIO
from obspy import read

#===================================================
# Write to ascii data
def write_ascii_source_signature(filename,nt,stf):
  print('write ascii source signature for source time function:')
  print('filename, nt = ',filename,nt)

  f=open(filename,'w')
  for i in range(3):
      np.savetxt(f,stf[i,:])
  f.close()
  return
#===================================================

if len(sys.argv) != 3:
  print("Usage: write_ascii_source_signature sacfile outfile ")
  exit()
else:
  print(str(sys.argv))

sacf=sys.argv[1]
outf=sys.argv[2]
st=read(sacf)
tr=st[0]
nt=tr.stats.npts
stf=np.zeros((3,nt),dtype='float32')
stf[0,:]=tr.data
stf[1,:]=tr.data
stf[2,:]=tr.data
write_ascii_source_signature(outf,nt,stf)
