#!/usr/bin/env phthon3
import sys
import numpy as np
from scipy.io import FortranFile

#=======================================
def read_binary_data( filename, nrec, nt, data ):
    print('Read binary data: ')
    print('filename, nrec, nt = ',filename,nrec,nt)
    f=open(filename,'rb')
    for i in range(nt):
        record=np.fromfile(f,dtype='float32',count=nrec)
        data[:,i]=record
    f.close()
    print(data.size)
    return
#=======================================

def main():

    if len(sys.argv) != 6:
        print("Usage: convert_ather_binary_to_sac.py filename comp nrec dt nt ")
        exit()
    else:
        print(str(sys.argv))

    fname=sys.argv[1]
    comp=sys.argv[2]
    nrec=int(sys.argv[3])
    dt=float(sys.argv[4])
    nt=int(sys.argv[5])
    data=np.zeros((nrec,nt),dtype="f4")
    read_binary_data(fname,nrec,nt,data)
    for i in range(nrec):
        filename="dat%04d%s.semd" % (i,comp)
        f=open(filename,'w')
        for j in range(nt):
            t=-10.+dt*(j-1)
            f.write("%f %f\n" %(t,data[i,j]))
        f.close()

if __name__ == '__main__':
    main()

