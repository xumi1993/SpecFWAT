#!/bin/bash -l

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

cat > grid.txt<<eof
0.2 -0.2 -80000
0.2 0.4 20000
7  1   4
eof

pert=0.1
mkdir -p model_pert
for name in vp vs rho;do
    mpirun -np $NPROC ../../bin/xadd_pert_trigo grid.txt $name ./DATABASES_MPI model_pert true $pert
    mpirun -np $NPROC ../../bin/xmodel_grid_cart 1000 1000 100 $name model_pert/ model_pert/ pert true
done


python << eof

import h5py

with h5py.File('model_pert/vp_pert.h5', 'r') as f:
    vp = f['vp'][:]
    x = f['x'][:]
    y = f['y'][:]
    z = f['z'][:]

with h5py.File('model_pert/vs_pert.h5', 'r') as f:
    vs = f['vs'][:]

with h5py.File('model_pert/rho_pert.h5', 'r') as f:
    rho = f['rho'][:]

with h5py.File('DATA/target_model.h5', 'w') as f:
    f.create_dataset('vp', data=vp)
    f.create_dataset('vs', data=vs)
    f.create_dataset('rho', data=rho)
    f.create_dataset('x', data=x)
    f.create_dataset('y', data=y)
    f.create_dataset('z', data=z)
eof
rm grid.txt