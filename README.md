# SpecFWAT

 An easy, fast, and powerful full-waveform adjoint tomography (FWAT) tool for multiple seismic data.

## Installation

See [Installation Guide](https://specfwat.xumijian.me/docs/installation/download) to build SpecFWAT on local machine and HPC system

## Quick Example

```bash
for it in `seq 0 9`; do
    model=`printf "M%02d" $it`
    if [ $it -eq 0 ]; then
        cp initial_model.h5 DATA/tomo_files/tomography_model.h5
    fi
    mpirun -np $NPROC ../../bin/xfwat_mesh_databases -s tele
    mpirun -np $NPROC ../../bin/xfwat_fwd_measure_adj -m $model -s tele -r 3
    mpirun -np $NPROC ../../bin/xfwat_post_proc -m $model
    mpirun -np $NPROC ../../bin/xfwat_optimize -m $model
done
```
