# SpecFWI

Full-waveform Inversion Based on SPECFEM3D

## Installation

### Dependencies

```bash
conda install -c conda-forge fortran-compiler c-compiler cmake openmpi
```

### For HPC (BSCC-N)
```bash
module load openmpi/4.1.5_gcc11.2_ucx1.14.1_cuda11.8 cmake/3.30.0 cuda/11.8  
```

```bash
CC=gcc FC=gfortran MPIFC=mpifort cmake .. -DUSE_CUDA=ON -DFORCE_DOWNLOAD_EXTERNAL_LIBS=ON
make -j4
```