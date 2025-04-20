# SpecFWI

Full-waveform Inversion Based on SPECFEM3D

## Installation

### For PC
```bash
conda install -c conda-forge fortran-compiler c-compiler cmake openmpi hdf5
```

```bash
cmake .. && make -j
```

### For HPC 

#### N40@BSCC (NVIDIA RTX 4090)
```bash
module load openmpi/4.1.5_gcc11.2_ucx1.14.1_cuda11.8 cmake/3.30.0 cuda/11.8  
```

```bash
CC=gcc FC=gfortran MPIFC=mpifort cmake .. -DUSE_CUDA=ON -DCUDA_ARCH=10
make -j4
```

#### Beluga@CCDB (NVIDIA V100)
```bash
module load cuda/12.2 hdf5/1.14.2
```

```bash
cmake .. -DUSE_CUDA=ON -DCUDA_ARCH=9
make -j4
```