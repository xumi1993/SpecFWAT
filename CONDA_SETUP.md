# Conda Environment Setup for SpecFWAT

This directory contains configuration files and scripts for setting up a conda environment to build and run SpecFWAT.

## Files

- `environment.yml` - Conda environment specification with all required dependencies
- `setup_conda_env.sh` - Interactive script to detect conda and set up the environment
- `.github/workflows/conda-detection.yml` - GitHub Action to automatically test conda setup

## Quick Start

### Option 1: Using the setup script (recommended)
```bash
./setup_conda_env.sh
```
This script will:
- Detect if conda is installed
- Offer to install Miniconda if needed
- Create the SpecFWAT conda environment
- Test all dependencies

### Option 2: Manual setup with environment.yml
```bash
# Create environment from file
conda env create -f environment.yml

# Activate environment
conda activate specfwat-env

# Verify installation
conda list
```

### Option 3: Manual setup with individual packages
```bash
# Create new environment
conda create -n specfwat-env python=3.11

# Activate environment
conda activate specfwat-env

# Install dependencies
conda install -c conda-forge cmake make gfortran openmpi hdf5 hdf5-parallel openblas lapack yaml-cpp numpy scipy matplotlib
```

## Building SpecFWAT

After setting up the conda environment:

```bash
# Activate the environment
conda activate specfwat-env

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DUSE_HDF5=ON \
         -DUSE_EXTERNAL_LIBS=OFF

# Build
make -j$(nproc)
```

## Dependencies

The conda environment includes:

### Required for SpecFWAT
- **CMake** (≥3.20) - Build system
- **Fortran compiler** (gfortran) - For Fortran code compilation
- **MPI** (OpenMPI) - Parallel processing support
- **HDF5** (with parallel support) - Data format support
- **BLAS/LAPACK** (OpenBLAS) - Linear algebra libraries
- **yaml-cpp** - YAML configuration file parser

### Optional/Development
- **Python** (3.11) with NumPy, SciPy, matplotlib - For analysis and visualization
- **h5py** - Python HDF5 interface
- **Jupyter** - Interactive development
- **VTK** - Visualization toolkit

## Testing

The GitHub Action automatically tests:
- Conda detection and installation
- Environment creation
- Dependency verification
- Basic CMake configuration

## Troubleshooting

### Conda not found
If conda is not installed, the setup script will offer to install Miniconda automatically.

### HDF5 issues
Make sure to use `hdf5-parallel` package for MPI support:
```bash
conda install -c conda-forge hdf5-parallel
```

### MPI conflicts
If you encounter MPI-related issues, try:
```bash
conda install -c conda-forge openmpi openmpi-mpifort
```

### Building issues
If CMake cannot find dependencies, try setting explicit paths:
```bash
cmake .. -DHDF5_ROOT=$CONDA_PREFIX \
         -DMPI_ROOT=$CONDA_PREFIX \
         -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
```

## Environment Management

### List environments
```bash
conda env list
```

### Update environment
```bash
conda env update -f environment.yml
```

### Remove environment
```bash
conda env remove -n specfwat-env
```

### Export current environment
```bash
conda env export > my-environment.yml
```