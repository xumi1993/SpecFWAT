# SPECFWAT

This is a package for SPECFEM3D driven Full-waveform adjoint Tomography

## Installation

To install the package, you can use the following command:

```bash
mkdir build && cd build
```

Then, you can use the following command to build the package:

```bash
cmake .. && make -j
```

Alternatively, you can specify the compiler and MPI compiler by using the following command:

```bash
CC=gcc-13 FC=gfortran-13 MPIFC=mpifort cmake .. && make -j
```

For Intel compilers, you can use the following command:

```bash
CC=icc FC=ifort MPIFC=mpiifort cmake .. && make -j
```

## Difference from the original SPECFEM3D

In this package, we have made some changes to the original SPECFEM3D code. The main changes are:

- Only the meshfem3D, generate databases, and forward/adjoint simulation part is kept.
- Use the CMake build system instead of the original Makefile system.
- Remove support for the ASDF and ADIOS format.

