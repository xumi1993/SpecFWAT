# Full Waveform Adjoint Tomography (FWAT)-v1.1

**Copyright (c) 2020  Kai Wang (wangkaim8@gmail.com)**

***Please do NOT distribute the codes without permission***


The FWAT package can perform both noise FWI, teleseismic FWI and their joint inversion.
The inversion workflow is illustrated in an older version named [SPECFEM3D_ANAT](
https://github.com/wangkaim8/SPECFEM3D_ANAT).

## Acknowledgements
For the use of the FWAT package, please cite as: 

We use the FWAT package (Wang et al., 2018; Wang et al., 2021) to conduct XXX inversions in this study. The FWAT package is available from https://gitlab.com/specfem_fwat/fwat.

1. Wang, K., Yang, Y., Basini, P., Tong P., Tape, C. and Liu Q., 2018.
    Refined crustal and uppermost structure of southern California by ambient noise 
    adjoint tomography. Geophysical Journal International, 215(3), 844-863.

2. Wang, K., Yang, Y., Jiang, C., Wang, Y., Tong, P., Liu, T., & Liu, Q. (2021). 
   Adjoint tomography of ambient noise data and teleseismic P waves: Methodology and applications to central California.
   Journal of Geophysical Research: Solid Earth, 126(6), e2021JB021648.

For **ambient noise adjoint tomography**, please also consider citing:

3. Wang, K., Liu, Q. and Yang, Y. , 2019. 
    Three‐dimensional sensitivity kernels for multicomponent empirical Green's functions
    from ambient noise: Methodology and application to Adjoint tomography. 
    Journal of Geophysical Research: Solid Earth, 124(6), 5794-5810.

4. Wang, K., Jiang, C., Yang, Y., Schulte‐Pelkum, V. and Liu, Q., 2020. 
    Crustal deformation in Southern California constrained by radial anisotropy from 
    ambient noise adjoint tomography. Geophysical Research Letters, 47(12), p.e2020GL088580.

For **teleseismic full-waveform inversion**, please also consider citing:

5. Wang, K., Wang, Y., Song, X., Tong, P., Liu, Q., & Yang, Y. (2022). 
    Full‐Waveform Inversion of High‐Frequency Teleseismic Body Waves Based on 
    Multiple Plane‐Wave Incidence: Methods and Practical Applications. 
    Bulletin of the Seismological Society of America, 112(1), 118-132.

Please contact Kai Wang (wangkaim8@gmail.com) if you have any suggestions.

## Installation

### Clone source codes
```
git clone --recursive https://github.com/geodynamics/specfem3d.git
cd specfem3d
git checkout -b fwat 42abac6d18
cd src
git clone --recursive https://gitlab.com/wangkaim8/fwat.git fullwave_adjoint_tomo
cd ..
```

### Change constant parameters

- Change `setup/constants_tomography.h.in`
```fortran
logical, parameter :: USE_ALPHA_BETA_RHO = .true.
logical, parameter :: USE_ALPHA_BETA_RHO_TISO = .false.  ! Kai added
```
- Change `src/shared/shared_par.F90`
```fortran
logical :: ANISOTROPIC_KL,SAVE_TRANSVERSE_KL,APPROXIMATE_HESS_KL,SAVE_MOHO_MESH
logical :: SAVE_AZIMUTH_KL ! Kai added
```
- Change `setup/constants.h.in`
```fortran
logical, parameter :: DO_BRUTE_FORCE_POINT_SEARCH = .true.
```
### Install FLEXWIN ttimes libary

```
cd src/fullwave_adjoint_tomo/preproc_flexwin/ttimes_mod

```

- For users of gfortran

```
make -f make_gfortran 
make -f make_gfortran install
```

-  For users of Intel ifort

```
make -f make_intel 
make -f make_intel install
```

### Configure and modify Makefile

#### Modify `flag.guess`
- For Intel ifort
```
ifort|*/ifort)
...

    DEF_FFLAGS="-heap-arrays 64 -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -std03 -diag-disable 6477 -implicitnone -gen-interfaces -warn all"
```

- For gfortran

```
gfortran|*/gfortran|f95|*/f95)
...

    DEF_FFLAGS="-std=gnu -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused -O3 -finline-functions"
```

#### Run `configure`
- For Intel ifort
```
./configure FC=ifort CC=icc
```

- For gfortran
```
./configure FC=gfortran CC=gcc MPIFC=mpif90 --with-mpi
```


##### 1. Set up additional libs

- For Intel ifort
```Makefile
MPILIBS = -Lsrc/fullwave_adjoint_tomo/preproc_flexwin/ttimes_mod -ltau -lm -lmkl
```
- For gfortran
  - install BLAS, LAPACK
    - In Ubuntu, run
      ```
      sudo apt-get install libblas-dev liblapack-dev
      ```
    - In linux clusters, try to module load openblas
    - In Mac OS, 
      ```
      brew install openblas
      brew install lapack
      ```

   - Change Makefile:
        ```makefile
        MPILIBS = -Lsrc/fullwave_adjoint_tomo/preproc_flexwin/ttimes_mod -ltau -lm -llapack -lblas
        ```

##### 2. Set sub-directory of FWAT to Makefile
Replace
```
SUBDIRS = \
        tomography \
        tomography/postprocess_sensitivity_kernels \
```
with
```
SUBDIRS = \
        fullwave_adjoint_tomo/optimization \
        fullwave_adjoint_tomo \
        fullwave_adjoint_tomo/post_processing \
```
##### 3. Add targets to Makefile for compilation 
```
DEFAULT = \
        fwat_postproc \
        optimization \
        fullwave_adjoint_tomo \
```
delete
```
all: postprocess tomography
```
add
```
all: fwat_postproc optimization fullwave_adjoint_tomo
```

#### Compilation
```
make all
```

## Input/Output
### Input
- `fwat_data`: seismic data for inversion
  - sac header should include `knetwk` `kstnm` `kcmpnm` `dist` `az` `baz` `stlo` `stla`
- `fwat_params`: Parameters of FWAT
- `src_rec`: Files to store source and receiver information

```
./fwat_data/---
    evtid1---
         net.stnm1.chan.sac
         net.stnm2.chan.sac
         ...
    evtid2---
         net.stnm1.chan.sac
         ...
    ...   ---

./fwat_params/---
           FWAT.PAR
           MEASUREMENT.PAR
./src_rec/---
       sources_set1.dat
       sources_set2.dat
       sources_set3.dat
       ....

       FORCESOLUTION_evtid
       STATIONS_evtid
```
### Output
- `solver`: forward waveforms, adjoint sources and event kernels
- `misfits`: misfits
- `optimize`: summed misfit kernels and updated models
```
./solver/---
	./M00.set1---
		./GRADIENT---
		./evtid1---
			./EKERNEL
			./OUTPUT_FILES
			./SEM
		./evtid2---
			./EKERNEL
			./OUTPUT_FILES
			./SEM
		...

	./M00.set2---
		./GRADIENT---
		./evtid?---
			./EKERNEL
			./OUTPUT_FILES
			./SEM
		...
	...
		./GRADIENT---
		./evtid?---
			./EKERNEL
			./OUTPUT_FILES
			./SEM
		...
./misfits/---
	M00.set1_T006_T015_evtid1_window_chi
	M00.set1_T010_T020_evtid1_window_chi
	M00.set1_T015_T030_evtid1_window_chi
	M00.set1_T020_T040_evtid1_window_chi
	M00.set2_T006_T015_evtid1_window_chi
	...
./optimize/---
	./SUM_KERNELS_M00---
	./SUM_KERNELS_M01---
	...
	./MODEL_M00---
	./MODEL_M01---
	...
```

