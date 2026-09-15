# SpecFWAT

 An easy, fast, and powerful full-waveform adjoint tomography (FWAT) tool for multiple seismic data.

<img width="3302" height="2968" alt="fwat" src="https://github.com/user-attachments/assets/493f3737-0d4c-4ca8-b153-8f0e96a2fcad" />


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

## Single-command GLL inversion

Build the `specfwat` target, then run from the case directory:

```bash
# In the repository, with your usual compiler/MPI/CUDA modules loaded:
cmake -S . -B build
cmake --build build --target specfwat -j 4

# In the case directory:
mpirun -np "$NPROC" /path/to/bin/xspecfwat -m M00 -s noise -n 10
```

This runs ten iterations (`M00` to `M10`) in one MPI job: mesh/database setup,
forward and adjoint simulations, kernel processing, and model updates. The mesh
is generated once; model-dependent databases are regenerated each iteration.
The original four commands remain available for regular-grid/joint inversion.

The first version supports one data type (`noise`, `tele`, or `leq`), a purely
elastic isotropic model (`MODEL_TYPE: 1`, `vp/vs/rho`), and SD or L-BFGS
(`OPT_METHOD: 1` or `2`). Use the standard internal mesher, binary SPECFEM
files, and one MPI group with the same number of ranks as `NPROC`.

### Configuration and initial model

In `DATA/fwat_params.yml`, retain the selected data section, measurement/output
settings, `POSTPROC`, and `MODEL_UPDATE`. Enable exactly the data type selected
by `-s` in `POSTPROC.INV_TYPE` (order: noise, tele, leq). For example:

```yaml
POSTPROC:
  INV_TYPE: [True, False, False]
  JOINT_WEIGHT: [1.0, 1.0, 1.0]
  IS_PRECOND: False
  # Existing taper settings may also be used.

MODEL_UPDATE:
  INIT_MODEL_PATH: initial_model.h5
  MODEL_TYPE: 1
  OPT_METHOD: 2
  ITER_START: 0
  LBFGS_M_STORE: 5
  MAX_SLEN: 0.02
  DO_LS: true
  MAX_SHRINK: 0.618
  MAX_SUB_ITER: 10
  C1: 0.01
  VPVS_RATIO_RANGE: [1.3, 2.5]
```

For a fresh run (`-m M00`), choose the initial model source with `MODEL` in
`DATA/Par_file`, following SPECFEM's `IMODEL` selection:

- `MODEL = external`: `INIT_MODEL_PATH` points to the initial HDF5 file, using
  the existing external-model format with `/x`, `/y`, `/z`, `/vp`, `/vs`, and
  `/rho` datasets. `generate_databases_fwat` interpolates this model onto the
  mesh in the solver database. At optimizer initialization, `vp/vs/rho` are
  obtained from its material arrays and saved as GLL model files. No separate
  projection or copy to `DATA/tomo_files/tomography_model.h5` is needed.
- `MODEL = gll`: provide `proc*_vp.bin`, `proc*_vs.bin`, and `proc*_rho.bin`
  directly in `LOCAL_PATH`. Database generation reads these files using
  SPECFEM's GLL loader. `INIT_MODEL_PATH` is unused and may be omitted.

At the start of optimization for `Mxx`, the current GLL model is archived in
`optimize/model_Mxx` for L-BFGS. The accepted update `Mxx+1` is written only to
`LOCAL_PATH`, ready for the next database build. Its history copy is written
when the next iteration reaches optimization. For example, after two iterations,
`model_M00` and `model_M01` exist in `optimize`, while `LOCAL_PATH` holds `M02`.
Subsequent iterations and line searches use `MODEL = gll` in memory.
`SAVE_MESH_FILES` stays false; the driver writes only the required model files.

Each model file contains one SPECFEM sequential unformatted record of
single-precision values in `(NGLLX, NGLLY, NGLLZ, NSPEC_local)` order. With
attenuation enabled, the first optimization also writes `proc*_qmu.bin`
and `proc*_qkappa.bin` using the existing external-model routine; GLL
initialization requires these files in `LOCAL_PATH`. They are archived in
`model_M00` and remain fixed during inversion.

`MODEL_GRID` is unused and may be omitted for `xspecfwat`; optimization vectors
stay on the GLL mesh. Use positive `SIGMA_H` and `SIGMA_V` and `PRECOND_TYPE: 1`,
`2`, or `3` in the selected data section.

### Files and restart

- `optimize/model_Mxx/proc*_vp.bin`, `vs.bin`, `rho.bin`: accepted model history.
- `optimize/SUM_KERNELS_Mxx/proc*_alpha_kernel_smooth.bin`, `beta_kernel_smooth.bin`,
  `rhop_kernel_smooth.bin`: processed gradient history read directly by L-BFGS.
- `optimize/SUM_KERNELS_Mxx/proc*_hess_inv.bin`: inverse diagonal preconditioner
  when `IS_PRECOND: False`; with `True`, preconditioning is applied to the kernels.
- `LOCAL_PATH/proc*_{vp,vs,rho}.bin`: latest accepted model, ready for the next
  database generation. Solver databases correspond to the last evaluated model.
- `output_optimize_gll_Mxx.log`: search direction and accepted step information.

L-BFGS uses log-model differences and MPI sums of GLL quadrature integrals
(including element Jacobians). Invalid curvature pairs are skipped and a
non-descent direction falls back to preconditioned SD. `DO_LS: true` enables
Armijo backtracking with forward evaluations; a failed search restores the
current model and database and exits with an error.

To continue after `M10` has been produced:

```bash
mpirun -np "$NPROC" /path/to/bin/xspecfwat -m M10 -s noise -n 5
```

A continuation reads the latest GLL model directly from `LOCAL_PATH`; it does
not reread the initial HDF5 file. Keep `LOCAL_PATH` and the L-BFGS history when
resuming, including the fixed quality factor files when attenuation is enabled.
To restart from an older archived model, first copy its `vp/vs/rho` files from
`optimize/model_Mxx` to `LOCAL_PATH` and pass the matching `-m Mxx`.
Keep the mesh, element ordering,
CPU/GPU mesh-coloring configuration, rank count, data selection, processing
settings, and optimization history consistent when resuming. Set `ITER_START`
to the restart index to reset L-BFGS history after changing processing settings.
Model indices are `M00`–`M99`.
