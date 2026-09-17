# Standalone FK wavefield generator

`xfwat_fk` computes incident boundary wavefields from the existing mesh using the
normal FWAT parameters, event list and FKmodels. It reuses the project's parameter
readers, mesh readers, logger and HDF5 cache writer. No SEM time stepping is run.
The shared `fk_coupling::compute_fk_wavefield(evtid)` entry point in
`src/fk/fk_coupling.f90` owns preparation, CPU/GPU dispatch and optional saving.
Both this driver and preproc call it. GPU code and the generated CPU reference are
compiled once into the shared FWAT library, outside the SPECFEM submodule.

`prepare_fk.f90` defines the `PrepareFK` type. Its `init()` reads FWAT parameters,
validates the selected event range and loads the mesh once. The small driver in
`fwat_fk_compute.f90` loops over `first_event:last_event`, calling
`prepare_for_event(ievt)` for each event, then `finalize()` closes logs and releases
acquisition resources. With `-e`, the range contains only that event. MPI startup
and shutdown remain in the driver; the executable is still named `xfwat_fk`.

## Build and run

Build with the project's usual MPI/Fortran/HDF5 environment:

```sh
cmake -S . -B build -DNO_UPDATE_SUBMODULE=ON -DUSE_CUDA=ON
cmake --build build --target fwat_fk -j 8
```

CUDA builds link cuFFT. `-DUSE_CUDA=OFF` builds the CPU backend alone and requires
`GPU_MODE = .false.` in `DATA/Par_file`. The Fortran compiler must match MPI and
HDF5; mesh databases and the consuming solver must use the same precision.

Run from the normal FWAT case directory:

```sh
# Compute all events in the configured teleseismic source list.
mpirun -np 8 /path/to/SpecFWI/bin/xfwat_fk

# Compute only the first event (the same 1-based index convention as FWAT).
mpirun -np 8 /path/to/SpecFWI/bin/xfwat_fk -e 1
```

`-e`/`--event` is the only computation-selection option. All other settings are
inherited through `fpar%read()`, `read_parameter_file()` and
`fpar%select_simu_type()`:

| Setting | Existing source |
| --- | --- |
| Time step and number of steps | `TELE.DT` and `TELE.NSTEP` in `DATA/fwat_params.yml`, overriding `Par_file` |
| CPU or GPU computation | `GPU_MODE` in `DATA/Par_file` |
| Mesh directory | FWAT's `local_path_fwat`, derived from `LOCAL_PATH` and joint-inversion settings |
| Event list | `TELE.TELE_TYPE`: `src_rec/sources_tele.dat`, `sources_rf.dat` or `sources_telecc.dat` |
| Incident-wave model | `src_rec/FKmodel_<event ID>`, read by the existing `read_fk_model()` |
| Save results | `TELE.SAVE_FK` |
| HDF5 compression | `TELE.COMPRESS_LEVEL` (0 through 9) |
| Cache root | `local_path_backup`, the original `LOCAL_PATH` from `Par_file` |
| GPU device | FWAT's `noderank` modulo the number of devices visible through `CUDA_VISIBLE_DEVICES` |

Set `TELE.INJECTION_TYPE` to FK (1, the default). No new YAML fields are required.
GPU batching is handled internally according to available device memory.

Mesh initialization uses `initialize_simulation_fwat()` and
`read_mesh_databases_fwat()`. Keep the normal inputs, including
`OUTPUT_FILES/values_from_mesher.h` and `DATA/CMTSOLUTION`. The MPI rank count must
match `NPROC` and the existing mesh partition. The standard reader allocates host
arrays; the standalone program avoids SEM GPU initialization. Station files are
not read by this program.

Initialization writes `OUTPUT_FILES/output_solver.txt`; the project logger records
timestamped progress in `OUTPUT_FILES/output_fk.log`. The CPU reference also emits
its small `OUTPUT_FILES/plot_FK_*` diagnostics.

When `TELE.SAVE_FK` is true, output is
`LOCAL_PATH/FK_wavefield_<event ID>/procXXXXXX_fk_wavefield.h5`. In joint inversion,
the mesh can be under `LOCAL_PATH/tele` while caches remain under `LOCAL_PATH`.
The shared writer preserves boundary-face/GLL order, padding, spline coefficients
and precision. Partitions too small for its fixed compression chunks use its
uncompressed branch. Each event's arrays are released before processing the next.
When `SAVE_FK` is false, the program computes the wavefields without writing caches.

Re-running an event with saving enabled replaces its cache. Regenerate after
changing the FK model, mesh/materials/partition, `TELE.DT` or `TELE.NSTEP`.
For CPU/GPU comparisons, run once with each `GPU_MODE` setting and copy the first
cache elsewhere before the second run.

## Preproc integration

`prepare_timerun_fwat()` uses the same shared interface when an event has no FK
cache. With `GPU_MODE = .true.`, both FK and the subsequent SEM simulation use
GPU execution. `TELE.SAVE_FK = .false.` still allows GPU FK computation; it only
suppresses writing a new cache. Existing caches are read through the original
`read_fk_coupling_file()` path regardless of this setting.

The shared computation leaves `Veloc_FK`, `Tract_FK` and `ipt_table` available for
SEM injection, replacing them for each new event. CUDA FK frees its temporary
workspace and restores the solver's active device before SEM preparation resumes.
The original FKmodel reader and HDF5 read/write routines are unchanged.

## Implementation and numerical compatibility

- CPU: CMake extracts the existing `FK3D`/`FK` routines into a generated build
  file and renames their entry points. It makes one numerical correction:
  initialize `eta_p = eta_s` for SV, since the original code uses uninitialized
  `eta_p` for its phase delay. P-wave calculations are unchanged. No submodule
  source is edited. The SV fix now applies to both standalone and preproc FK.
- GPU: the small, frequency-dependent reflection solve and layer-interface
  states are prepared on CPU. CUDA evaluates each point/frequency's remaining
  propagation, followed by batched double-complex cuFFT and per-trace spline
  filtering. Full-layer propagation is reused across points. Results are copied
  back into the existing Fortran arrays and passed to the shared HDF5 writer.
- GPU batches are internal workspace management: at most 1024 points, reduced
  automatically if the reported free device memory is insufficient. The full output remains in host
  memory; cuFFT additionally needs plan workspace.
- FFT direction, `1/(N*DT)` normalization, zero Nyquist, time shift and taper
  follow the reference. Elastic traction spline filtering intentionally retains
  the reference's vertical-velocity scratch tail; changing that numerical
  convention is outside this implementation.
- Tiny cases can be slower on GPU because context/FFT-plan setup dominates.
  Speedups must be measured at your real boundary-point and frequency counts.

Mesh format selection follows the existing reader and `Par_file` flags. The regression
suite covers binary `procXXXXXX_external_mesh.bin`. HDF5/ADIOS mesh input depends
on the existing SPECFEM build support; it is not verified here. Dedicated HDF5 I/O ranks are not supported (`HDF5_IO_NODES` must be zero).
FK supports elastic or acoustic/elastic domains, P or SV incidence and propagating
slowness. Positive `ORIGIN_TIME`, critical/evanescent slowness and
poroelastic meshes are rejected, because the reference does not handle these
reliably. Negative time shifts must fit within both `NSTEP` and the FFT length.

## Validation

Compare two saved runs with NumPy and h5py:

```sh
python tests/fk/compare.py fk_cpu/FK_wavefield_0 fk_gpu/FK_wavefield_0
```

The comparison checks partitions, shape, dtype, finite values, relative L2 and
peak errors. Its default tolerance is `1e-4`, not a claim of bitwise identity.
Even double-storage reference FFTs include default-kind complex rounding.

Run the regression suite on complete meshes generated by the standard tools:

```sh
cmake --build build --target fwat_fk fwat_fwd_measure_adj meshfem3D generate_databases -j 8
python tests/fk/run.py --exe bin/xfwat_fk --gpu --preproc-exe bin/xfwat_fwd_measure_adj
# Without CUDA:
python tests/fk/run.py --exe bin/xfwat_fk
```

Fixtures deliberately give `Par_file` different time settings from the TELE
section to check inheritance. Tests cover event selection and consecutive events,
joint-inversion paths, all three teleseismic source lists, saving and compression
settings, CPU/GPU selection, P/SV waves, fluid/solid interfaces, negative origin
time, partial GPU batches, multiple MPI ranks, zero padding and invalid options.
Use `--mesh-bin-dir` if the mesh executables are not beside `xfwat_fk`.

With `--preproc-exe`, tests also run complete CPU/GPU forward simulations over two
events, compare standalone and preproc caches, verify cache reuse, and verify
matching receiver traces with saving disabled. FK caches are compared exactly;
GPU SEM traces allow `1e-4` relative L2/peak error because atomic accumulation can
change floating-point rounding between runs.
