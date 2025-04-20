# Toy example for teleseismic FWI

## Background (starting) model

The background model is a two-layer 1D model with a velocity as follows.

| Depth (m) | Vp (km/s) | Vs (km/s) | Density (g/cm³) |
|-----------|-----------|-----------|-----------------|
| 0         | 6375.     | 3750.     | 2625.           |
| -40000    | 7616.     | 4480.     | 3136            |


## Target model

We add trigonometric perturbations to the background model to create a target model with a maximum perturbation of 10% in Vp, Vs and Rho.
![image](https://github.com/user-attachments/assets/809a5f97-0a74-42df-a809-d9bc2dd03fcd)


## Forward modeling

Run `00_forward.sh` to generate the synthetic data. The script will run the following steps:

1. Copy the target model `target_model.h5` to the `DATA/tomo_files/tomography_model.h5`.
2. Create mesh and databases for SPECFEM3D.
3. Do the forward modeling for events in the `src_rec/sources_telecc.dat`

## Inversion

Run `01_run_this_test.sh` to run the inversion for 9 iterations.
![image](https://github.com/user-attachments/assets/26949152-86f4-4040-b7a5-642cc493adfc)
