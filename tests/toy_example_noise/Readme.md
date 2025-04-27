# Toy example for ambient noise adjoint tomography

This is a toy example for the ambient noise adjoint tomography.

## 1. Create initial model and target model

Run `create_ckb.ipynb` to create the initial model and target model. The initial model is saved as `initial_model.h5`, and the target model is saved as `target_model.h5`.

## 2. Generate synthetic data

Run `00_forward.sh` to generate synthetic data for target model. The synthetic data is saved in the `fwat_data` folder.

## 3. Run inversion

Run `01_run_this_test.sh` for the inversion. The inversion will run for 9 iterations. Then change `ITER_START` to 9 to restart L-BFGS. setup `for it in ``seq 9 12``; do`

## 4. Check the result

Run `plot_model.ipynb` to plot the final model.

