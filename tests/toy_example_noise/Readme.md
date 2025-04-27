# Toy example for ambient noise adjoint tomography

This is a toy example for the ambient noise adjoint tomography.

## 1. Create initial model and target model

Run `create_ckb.ipynb` to create the initial and target models. The initial model is saved as `initial_model.h5`, and the target model is saved as `target_model.h5`.

![image](https://github.com/user-attachments/assets/865a1a2b-d6fc-47aa-a93d-3d595dab7f2f)

- Dots are stations.
- Blue dots are virtual sources.

## 2. Generate synthetic data

Run `00_forward.sh` to generate synthetic data for the target model. The synthetic data is saved in the `fwat_data` folder.

## 3. Run inversion

Run `01_run_this_test.sh` for the inversion. The inversion will run for 9 iterations. Then change `ITER_START` to 9 to restart L-BFGS. setup `for it in ``seq 9 12``; do` in the `01_run_this_test.sh` to run 9th to 12th iterations.

## 4. Check the result

Run `plot_model.ipynb` to plot the final model.

![image](https://github.com/user-attachments/assets/c1718a80-4c93-466c-85b1-dd7b040163ef)

![image](https://github.com/user-attachments/assets/5ca41223-63b9-4e46-a0c9-45b3a2052255)

![image](https://github.com/user-attachments/assets/a5c949bf-1a1a-42b3-b377-6e766ed17e5d)
