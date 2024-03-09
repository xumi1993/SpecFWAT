!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! add_model_iso
!
! this program can be used to update ISOTROPIC model files with
! (smoothed & summed) event kernels.
! the kernels are given for isotropic parameters (alpha,beta,rho) or ( bulk_c,beta,rho).
!
! the algorithm uses a steepest descent method with a step length
! determined by the given maximum update percentage.
!
! input:
!    - step_fac : step length to update the models, f.e. 0.03 for plusminus 3%
!
! setup:
!
!- INPUT_MODEL/  contains:
!       proc000***_vs.bin &
!       proc000***_vp.bin &
!       proc000***_rho.bin
!
!- INPUT_GRADIENT/ contains:
!       proc000***_bulk_c_kernel_smooth.bin &
!       proc000***_bulk_beta_kernel_smooth.bin &
!       proc000***_rho_kernel_smooth.bin
!     or
!       proc000***_alpha_kernel_smooth.bin &
!       proc000***_beta_kernel_smooth.bin &
!       proc000***_rho_kernel_smooth.bin
!
!- topo/ contains:
!       proc000***_solver_data.bin
!
! new models are stored in
!- OUTPUT_MODEL/ as
!   proc000***_vp_new.bin and
!   proc000***_vs_new.bin and
!   proc000***_rho_new.bin and
!
! USAGE: e.g. ./add_model_iso 0.3


program add_model

  use tomography_model_tiso
  use tomography_kernels_tiso
  use tomography_kernels_tiso_cg

  use specfem_par, only: NPROC
  ! use fullwave_adjoint_tomo_par, only: USE_RHO_SCALING_FWAT
  use fwat_input, only: tomo_par

  implicit none

  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: betav1,betav0,betah1,betah0,rho1,rho0, &
      betaiso1,betaiso0,alphav1,alphav0,alphah1,alphah0,eta1,eta0
  real(kind=CUSTOM_REAL) :: dbetaiso,dbulk


  ! ============ program starts here =====================

  ! initializes arrays
  call initialize()

  ! reads in parameters needed
  call read_parameters_tomo()

  ! user output
  if (myrank == 0) then
    print *
    print *,'***********'
    print *,'program add_model_tiso_cg: '
    print *,'  NPROC: ',NPROC
    print *
    print *,'model update for vsv, vsh & vpv, vph & rho & eta'
    print *,'  step_fac = ',step_fac
    print *
    if (USE_ALPHA_BETA_RHO_TISO) then
      print *,'kernel parameterization: (alphav,alphah,betav,betah,rho,eta)'
    else
      print *,'kernel parameterization: (bulk,betav,betah,rho,eta)'
    endif
    print *
    if (tomo_par%USE_RHO_SCALING_FWAT) then
      print *,'scaling rho perturbations'
      print *
    endif
    print *,'***********'
    print *
  endif

  ! reads in current isotropic model files: vp & vs & rho
  call read_model_tiso()

  ! reads in smoothed kernels: bulk, beta, rho
  call read_kernels_tiso()

  ! reads in old (former inversion) smoothed kernels: alphav, alphah, betav, betah, eta
  call read_kernels_cg_tiso_old()

  ! calculates gradient
  ! conjungate gradient method
  call get_cg_direction_tiso()

  ! computes new model values for alpha, beta and rho
  ! and stores new model files
  ! allocate new model arrays
  allocate(model_vpv_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_vph_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_vsv_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_vsh_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_eta_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_rho_new(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating model arrays'

  ! initializes arrays
  model_vpv_new = 0.0_CUSTOM_REAL
  model_vph_new = 0.0_CUSTOM_REAL
  model_vsv_new = 0.0_CUSTOM_REAL
  model_vsh_new = 0.0_CUSTOM_REAL
  model_eta_new = 0.0_CUSTOM_REAL
  model_rho_new = 0.0_CUSTOM_REAL

  ! model update:
  !   isotropic update everywhere
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! initial model values
          betav0 = model_vsv(i,j,k,ispec)
          betah0 = model_vsh(i,j,k,ispec)
          rho0 = model_rho(i,j,k,ispec)
          eta0 = model_eta(i,j,k,ispec)
          alphav0 = model_vpv(i,j,k,ispec)
          alphah0 = model_vph(i,j,k,ispec)

          betav1 = 0._CUSTOM_REAL
          betah1 = 0._CUSTOM_REAL
          eta1 = 0._CUSTOM_REAL
          rho1 = 0._CUSTOM_REAL
          alphav1 = 0._CUSTOM_REAL
          alphah1 = 0._CUSTOM_REAL

          ! transverse isotropic model update

          ! eta value : limits updated values for eta range constraint
          eta1 = eta0 * exp( model_deta(i,j,k,ispec) )
          if (eta1 < LIMIT_ETA_MIN ) eta1 = LIMIT_ETA_MIN
          if (eta1 > LIMIT_ETA_MAX ) eta1 = LIMIT_ETA_MAX

 
          ! shear values
          betav1 = betav0 * exp( model_dbetav(i,j,k,ispec) )
          betah1 = betah0 * exp( model_dbetah(i,j,k,ispec) )

          ! density: uses scaling relation with Voigt average of shear perturbations
          betaiso0 = sqrt(  ( 2.0 * betav0**2 + betah0**2 ) / 3.0 )
          betaiso1 = sqrt(  ( 2.0 * betav1**2 + betah1**2 ) / 3.0 )
          dbetaiso = log( betaiso1 / betaiso0 )
          rho1 = rho0 * exp( RHO_SCALING * dbetaiso )

         ! alpha values
          if (USE_ALPHA_BETA_RHO_TISO) then
             alphav1 = alphav0 * exp( model_dalphav(i,j,k,ispec) )
             alphah1 = alphah0 * exp( model_dalphah(i,j,k,ispec) )
          else
             dbulk = model_dbulk(i,j,k,ispec)
             alphav1 = sqrt( alphav0**2 * exp(2.0*dbulk) &
                            + FOUR_THIRDS * betav0**2 * ( &
                                exp(2.0*model_dbetav(i,j,k,ispec)) -exp(2.0*dbulk) ) )
             alphah1 = sqrt( alphah0**2 * exp(2.0*dbulk) &
                            + FOUR_THIRDS * betah0**2 * ( &
                                exp(2.0*model_dbetah(i,j,k,ispec)) - exp(2.0*dbulk) ) )
          endif
             

          ! stores new model values
          model_vpv_new(i,j,k,ispec) = alphav1
          model_vph_new(i,j,k,ispec) = alphah1
          model_vsv_new(i,j,k,ispec) = betav1
          model_vsh_new(i,j,k,ispec) = betah1
          model_eta_new(i,j,k,ispec) = eta1
          model_rho_new(i,j,k,ispec) = rho1

        enddo
      enddo
    enddo
  enddo
  call synchronize_all()

  ! stores new model in files
  call write_new_model_tiso()

  ! stores relative model perturbations
  call write_new_model_perturbations_tiso()

  ! computes volume element associated with points, calculates kernel integral for statistics
  call compute_kernel_integral_tiso()

  ! stop all the MPI processes, and exit
  call finalize_mpi()

end program add_model

!
!-------------------------------------------------------------------------------------------------
!

subroutine initialize()

! initializes arrays

  use tomography_par

  use specfem_par, only: NSPEC_AB,NGLOB_AB,NPROC,ADIOS_ENABLED

  implicit none

  logical :: BROADCAST_AFTER_READ

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! reads the parameter file
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  if (ADIOS_ENABLED) stop 'Flag ADIOS_ENABLED not supported yet for xadd_model, please rerun program...'

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *, 'Error number of processors supposed to run on: ',NPROC
      print *, 'Error number of MPI processors actually run on: ',sizeprocs
      print *
      print *, 'please rerun with: mpirun -np ',NPROC,' bin/xadd_model .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif

  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some array sizes below
  call read_mesh_for_init()

  ! sets tomography array dimensions
  NSPEC = NSPEC_AB
  NGLOB = NGLOB_AB

end subroutine initialize

