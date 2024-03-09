!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
subroutine model_update_opt()


  ! use fullwave_adjoint_tomo_par, only: OPT_METHOD, USE_EMPIRICAL_VP, USE_RHO_SCALING_FWAT
  use fwat_input, only: tomo_par
  use tomography_model_iso 
  use tomography_kernels_iso
  use tomography_kernels_iso_cg
 ! use tomography_model_tiso
  !use tomography_kernels_tiso
  !use tomography_kernels_tiso_cg

  use specfem_par, only: NPROC,MAX_STRING_LEN,NSPEC_AB,NGLOB_AB

  implicit none

  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: beta1,beta0,rho1,rho0,alpha1,alpha0,beta_km
  real(kind=CUSTOM_REAL) :: dbetaiso,dbulk
  !logical :: BROADCAST_AFTER_READ

  ! ============ program starts here =====================

  ! initializes arrays
  !call initialize()
  ! reads the parameter file
  !BROADCAST_AFTER_READ = .true.
  !call read_parameter_file(myrank,BROADCAST_AFTER_READ)
  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some array sizes below
  !call read_mesh_for_init()

  ! sets tomography array dimensions
  NSPEC = NSPEC_AB
  NGLOB = NGLOB_AB


  ! reads in parameters needed
  !call read_parameters_tomo()

  ! user output
  if (myrank == 0) then
    print *
    print *,'***********'
    print *,'program add_mode_iso_'//trim(tomo_par%OPT_METHOD)//': '
    print *,'  NPROC: ',NPROC
    print *
    print *,'model update for vs & vp & rho'
    print *,'  step_fac = ',step_fac
    print *,'  iter_start = ',iter_start
    print *,'INPUT_MODEL_DIR=',trim(INPUT_MODEL_DIR)
    print *,'OUTPUT_MODEL_DIR=',trim(OUTPUT_MODEL_DIR)
    print *,'KERNEL_OLD_DIR=',trim(KERNEL_OLD_DIR)
    print *,'INPUT_KERNELS_DIR=',trim(INPUT_KERNELS_DIR)
    print *
    if (USE_ALPHA_BETA_RHO) then
      print *,'kernel parameterization: (alpha,beta,rho)'
    else
      print *,'kernel parameterization: (bulk,beta,rho)'
    endif
    print *
    if (tomo_par%USE_RHO_SCALING_FWAT) then
      print *,'scaling rho perturbations'
      print *
    endif
    ! if (USE_EMPIRICAL_VP) then
    !   print *, 'use empirical vp'
    !   print *
    ! endif
    print *,'***********'
    print *
  endif

  ! reads in current isotropic model files: vp & vs & rho
  call read_model_iso()

  ! reads in smoothed kernels: bulk, beta, rho
  call read_kernels_iso()
  if (trim(tomo_par%OPT_METHOD)=='CG'.and. iter_current-iter_start>0) then
     call read_kernels_cg_iso_old()
  endif

  ! calculates gradient
  ! steepest descent method
  if (trim(tomo_par%OPT_METHOD)=='SD' .or. iter_current-iter_start==0) then
     call get_sd_direction_iso()
  elseif (trim(tomo_par%OPT_METHOD)=='CG') then 
     call get_cg_direction_iso()
  else
     call get_lbfgs_direction_iso()
  endif
  ! computes new model values for alpha, beta and rho
  ! and stores new model files
  ! allocate new model arrays
  allocate(model_vp_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_vs_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_rho_new(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating model arrays'

  ! initializes arrays
  model_vp_new = 0.0_CUSTOM_REAL
  model_vs_new = 0.0_CUSTOM_REAL
  model_rho_new = 0.0_CUSTOM_REAL

  ! model update:
  !   isotropic update everywhere
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! initial model values
          beta0 = model_vs(i,j,k,ispec)
          rho0 = model_rho(i,j,k,ispec)
          alpha0 = model_vp(i,j,k,ispec)

          beta1 = 0._CUSTOM_REAL
          rho1 = 0._CUSTOM_REAL
          alpha1 = 0._CUSTOM_REAL

          ! isotropic model update

          ! shear values
          dbetaiso = model_dbeta(i,j,k,ispec)
          beta1 = beta0 * exp( dbetaiso )

          ! density
          rho1 = rho0 * exp( model_drho(i,j,k,ispec) )

          ! alpha values
          dbulk = model_dbulk(i,j,k,ispec)
          ! if (USE_EMPIRICAL_VP) then
          !   beta_km = beta1/1000
          !   alpha1 = 0.9409 + 2.0947 * beta_km - 0.8206 * beta_km**2+ &
          !            0.2683 * beta_km**3 - 0.0251 * beta_km**4
          !   alpha1 = alpha1 *1000
          ! else 
          if (USE_ALPHA_BETA_RHO) then
            ! new vp values use alpha model update
            alpha1 = alpha0 * exp( dbulk )
          else
            ! new vp values use bulk model update:
            ! this is based on vp_new = sqrt( bulk_new**2 + 4/3 vs_new**2)
            alpha1 = sqrt( alpha0**2 * exp(2.0*dbulk) + FOUR_THIRDS * beta0**2 * ( &
                              exp(2.0*dbetaiso) - exp(2.0*dbulk) ) )
          endif
          call alpha_scaling(alpha1, beta1)
          ! stores new model values
          model_vp_new(i,j,k,ispec) = alpha1
          model_vs_new(i,j,k,ispec) = beta1
          model_rho_new(i,j,k,ispec) = rho1

        enddo
      enddo
    enddo
  enddo
  call synchronize_all()

  ! stores new model in files
  call write_new_model_iso()

  ! stores relative model perturbations
  !call write_new_model_perturbations_iso()

  ! computes volume element associated with points, calculates kernel integral for statistics
  call compute_kernel_integral_iso()


  !!! WK deallocate before allocate next time
  if(allocated(model_vp)) deallocate(model_vp) 
  if(allocated(model_vp_new)) deallocate(model_vp_new) 
  !if(allocated(model_vpv)) deallocate(model_vpv) 
  !if(allocated(model_vph)) deallocate(model_vph) 
  if(allocated(model_vs)) deallocate(model_vs) 
  if(allocated(model_vs_new)) deallocate(model_vs_new) 
  !if(allocated(model_vsv)) deallocate(model_vsv) 
  !if(allocated(model_vsh)) deallocate(model_vsh) 
  if(allocated(model_rho)) deallocate(model_rho) 
  if(allocated(model_rho_new)) deallocate(model_rho_new) 
  !if(allocated(model_eta)) deallocate(model_eta) 

  !
  if(allocated(model_dbulk)) deallocate(model_dbulk) 
  if(allocated(model_dbulk_old)) deallocate(model_dbulk_old) 
  if(allocated(model_dbeta)) deallocate(model_dbeta) 
  if(allocated(model_dbeta_old)) deallocate(model_dbeta_old) 
  if(allocated(model_drho)) deallocate(model_drho) 
  if(allocated(model_drho_old)) deallocate(model_drho_old) 

  if(allocated(kernel_bulk)) deallocate(kernel_bulk) 
  if(allocated(kernel_bulk_old)) deallocate(kernel_bulk_old) 
  !if(allocated(kernel_alphav)) deallocate(kernel_alphav) 
  !if(allocated(kernel_alphah)) deallocate(kernel_alphah) 
  if(allocated(kernel_beta)) deallocate(kernel_beta) 
  if(allocated(kernel_beta_old)) deallocate(kernel_beta_old) 
  !if(allocated(kernel_betav)) deallocate(kernel_betav) 
  !if(allocated(kernel_betah)) deallocate(kernel_betah) 
  if(allocated(kernel_rho)) deallocate(kernel_rho) 
  if(allocated(kernel_rho_old)) deallocate(kernel_rho_old) 
  deallocate(x,y,z)

  
  ! stop all the MPI processes, and exit
  !call finalize_mpi()

end subroutine model_update_opt

subroutine alpha_scaling(alpha, beta)
  use fwat_input, only: tomo_par
  use specfem_par, only: CUSTOM_REAL
  implicit none
  
  real(kind=CUSTOM_REAL), intent(inout) :: alpha
  real(kind=CUSTOM_REAL), intent(in) :: beta

  if (alpha/beta < tomo_par%VPVS_RATIO_RANGE(1)) then
    alpha = beta * tomo_par%VPVS_RATIO_RANGE(1)
  elseif (alpha/beta > tomo_par%VPVS_RATIO_RANGE(2)) then
    alpha = beta * tomo_par%VPVS_RATIO_RANGE(2)
  endif

end subroutine alpha_scaling
