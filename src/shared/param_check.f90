module param_check
  ! Validation of fwat_params.yml against SPECFEM's Par_file, shared by every
  ! step of the workflow. Each step has one entry point that groups the checks
  ! it needs, so xfwat_fwd_measure_adj, xfwat_post_proc, xfwat_optimize and the
  ! xspecfwat driver all reject the same inconsistent setup, on the GLL mesh as
  ! well as on the regular model grid.
  !
  ! The stage routines read parameters only, so they can run as soon as both
  ! files have been read and fpar%select_simu_type() has attached fpar%sim.
  ! Checks that need the source list, the mesh databases or a misfit value are
  ! separate, and are called where that state becomes available.
  use config
  use input_params, only: fpar => fwat_par_global
  use shared_parameters, only: NPROC, NUMBER_OF_SIMULTANEOUS_RUNS, HDF5_IO_NODES, &
    ADIOS_ENABLED, HDF5_ENABLED, MESH_A_CHUNK_OF_THE_EARTH, IMODEL, ANISOTROPY, &
    ANISOTROPIC_KL, ELASTIC_SIMULATION, ACOUSTIC_SIMULATION, POROELASTIC_SIMULATION

  implicit none
  private

  ! One entry point per step of the workflow.
  public :: check_model_name, check_inv_type, check_preproc_params, &
            check_postproc_params, check_optimize_params, check_inversion_params
  ! Checks on state that only exists once a step is under way.
  public :: check_nevents, check_elastic_mesh, check_misfit

contains

!===============================================================================
! Stage entry points
!===============================================================================

  !> Forward, measurement and adjoint step (xfwat_fwd_measure_adj).
  subroutine check_preproc_params(iter)
    integer, intent(in) :: iter

    call check_common_params(iter)
    call check_data_type()
    call check_waveform_params()
  end subroutine check_preproc_params

  !> Kernel summation, smoothing and tapering step (xfwat_post_proc).
  subroutine check_postproc_params(iter)
    integer, intent(in) :: iter

    call check_common_params(iter)
    call check_data_type()
    call check_kernel_postproc()
    if (.not. use_gll) call check_model_grid()
  end subroutine check_postproc_params

  !> Model update step (xfwat_optimize), on the GLL mesh or on the grid.
  subroutine check_optimize_params(iter)
    integer, intent(in) :: iter

    call check_common_params(iter)
    call check_iteration_start(iter)
    call check_optimizer()
    if (use_gll) then
      call check_gll_model_source(iter)
    else
      call check_model_grid()
    endif
  end subroutine check_optimize_params

  !> The xspecfwat driver runs all three steps in one process, so it validates
  !! every one of them up front rather than failing halfway through the first
  !! iteration. It also owns the command line, hence the extra checks on the
  !! range of iterations and on the data type selected with -s.
  subroutine check_inversion_params(first, iteration_count, itype)
    integer, intent(in) :: first, iteration_count
    integer, intent(out) :: itype

    call check_preproc_params(first)
    call check_postproc_params(first)
    call check_optimize_params(first)
    ! Models are written up to M<first+iteration_count-1>.
    call check_iteration_bounds(first, first+iteration_count-1)
    call check_inv_type(itype)
    if (itype < 1) call abort_check('xspecfwat runs one POSTPROC.INV_TYPE at a time')
    if (INV_TYPE_NAMES(itype) /= simu_type) &
      call abort_check('-s must match the enabled POSTPROC.INV_TYPE')
  end subroutine check_inversion_params

!===============================================================================
! Parameter checks
!===============================================================================

  !> What every step depends on: a usable Mxx label, the data types the
  !! inversion combines, and the parameterization of the model.
  subroutine check_common_params(iter)
    integer, intent(in) :: iter

    call check_iteration_bounds(iter, iter)
    call check_inv_type()
    call check_model_space()
  end subroutine check_common_params

  !> Accept only the Mxx model names the model, gradient and history files are
  !! keyed on, and return the iteration index they encode.
  subroutine check_model_name(iter)
    integer, intent(out) :: iter
    integer :: ios

    if (len_trim(model_name) /= 3 .or. model_name(1:1) /= 'M' .or. &
        verify(model_name(2:3), '0123456789') /= 0) &
      call abort_check('Model must have the form M00 through M98')
    read(model_name(2:3), '(I2)', iostat=ios) iter
    if (ios /= 0) call abort_check('Invalid model index')
  end subroutine check_model_name

  !> Iteration indices have to stay inside the two-digit Mxx naming, since the
  !! next model is always written as M<iter+1>.
  subroutine check_iteration_bounds(iter_first, iter_last)
    integer, intent(in) :: iter_first, iter_last

    if (fpar%update%ITER_START < 0) call abort_check('ITER_START must be non-negative')
    if (iter_first < 0) call abort_check('Model index must be non-negative')
    if (iter_last > 98) call abort_check('Models cannot run past M98')
  end subroutine check_iteration_bounds

  !> Only the optimizer needs a gradient history, so only it requires the model
  !! to be at or after the iteration the inversion started from.
  subroutine check_iteration_start(iter)
    integer, intent(in) :: iter

    if (iter < fpar%update%ITER_START) &
      call abort_check('Model to update must be at or after ITER_START')
  end subroutine check_iteration_start

  !> Data types combined in this inversion. The GLL flow updates one type at a
  !! time; the grid flow may sum several into a joint kernel. itype returns the
  !! single enabled index, or zero for a joint inversion.
  subroutine check_inv_type(itype)
    integer, intent(out), optional :: itype
    integer :: i, ntype

    ntype = count(fpar%postproc%INV_TYPE)
    if (ntype < 1) call abort_check('At least one POSTPROC.INV_TYPE must be enabled')
    if (use_gll .and. ntype /= 1) &
      call abort_check('GLL inversion requires exactly one POSTPROC.INV_TYPE')
    if (ntype > 1) then
      do i = 1, NUM_INV_TYPE
        if (fpar%postproc%INV_TYPE(i) .and. fpar%postproc%JOINT_WEIGHT(i) <= 0.0_cr) &
          call abort_check('JOINT_WEIGHT must be positive for every enabled POSTPROC.INV_TYPE')
      enddo
    endif
    if (present(itype)) then
      itype = 0
      if (ntype == 1) then
        do i = 1, NUM_INV_TYPE
          if (fpar%postproc%INV_TYPE(i)) itype = i
        enddo
      endif
    endif
  end subroutine check_inv_type

  !> The data type selected on the command line, and its teleseismic flavour.
  subroutine check_data_type()
    if (.not. any(INV_TYPE_NAMES == simu_type)) &
      call abort_check('Unknown simulation type; -s must be '//trim(SIMU_TYPE_NOISE)//', '// &
                       trim(SIMU_TYPE_TELE)//' or '//trim(SIMU_TYPE_LEQ))
    if (simu_type == SIMU_TYPE_TELE .and. &
        (fpar%sim%TELE_TYPE < 1 .or. fpar%sim%TELE_TYPE > 3)) &
      call abort_check('TELE_TYPE must be 1 (tele), 2 (rf) or 3 (telecc)')
  end subroutine check_data_type

  !> Parameterization of the update, and the I/O layout the GLL flow needs to
  !! keep its per-rank binaries aligned from one iteration to the next.
  subroutine check_model_space()
    if (fpar%update%MODEL_TYPE /= 1 .and. fpar%update%MODEL_TYPE /= 2) &
      call abort_check('MODEL_TYPE must be 1 (vp/vs/rho) or 2 (azimuthal anisotropy)')
    if (.not. use_gll) return
    if (fpar%update%MODEL_TYPE /= 1) &
      call abort_check('GLL inversion supports MODEL_TYPE: 1 (vp/vs/rho)')
    if (NPROC /= worldsize .or. NUMBER_OF_SIMULTANEOUS_RUNS /= 1 .or. HDF5_IO_NODES /= 0) &
      call abort_check('GLL inversion requires one MPI group with exactly NPROC ranks and no separate I/O ranks')
    if (ADIOS_ENABLED .or. HDF5_ENABLED .or. MESH_A_CHUNK_OF_THE_EARTH) &
      call abort_check('GLL inversion requires binary databases and the standard internal mesher')
  end subroutine check_model_space

  !> Geometry of the regular grid the model and the kernels are projected onto.
  subroutine check_model_grid()
    if (any(fpar%grid%regular_grid_size < 2)) &
      call abort_check('REGULAR_GRID_SIZE needs at least two points along each axis')
    if (any(fpar%grid%regular_grid_interval <= 0.0_cr)) &
      call abort_check('REGULAR_GRID_INTERVAL must be positive along each axis')
    if (len_trim(fpar%update%INIT_MODEL_PATH) == 0) &
      call abort_check('INIT_MODEL_PATH is required for a grid model')
  end subroutine check_model_grid

  !> Time stepping and the period bands the misfit is measured in.
  subroutine check_waveform_params()
    integer :: i

    if (fpar%sim%NSTEP < 1 .or. fpar%sim%DT <= 0.0_cr) &
      call abort_check('NSTEP and DT must be positive')
    if (fpar%sim%NRCOMP < 1) &
      call abort_check('RCOMPS must list at least one component')
    if (simu_type == SIMU_TYPE_NOISE .and. fpar%sim%NSCOMP < 1) &
      call abort_check('SCOMPS must list at least one component')
    if (fpar%sim%NUM_FILTER < 1) call abort_check('At least one period band is required')
    if (.not. allocated(fpar%sim%SHORT_P) .or. .not. allocated(fpar%sim%LONG_P)) &
      call abort_check('SHORT_P and LONG_P are not set for this data type')
    if (size(fpar%sim%SHORT_P) < fpar%sim%NUM_FILTER .or. &
        size(fpar%sim%LONG_P) < fpar%sim%NUM_FILTER) &
      call abort_check('SHORT_P and LONG_P must cover every period band')
    do i = 1, fpar%sim%NUM_FILTER
      if (fpar%sim%SHORT_P(i) <= 0.0_cr .or. fpar%sim%LONG_P(i) <= fpar%sim%SHORT_P(i)) &
        call abort_check('Each period band needs 0 < SHORT_P < LONG_P')
    enddo
  end subroutine check_waveform_params

  !> Smoothing, tapering, preconditioning and normalization of the summed
  !! kernels. Shared by the GLL and the grid post-processing.
  subroutine check_kernel_postproc()
    if (fpar%sim%SIGMA_H <= 0.0_cr .or. fpar%sim%SIGMA_V <= 0.0_cr) &
      call abort_check('PDE smoothing requires positive SIGMA_H and SIGMA_V')
    if (min(fpar%postproc%TAPER_H_SUPPRESS, fpar%postproc%TAPER_H_BUFFER, &
            fpar%postproc%TAPER_V_SUPPRESS, fpar%postproc%TAPER_V_BUFFER) < 0.0_cr) &
      call abort_check('Taper distances must be non-negative')
    if (fpar%sim%PRECOND_TYPE < DEFAULT_PRECOND .or. fpar%sim%PRECOND_TYPE > Z_SQRT_PRECOND) &
      call abort_check('PRECOND_TYPE must be 1 (default), 2 (depth) or 3 (sqrt of depth)')
    if (fpar%postproc%NORM_TYPE < 1 .or. fpar%postproc%NORM_TYPE > 3) &
      call abort_check('NORM_TYPE must be 1, 2 or 3')
  end subroutine check_kernel_postproc

  !> Search direction, step length and the bounds the update is clipped to.
  subroutine check_optimizer()
    if (use_gll) then
      if (fpar%update%OPT_METHOD /= 1 .and. fpar%update%OPT_METHOD /= 2) &
        call abort_check('GLL inversion supports OPT_METHOD: 1 (SD) or 2 (L-BFGS)')
    else
      if (fpar%update%OPT_METHOD < 1 .or. fpar%update%OPT_METHOD > 3) &
        call abort_check('OPT_METHOD must be 1 (SD), 2 (L-BFGS) or 3 (CG)')
    endif
    if (fpar%update%LBFGS_M_STORE < 1 .or. fpar%update%MAX_SLEN <= 0.0_cr) &
      call abort_check('LBFGS_M_STORE and MAX_SLEN must be positive')
    if (fpar%update%VPVS_RATIO_RANGE(1) <= sqrt(4.0_cr/3.0_cr) .or. &
        fpar%update%VPVS_RATIO_RANGE(2) < fpar%update%VPVS_RATIO_RANGE(1)) &
      call abort_check('VPVS_RATIO_RANGE must be an increasing range above sqrt(4/3)')
    if (fpar%update%DO_LS) then
      if (fpar%update%MAX_SUB_ITER < 1 .or. fpar%update%MAX_SHRINK <= 0.0_cr .or. &
          fpar%update%MAX_SHRINK >= 1.0_cr .or. fpar%update%C1 <= 0.0_cr .or. &
          fpar%update%C1 >= 1.0_cr) &
        call abort_check('Line search needs MAX_SUB_ITER >= 1 and MAX_SHRINK, C1 inside (0,1)')
    endif
  end subroutine check_optimizer

  !> Where the GLL optimizer picks up the model it updates: an external HDF5
  !! file for the first iteration, the previous GLL binaries afterwards.
  subroutine check_gll_model_source(iter)
    integer, intent(in) :: iter

    if (IMODEL /= IMODEL_USER_EXTERNAL .and. IMODEL /= IMODEL_GLL) &
      call abort_check('GLL inversion requires MODEL = external or gll in Par_file')
    if (IMODEL == IMODEL_USER_EXTERNAL .and. len_trim(fpar%update%INIT_MODEL_PATH) == 0) &
      call abort_check('INIT_MODEL_PATH is required for an external initial model')
    if (iter > 0 .and. IMODEL /= IMODEL_GLL) &
      call abort_check('Set MODEL = gll to update a model after M00')
    if (ANISOTROPY .or. ANISOTROPIC_KL) &
      call abort_check('GLL inversion requires isotropic model databases')
  end subroutine check_gll_model_source

!===============================================================================
! Checks on state built up while a step runs
!===============================================================================

  !> Called once fpar%acqui%read() has loaded the source list.
  subroutine check_nevents()
    if (fpar%acqui%nevents < 1) &
      call abort_check('No events in the selected source list')
  end subroutine check_nevents

  !> Material flags are only known once the mesh databases have been read.
  subroutine check_elastic_mesh()
    if (.not. ELASTIC_SIMULATION .or. ACOUSTIC_SIMULATION .or. POROELASTIC_SIMULATION) &
      call abort_check('GLL inversion requires a purely elastic mesh')
  end subroutine check_elastic_mesh

  !> Guard the line-search acceptance test against a corrupted objective.
  subroutine check_misfit(misfit)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(kind=dp), intent(in) :: misfit

    if (.not. ieee_is_finite(misfit)) call abort_check('Non-finite current misfit')
  end subroutine check_misfit

  subroutine abort_check(msg)
    character(len=*), intent(in) :: msg

    call exit_MPI(worldrank, trim(msg))
  end subroutine abort_check

end module param_check
