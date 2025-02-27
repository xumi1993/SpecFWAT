module preproc_fwd
  use config
  use fwat_mpi
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use imput_params, fpar => fwat_par_global
  use fk_coupling

  implicit none

  integer :: ier

  type :: PrepareFWD
    integer :: ievt, simu_opt
    contains
    procedure :: init, calc_fk_wavefield
  end type PrepareFWD

contains

  subroutine init(this, simu_opt)
    class(PrepareFWD), intent(inout) :: this
    integer, intent(in) :: simu_opt

    this%simu_opt = simu_opt

    call force_ftz()

    ! reads in parameters
    call initialize_simulation_fwat()
    
    SAVE_MESH_FILES=.false.
    PRINT_SOURCE_TIME_FUNCTION=.false.
    MOVIE_VOLUME=.false.
    local_path_backup = LOCAL_PATH
    if (simu_opt < 3) then
      SIMULATION_TYPE = 1
      if (simu_opt == 1) then
        SAVE_FORWARD = .true.
      else
        SAVE_FORWARD = .false.
      endif
      APPROXIMATE_HESS_KL=.false.
    else
      SIMULATION_TYPE = 3
      APPROXIMATE_HESS_KL=.true.
    endif

    ! copy mass matrixs for prepare_timerun()
    if (ACOUSTIC_SIMULATION) then
      allocate(rmass_acoustic_copy(NGLOB_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating array rmass_acoustic_copy'
      rmass_acoustic_copy=rmass_acoustic
    endif
    if (ELASTIC_SIMULATION) then
      allocate(rmass_copy(NGLOB_AB),stat=ier)
      allocate(rmassx_copy(NGLOB_AB),stat=ier)
      allocate(rmassy_copy(NGLOB_AB),stat=ier)
      allocate(rmassz_copy(NGLOB_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating array rmass_copy'
      rmass_copy=rmass
      rmassx_copy=rmassx
      rmassy_copy=rmassy
      rmassz_copy=rmassz
    endif

    call read_mesh_databases_fwat()

  end subroutine init

  subroutine calc_fk_wavefield(this, evtid)
    class(PrepareFWD), intent(inout) :: this
    character(len=*), intent(in) :: evtid
    integer :: iev

    ! do iev = 1, fpar%acqui%nevents
    if (.not. check_fk_files(evtid)) then
      call couple_with_injection_prepare_boundary_fwat(evtid)
    endif
    ! enddo

  end subroutine calc_fk_wavefield


end module