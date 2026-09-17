! Prepare incident FK boundary wavefields using the normal FWAT configuration.
module prepare_fk
  use specfem_par
  use specfem_par_coupling, only: Veloc_FK, Tract_FK, ipt_table
  use config, only: worldrank, worldsize, local_path_backup, local_path_fwat, compress_level, &
                    simu_type, model_name, event_index, use_gll, FWAT_PAR_FILE, SIMU_TYPE_TELE, INJECTION_FK
  use input_params, only: fpar => fwat_par_global
  use param_check, only: check_nevents
  use fk_coupling, only: compute_fk_wavefield, fk_gpu_available
  use common_lib, only: mkdir, get_dat_type
  use logger, only: log

  implicit none
  private
  public :: PrepareFK

  ! The event range follows FWAT's event_index; all events share one loaded mesh.
  type :: PrepareFK
    integer :: first_event = 1, last_event = 0
  contains
    procedure :: init
    procedure :: prepare_for_event
    procedure :: finalize
  end type PrepareFK

contains

  subroutine init(this)
    class(PrepareFK), intent(inout) :: this
    integer :: ievt
    logical :: exists

    ! MPI and command-line parsing belong to the caller.
    myrank = worldrank
    simu_type = SIMU_TYPE_TELE
    model_name = ''
    ! FK uses the existing GLL mesh and needs no regular inversion grid.
    use_gll = .true.
    call fpar%read(FWAT_PAR_FILE)
    call read_parameter_file(.true.)
    local_path_backup = LOCAL_PATH
    call fpar%select_simu_type()

    ! DT/NSTEP and compression now come from TELE; mesh paths follow FWAT's joint-inversion rules.
    if (worldsize /= NPROC) call exit_MPI(myrank, 'MPI rank count must equal NPROC in DATA/Par_file')
    if (NUMBER_OF_SIMULTANEOUS_RUNS /= 1) call exit_MPI(myrank, 'Only one MPI group is supported')
    if (HDF5_IO_NODES > 0) call exit_MPI(myrank, 'Dedicated HDF5 I/O ranks are not supported by standalone FK')
    if (fpar%sim%INJECTION_TYPE /= INJECTION_FK) call exit_MPI(myrank, 'TELE.INJECTION_TYPE must select FK')
    if (DT <= 0 .or. NSTEP < 1) call exit_MPI(myrank, 'TELE.DT and TELE.NSTEP must be positive')
    if (compress_level < 0 .or. compress_level > 9) call exit_MPI(myrank, 'TELE.COMPRESS_LEVEL must be 0..9')
    if (GPU_MODE .and. .not. fk_gpu_available) call exit_MPI(myrank, 'GPU_MODE requires a USE_CUDA build for FK')

    ! TELE_TYPE chooses sources_tele/rf/telecc.dat through the existing acquisition reader.
    call get_dat_type()
    call fpar%acqui%read()
    call check_nevents()
    this%first_event = 1
    this%last_event = fpar%acqui%nevents
    if (event_index > 0) then
      if (event_index > this%last_event) call exit_MPI(myrank, 'Event index exceeds the FWAT source list')
      this%first_event = event_index
      this%last_event = event_index
    endif

    ! Use the acquisition reader's model paths to reject missing files before the legacy FK reader runs.
    do ievt = this%first_event, this%last_event
      inquire(file=trim(fpar%acqui%fkmodel_file(ievt)), exist=exists)
      if (.not. exists) call exit_MPI(myrank, 'Missing FK model: '//trim(fpar%acqui%fkmodel_file(ievt)))
    enddo

    call read_mesh()
  end subroutine init

  subroutine read_mesh()
    logical :: use_gpu
    character(len=MAX_STRING_LEN) :: message

    ! Initialize host mesh arrays once, then restore the configured GPU_MODE for FK computation.
    LOCAL_PATH = local_path_fwat
    SIMULATION_TYPE = 1
    SAVE_FORWARD = .false.
    ANISOTROPY = fpar%update%MODEL_TYPE > 1
    use_gpu = GPU_MODE
    GPU_MODE = .false.
    t0 = 0.d0
    call mkdir(OUTPUT_FILES)
    call log%init(trim(OUTPUT_FILES)//'/output_fk.log')
    call log%write('Reading mesh databases from '//trim(LOCAL_PATH), .true.)
    call initialize_simulation_fwat()
    call read_mesh_databases_fwat()
    GPU_MODE = use_gpu
    if (POROELASTIC_SIMULATION) call exit_MPI(myrank, 'FK does not support poroelastic meshes')
    call log%write('Finished reading mesh databases', .true.)
    write(message, '(a,i0,a,es12.5,a,l1)') 'NSTEP = ', NSTEP, ', DT = ', DT, ', GPU_MODE = ', GPU_MODE
    call log%write(trim(message))
  end subroutine read_mesh

  subroutine prepare_for_event(this, ievt)
    class(PrepareFK), intent(in) :: this
    integer, intent(in) :: ievt

    if (ievt < this%first_event .or. ievt > this%last_event) &
      call exit_MPI(myrank, 'Event index is outside the selected FK event range')

    ! Reuse preproc's CPU/GPU dispatch and SAVE_FK policy for this event.
    call compute_fk_wavefield(fpar%acqui%evtid_names(ievt))
    ! No SEM time stepping follows, so release this event's injection arrays.
    deallocate(Veloc_FK, Tract_FK, ipt_table)
  end subroutine prepare_for_event

  subroutine finalize(this)
    class(PrepareFK), intent(inout) :: this

    ! Release acquisition shared memory and close project and solver logs.
    call fpar%acqui%finalize()
    call log%finalize()
    if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)
    this%first_event = 1
    this%last_event = 0
  end subroutine finalize

end module prepare_fk
