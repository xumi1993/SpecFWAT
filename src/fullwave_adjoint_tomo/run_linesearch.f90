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
!
subroutine run_linesearch(model,evtset,simu_type)

  use fullwave_adjoint_tomo_par
  use fwat_input
  use fwat_utils
  use my_mpi
  use specfem_par
  use specfem_interface
  use ma_variables, only: OUT_DIR


  character(len=MAX_STRING_LEN)                   :: model 
  character(len=MAX_STRING_LEN)                   :: evtset 
  character(len=MAX_STRING_LEN)                   :: simu_type
  character(len=MAX_STRING_LEN)                   :: old_local_path
  ! character(len=MAX_STRING_LEN)                                 :: evtnm 
  character(len=MAX_STRING_LEN)                                 :: fwd_dir 
  integer                                                       :: ier
  integer                                                       :: ievt
  logical                                         :: VERBOSE_MODE_copy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_copy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx_copy,rmassy_copy,rmassz_copy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_acoustic_copy

  
  !**************** Build directories for the storage ****************************
  call world_rank(myrank)  
  call setup_common_variables(simu_type)
  VERBOSE_MODE_copy = VERBOSE_MODE
  VERBOSE_MODE = .false. ! turn off verbose mode for line search
  if (myrank == 0) then
    open(unit=OUT_FWAT_LOG,file='output_fwat3_log_'//trim(model)//'.'//trim(evtset)//'.txt')
    write(OUT_FWAT_LOG,*) 'This is run_linesearch !!!'
    call system('mkdir -p solver')
    call system('mkdir -p misfits')
    call system('mkdir -p solver/'//trim(model)//'.'//trim(evtset))
  endif
  call acqui_par%read_source_set(model,evtset,simu_type)
  call fwat1_out_log(model, evtset,simu_type)
  if (myrank==0) then
    do ievt=1,acqui_par%nevents
      fwd_dir='solver/'//trim(model)//'.'//trim(evtset)//'/'//trim(acqui_par%evtid_names(ievt))
      call system('mkdir -p '//trim(fwd_dir))
      call system('mkdir -p '//trim(fwd_dir)//'/EKERNEL')
      call system('mkdir -p '//trim(fwd_dir)//'/OUTPUT_FILES')
    enddo
  endif
  call synchronize_all()

  !*************** end of building directories ******************************
     
  ! ************** From specfem3D() **************

  ! force Flush-To-Zero if available to avoid very slow Gradual Underflow trapping
  call force_ftz()

  ! reads in parameters
  call initialize_simulation_fwat()

  if (ATTENUATION) then
    COMPUTE_AND_STORE_STRAIN = .true.
  else
    COMPUTE_AND_STORE_STRAIN = .false.
  endif

  !!enforce to allocate arrays in read_mesh_databases_adjoint()
  APPROXIMATE_HESS_KL=.false. !! test preconditionner
  PRINT_SOURCE_TIME_FUNCTION=.false.
  MOVIE_VOLUME=.false.
  SAVE_MESH_FILES=.false.
  ! LOCAL_PATH='optimize/MODEL_'//trim(model)
  old_local_path = LOCAL_PATH
  if (myrank==0) write(*,*) 'Read database from ',trim(LOCAL_PATH)
  
  ! reads in external mesh
  if (ADIOS_FOR_MESH) then
    call read_mesh_databases_adios()
  else
    call read_mesh_databases_fwat()
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

  ! reads in moho mesh
  if (ADIOS_FOR_MESH) then
    call read_mesh_databases_moho_adios()
  else
    call read_mesh_databases_moho_fwat()
  endif

  ! reads adjoint parameters
  call read_mesh_databases_adjoint_fwat()

  ! for coupling with external codes
  call couple_with_injection_setup()
  
  ! sets up reference element GLL points/weights/derivatives
  call setup_GLL_points()

  ! detects surfaces
  call detect_mesh_surfaces()

  !=================================================================================
  do ievt=1,acqui_par%nevents 
    if (simu_type=='noise') then
      USE_FORCE_POINT_SOURCE = .true.
    else
      USE_FORCE_POINT_SOURCE = .false.
    endif
    if (simu_type=='noise' .or. simu_type=='leq') then
      source_fname=acqui_par%src_solution_file(ievt)
      COUPLE_WITH_INJECTION_TECHNIQUE=.false.
    else
      source_fname='DATA/CMTSOLUTION'
      COUPLE_WITH_INJECTION_TECHNIQUE=.true.
      INJECTION_TECHNIQUE_TYPE=3
      FKMODEL_FILE='src_rec/FKmodel_'//trim(acqui_par%evtid_names(ievt))
    endif
    station_fname=acqui_par%station_file(ievt)
    !!! WK: OUTPUT_FILES is where mesh files, seismograms are saved
    OUTPUT_FILES=trim(acqui_par%out_fwd_path(ievt))//'/OUTPUT_FILES'
    !!! WK: prname is where forward fields, kernels are saved
    LOCAL_PATH=trim(acqui_par%out_fwd_path(ievt))//'/EKERNEL'
    call create_name_database(prname,myrank,LOCAL_PATH) 
    if (myrank==0) then 
      write(OUT_FWAT_LOG,*) 'ievt, source_fname, station_fname= ',ievt,trim(source_fname),trim(station_fname)
      write(OUT_FWAT_LOG,*) 'out_fwd_path= ',trim(acqui_par%out_fwd_path(ievt))
    endif

    !******************* forward **********************************
    if(myrank==0) then
      call write_timestamp_log(OUT_FWAT_LOG, 'This is forward simulations ...')
      flush(OUT_FWAT_LOG)
    endif
    SIMULATION_TYPE=1
    SAVE_FORWARD=.false.

    ! Initialize variable for solver
    call InitSpecfem()
    ! prepares sources and receivers
    call setup_sources_receivers_fwat(source_fname,station_fname)

    ! sets up and precomputes simulation arrays

    if (ACOUSTIC_SIMULATION) then
      if( .not. allocated(rmass_acoustic)) allocate(rmass_acoustic(NGLOB_AB))
      rmass_acoustic=rmass_acoustic_copy
    endif
    if (ELASTIC_SIMULATION) then
      if( .not. allocated(rmass)) allocate(rmass(NGLOB_AB))
      if( .not. allocated(rmassx)) allocate(rmassx(NGLOB_AB))
      if( .not. allocated(rmassy)) allocate(rmassy(NGLOB_AB))
      if( .not. allocated(rmassz)) allocate(rmassz(NGLOB_AB))
      rmass=rmass_copy
      rmassx=rmassx_copy
      rmassy=rmassy_copy
      rmassz=rmassz_copy
    endif
    call prepare_timerun_fwat(old_local_path) !!! WK: all wavefields and kernels are reset

    ! steps through time linesearch
    call iterate_time()

    ! saves last time frame and finishes kernel calculations
    call FinalizeSpecfem() 
    !************************** measure ******************************
    OUT_DIR=trim(OUTPUT_FILES)
    call run_preprocessing(model,evtset,ievt,simu_type,0)

    LOCAL_PATH = old_local_path
    VERBOSE_MODE = VERBOSE_MODE_copy
    if(myrank==0) then
      write(OUT_FWAT_LOG,*) '----------------------------------------'
      flush(OUT_FWAT_LOG)
    endif
  enddo ! end loop of ievt

  ! close misfit file
  block
    integer :: nflt, iband, chi_fileid, win_fileid
    if(myrank==0) then
      if (simu_type == 'rf') then
        nflt = rf_par%NGAUSS
      else
        nflt = NUM_FILTER
      endif
      do iband=1,NUM_FILTER
        chi_fileid=30+iband;win_fileid=100+iband
        close(chi_fileid)
        close(win_fileid)
      enddo
    endif
  end block
  !=================================================================================
 
  !********************* end  ******************************************************
  if(myrank==0)  write(OUT_FWAT_LOG,*) '*******************************************************'
  if(myrank==0)  write(OUT_FWAT_LOG,*) 'Finished simulations here!!!'
  if(myrank==0) close(OUT_FWAT_LOG)
  call synchronize_all()
  if(allocated(rmass_copy)) deallocate(rmass_copy)
  if(allocated(rmassx_copy)) deallocate(rmassx_copy)
  if(allocated(rmassy_copy)) deallocate(rmassy_copy)
  if(allocated(rmassz_copy)) deallocate(rmassz_copy)
  if(allocated(rmass_acoustic_copy)) deallocate(rmass_acoustic_copy)
  call fwat_parameters_free()

end subroutine run_linesearch
