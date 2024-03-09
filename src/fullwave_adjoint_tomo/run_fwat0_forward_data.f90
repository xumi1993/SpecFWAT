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
subroutine run_forward_data(model,evtset,simu_type)

  use fullwave_adjoint_tomo_par
  use fwat_input
  use fwat_utils
  use my_mpi
  use specfem_par
  use specfem_interface


  character(len=MAX_STRING_LEN)                   :: model 
  character(len=MAX_STRING_LEN)                   :: evtset 
  character(len=MAX_STRING_LEN)                   :: simu_type
  character(len=MAX_STRING_LEN)                                 :: fwd_dir 
  integer                                                       :: ier
  integer                                                       :: ievt
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_copy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx_copy,rmassy_copy,rmassz_copy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_acoustic_copy

  
  !**************** Build directories for the storage ****************************
  call world_rank(myrank)  
  call read_fwat_par_file()
  call setup_common_variables(simu_type)
  if (myrank == 0) then
    open(unit=OUT_FWAT_LOG,file='output_fwat0_log_'//trim(model)//'.'//trim(evtset)//'.txt')
    write(OUT_FWAT_LOG,*) 'This is run_forward_data !!!'
    call system('mkdir -p solver')
    call system('mkdir -p fwat_data')
    call system('mkdir -p solver/'//trim(model)//'.'//trim(evtset))
  endif
  call acqui_par%read_source_set(model, evtset, simu_type)
  call fwat1_out_log(model, evtset, simu_type)
  if (myrank == 0) then
    do ievt=1,acqui_par%nevents
      fwd_dir='fwat_data/'//trim(acqui_par%evtid_names(ievt))
      call system('mkdir -p '//trim(fwd_dir))
      fwd_dir=acqui_par%out_fwd_path(ievt)
      call system('mkdir -p '//trim(fwd_dir))
      call system('mkdir -p '//trim(fwd_dir)//'/EKERNEL')
      call system('mkdir -p '//trim(fwd_dir)//'/OUTPUT_FILES')
    enddo
  endif
  !*************** end of building directories ******************************
     
  ! ************** From specfem3D() **************

  ! force Flush-To-Zero if available to avoid very slow Gradual Underflow trapping
  call force_ftz()

  ! reads in parameters
  call initialize_simulation_fwat()

  PRINT_SOURCE_TIME_FUNCTION=.true.
  SAVE_FORWARD=.false.
  SAVE_MESH_FILES=.false.

  ! reads in external mesh
  if (ADIOS_FOR_MESH) then
    call read_mesh_databases_adios()
  else
    call read_mesh_databases()
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
    call read_mesh_databases_moho()
  endif

  ! reads adjoint parameters
  call read_mesh_databases_adjoint()

  ! for coupling with external codes
  call couple_with_injection_setup()

  ! sets up reference element GLL points/weights/derivatives
  call setup_GLL_points()

  ! detects surfaces
  call detect_mesh_surfaces()

  !=================================================================================
  do ievt=1,acqui_par%nevents 
    source_fname=acqui_par%src_solution_file(ievt)
    station_fname=acqui_par%station_file(ievt)
    !!! WK: OUTPUT_FILES is where mesh files, seismograms are saved
    OUTPUT_FILES=trim(acqui_par%out_fwd_path(ievt))//'/OUTPUT_FILES'
    !!! WK: prname is where forward fields, kernels are saved
   !   LOCAL_PATH=trim(out_fwd_path(ievt))//'/EKERNEL'
    if (myrank==0) then 
      write(OUT_FWAT_LOG,*) 'ievt, source_fname, station_fname: ',&
          ievt,', ',trim(source_fname),', ',trim(station_fname)
      write(OUT_FWAT_LOG,*) 'out_fwd_path = ',trim(acqui_par%out_fwd_path(ievt))
    endif

    !******************* forward **********************************
    if(myrank==0) then
      ! write(OUT_FWAT_LOG,*) 'This is forward simulations ...'
      call write_timestamp_log(OUT_FWAT_LOG, 'This is forward simulations ...')
      flush(OUT_FWAT_LOG)
    endif
    SIMULATION_TYPE=1
    SAVE_FORWARD=.false.
    if (ATTENUATION) then
      COMPUTE_AND_STORE_STRAIN = .true.
    else
      COMPUTE_AND_STORE_STRAIN = .false.
    endif
    if (simu_type=='noise') then
      USE_FORCE_POINT_SOURCE = .true.
    else
      USE_FORCE_POINT_SOURCE = .false.
    endif
    if(simu_type=='tele' .or. simu_type=='rf') then
      COUPLE_WITH_INJECTION_TECHNIQUE=.true.
      INJECTION_TECHNIQUE_TYPE=3
      FKMODEL_FILE='src_rec/FKmodel_'//trim(acqui_par%evtid_names(ievt))
    else
      COUPLE_WITH_INJECTION_TECHNIQUE=.false.
    endif

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
  !   call prepare_timerun_fwat(LOCAL_PATH) !!! WK: all wavefields and kernels are reset
    call prepare_timerun() !!! WK: all wavefields and kernels are reset

    ! steps through time linesearch
    call iterate_time()

    ! saves last time frame and finishes kernel calculations
    call FinalizeSpecfem() 
    !************************** measure ******************************
    call run_semd2sac(ievt,simu_type)
    
    open(unit=1234, iostat=ier, file=trim(prname)//'absorb_field.bin', status='old')
    if (ier == 0) close(1234, status='delete') 

    if(myrank==0)  write(OUT_FWAT_LOG,*) '----------------------------------------'
  enddo ! end loop of ievt
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

end subroutine run_forward_data
