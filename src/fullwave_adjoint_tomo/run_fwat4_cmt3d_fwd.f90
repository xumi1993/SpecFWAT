!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.0
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
!
subroutine run_cmt3d_fwd(model,evtset,simu_type)

  use fullwave_adjoint_tomo_par
  use fwat_input
  use my_mpi
  use specfem_par
  use specfem_interface
  use ma_variables, only: OUT_DIR


  character(len=MAX_STRING_LEN)                   :: model 
  character(len=MAX_STRING_LEN)                   :: evtset 
  character(len=MAX_STRING_LEN)                   :: simu_type
  character(len=MAX_STRING_LEN)                                 :: evtset_file 
  character(len=MAX_STRING_LEN)                                 :: evtnm 
  character(len=MAX_STRING_LEN)                                 :: fwd_dir 
  integer                                                       :: ier
  integer                                                       :: ievt,nevt,icmt
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_copy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx_copy,rmassy_copy,rmassz_copy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_acoustic_copy

  
  leq_par%CMT3D_INV=.true.
  call exit_mpi(myrank,'****** FWAT4 for CMT inversion is under development by Kai Wang. Please do NOT use this !!! ******')
  leq_par%CMT_PAR=["Mrr","Mtt","Mpp","Mrt","Mrp","Mtp","dep","lat","lon"]
  !**************** Build directories for the storage ****************************
  call world_rank(myrank)  
  call read_fwat_par_file()
  if (myrank == 0) then
    open(unit=OUT_FWAT_LOG,file='output_fwat4_log_'//trim(model)//'.'//trim(evtset)//'.txt')
    write(OUT_FWAT_LOG,*) 'This is run_fwat4_cmt3d_fwd !!!'
    write(OUT_FWAT_LOG,*) 'model,evtset,simu_type: ',trim(model),'.', trim(evtset),' ',trim(simu_type)
    evtset_file='src_rec/sources_'//trim(evtset)//'.dat'
    write(OUT_FWAT_LOG,*) 'evtset_file: ',trim(evtset_file)
    write(OUT_FWAT_LOG,*) 'SAVE_OUTPUT_EACH_EVENT: ',tomo_par%SAVE_OUTPUT_EACH_EVENT
    call system('mkdir -p solver')
    call system('mkdir -p misfits')
    call system('mkdir -p solver/'//trim(model)//'.'//trim(evtset))

    open(unit=400,file=trim(evtset_file),status='old',iostat=ier)
    if (ier /=0) then
      print *,'Error could not open source subset file: ',trim(evtset_file) 
      call exit_mpi(myrank,'Error opening source subset file')
    endif
    ievt=0
    do 
      read(400,*,iostat=ier) evtnm 
      if (ier /=0) exit
      fwd_dir='solver/'//trim(model)//'.'//trim(evtset)//'/'//trim(evtnm)
      call system('mkdir -p '//trim(fwd_dir))
      call system('mkdir -p '//trim(fwd_dir)//'/EKERNEL')
      call system('mkdir -p '//trim(fwd_dir)//'/OUTPUT_FILES')
      ievt=ievt+1
    enddo
    close(400)
    nevt=ievt
    allocate(acqui_par%station_file(nevt))
    allocate(acqui_par%src_solution_file(nevt))
    allocate(acqui_par%evtid_names(nevt))
    allocate(acqui_par%out_fwd_path(nevt))
    allocate(acqui_par%in_dat_path(nevt))

    open(unit=400,file=trim(evtset_file),status='old',iostat=ier)
    do ievt=1,nevt
      read(400,*,iostat=ier) evtnm
      if (ier /=0) exit
      acqui_par%src_solution_file(ievt)='src_rec/CMTSOLUTION_'//trim(evtnm)
      acqui_par%evtid_names(ievt)=trim(evtnm)
      acqui_par%station_file(ievt)='src_rec/STATIONS_'//trim(evtnm)
      acqui_par%out_fwd_path(ievt)='solver/'//trim(model)//'.'//trim(evtset)//'/'//trim(evtnm)
      acqui_par%in_dat_path(ievt)='fwat_data/'//trim(evtnm)
    enddo
    write(OUT_FWAT_LOG,*) 'nevt= ',nevt
    write(OUT_FWAT_LOG,*) '*******************************************************'
    close(400)
  endif
  call bcast_all_singlei(nevt)
  if (myrank /= 0) then
     allocate(acqui_par%station_file(nevt))
     allocate(acqui_par%src_solution_file(nevt))
     allocate(acqui_par%evtid_names(nevt))
     allocate(acqui_par%out_fwd_path(nevt))
     allocate(acqui_par%in_dat_path(nevt))
  endif
  call bcast_all_ch_array(acqui_par%src_solution_file,nevt,MAX_STRING_LEN)
  call bcast_all_ch_array(acqui_par%evtid_names,nevt,MAX_STRING_LEN)
  call bcast_all_ch_array(acqui_par%station_file,nevt,MAX_STRING_LEN)
  call bcast_all_ch_array(acqui_par%out_fwd_path,nevt,MAX_STRING_LEN)
  call bcast_all_ch_array(acqui_par%in_dat_path,nevt,MAX_STRING_LEN)
  call bcast_all_ch_array(model,1,MAX_STRING_LEN)
  !*************** end of building directories ******************************
     
  ! ************** From specfem3D() **************

  ! force Flush-To-Zero if available to avoid very slow Gradual Underflow trapping
  call force_ftz()

  ! reads in parameters
  call initialize_simulation()

  ! reset this value for SIMULATION_TYPE=3
  COMPUTE_AND_STORE_STRAIN = .true.
  NSPEC_STRAIN_ONLY = NSPEC_AB
  NGLOB_ADJOINT = NGLOB_AB
  NSPEC_ADJOINT = NSPEC_AB

  !!enforce to allocate arrays in read_mesh_databases_adjoint()
  SIMULATION_TYPE=3
  APPROXIMATE_HESS_KL=.true. !! test preconditionner
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

  ! sets up reference element GLL points/weights/derivatives
  call setup_GLL_points()

  ! detects surfaces
  call detect_mesh_surfaces()

  !=================================================================================
  do ievt=1,nevt 
   do icmt=0,9
     if ( icmt==0) then
       source_fname=acqui_par%src_solution_file(ievt)
     else
       source_fname=trim(acqui_par%src_solution_file(ievt))//'_der/CMTSOLUTION_'//trim(leq_par%CMT_PAR(icmt))
     endif
     station_fname=acqui_par%station_file(ievt)
     !!! WK: OUTPUT_FILES is where mesh files, seismograms are saved
     OUTPUT_FILES=trim(acqui_par%out_fwd_path(ievt))//'/OUTPUT_FILES'
     !!! WK: prname is where forward fields, kernels are saved
     LOCAL_PATH=trim(acqui_par%out_fwd_path(ievt))//'/EKERNEL'
     call create_name_database(prname,myrank,LOCAL_PATH) 
     if (myrank==0) then 
        write(OUT_FWAT_LOG,*) 'ievt,icmt, source_fname, station_fname= ',ievt,icmt,trim(source_fname),trim(station_fname)
        write(OUT_FWAT_LOG,*) 'out_fwd_path= ',ievt,trim(acqui_par%out_fwd_path(ievt))
     endif

     !******************* forward **********************************
     if(myrank==0) write(OUT_FWAT_LOG,*) 'This is forward simulations ...'
     SIMULATION_TYPE=1
     SAVE_FORWARD=.false.
     COMPUTE_AND_STORE_STRAIN = .false.
     if(simu_type=='tele') then
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
     call prepare_timerun() !!! WK: all wavefields and kernels are reset

     ! steps through time linesearch
     call iterate_time()

     ! saves last time frame and finishes kernel calculations
     call FinalizeSpecfem() 
     !************************** measure ******************************
     OUT_DIR=trim(OUTPUT_FILES)
     call run_preprocessing(model,evtset,ievt,simu_type,icmt)
    
     if(myrank==0)  write(OUT_FWAT_LOG,*) '----------------------------------------'
   enddo ! end loop of icmt
   open(unit=1234, iostat=ier, file=trim(prname)//'absorb_field.bin', status='old')
   if (ier == 0) close(1234, status='delete') 
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

end subroutine run_cmt3d_fwd
