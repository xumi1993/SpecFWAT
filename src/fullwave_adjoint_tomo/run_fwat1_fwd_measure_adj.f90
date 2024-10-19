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
subroutine run_fwat1_fwd_measure_adj(model,evtset,simu_type,run_opt_num)

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
  character(len=MAX_STRING_LEN)                   :: run_opt 
  integer                                         :: run_opt_num
  ! character(len=MAX_STRING_LEN)                                 :: evtset_file 
  ! character(len=MAX_STRING_LEN)                                 :: evtnm
  character(len=MAX_STRING_LEN)                                 :: fwd_dir, old_local_path
  integer                                                       :: ier
  integer                                                       :: ievt
  logical                                          :: use_attenuation
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_copy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx_copy,rmassy_copy,rmassz_copy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_acoustic_copy
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sum_rho_kl 
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: sum_cijkl_kl 
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sum_mu_kl 
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sum_kappa_kl 
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sum_hess_kl 
  integer :: imod_current,iband,chi_fileid,win_fileid

  
  !**************** Build directories for the storage ****************************
  call world_rank(myrank)  
  call setup_common_variables(simu_type)
  if (myrank == 0) then
    ! read(run_opt,*) run_opt_num
    open(unit=OUT_FWAT_LOG,file='output_fwat1_log_'//trim(model)//'.'//trim(evtset)//'.txt')
    write(OUT_FWAT_LOG,*) 'This is run_fwat1_fwd_measure_adj !!!'
    call system('mkdir -p solver')
    call system('mkdir -p misfits')
    call system('mkdir -p solver/'//trim(model)//'.'//trim(evtset))
    call system('mkdir -p solver/'//trim(model)//'.'//trim(evtset)//'/GRADIENT')
  endif
  call acqui_par%read_source_set(model,evtset,simu_type)
  call fwat1_out_log(model,evtset,simu_type)
  if (myrank==0) then
    do ievt=1,acqui_par%nevents
      fwd_dir=acqui_par%out_fwd_path(ievt)
      call system('mkdir -p '//trim(fwd_dir))
      call system('mkdir -p '//trim(fwd_dir)//'/EKERNEL')
      call system('mkdir -p '//trim(fwd_dir)//'/OUTPUT_FILES')
      call system('mkdir -p '//trim(fwd_dir)//'/SEM')
    enddo
  endif
  call synchronize_all()
  ! call bcast_all_ch_array(model,1,MAX_STRING_LEN)
  ! call bcast_all_singlei(run_opt_num)

  !*************** end of building directories ******************************
     
  ! ************** From specfem3D() **************

  ! force Flush-To-Zero if available to avoid very slow Gradual Underflow trapping
  call force_ftz()

  ! reads in parameters
  call initialize_simulation_fwat()

  ! reset this value for SIMULATION_TYPE=3
  COMPUTE_AND_STORE_STRAIN = .true.
  NSPEC_STRAIN_ONLY = NSPEC_AB
  NGLOB_ADJOINT = NGLOB_AB
  NSPEC_ADJOINT = NSPEC_AB
  use_attenuation = ATTENUATION

  !!enforce to allocate arrays in read_mesh_databases_adjoint()
  SIMULATION_TYPE=3
  APPROXIMATE_HESS_KL=.true. !! test preconditionner
  PRINT_SOURCE_TIME_FUNCTION=.true.
  SAVE_FORWARD=.false.
  MOVIE_VOLUME=.false.
  SAVE_MESH_FILES=.false.
  ! read(model(2:),'(I2.2)') imod_current
  ! if (tomo_par%DO_LS .and. imod_current>0) then
    ! LOCAL_PATH='optimize/MODEL_'//trim(model)
  ! endif
  old_local_path = LOCAL_PATH
  if(myrank==0) write(OUT_FWAT_LOG,*) 'Read database from ',trim(LOCAL_PATH)
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

  ! allocate arrays
  if (ELASTIC_SIMULATION) then
    allocate(sum_rho_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    sum_rho_kl = 0.
    if (ier /= 0) stop 'Error allocating array sum_rho_kl'
    if (ANISOTROPIC_KL) then
      allocate(sum_cijkl_kl(21,NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      sum_cijkl_kl = 0.
      if (ier /= 0) stop 'Error allocating array sum_cijkl_kl'
    else
      allocate(sum_mu_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      sum_mu_kl = 0.
      if (ier /= 0) stop 'Error allocating array sum_mu_kl'
      allocate(sum_kappa_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      sum_kappa_kl = 0.
      if (ier /= 0) stop 'Error allocating array sum_kappa_kl'
    endif
    if (APPROXIMATE_HESS_KL) then
      allocate(sum_hess_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      sum_hess_kl = 0.
      if (ier /= 0) stop 'Error allocating array sum_hess_kl'
     endif
  endif
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
      write(OUT_FWAT_LOG,*) 'ievt, source_fname, station_fname: ',ievt,trim(source_fname),trim(station_fname)
      write(*,*) 'ievt, source_fname, station_fname = ',ievt,trim(source_fname),trim(station_fname)
      write(OUT_FWAT_LOG,*) 'Weight=',src_weight(ievt)
      write(OUT_FWAT_LOG,*) 'out_fwd_path = ',trim(acqui_par%out_fwd_path(ievt))
    endif

    !******************* forward **********************************
    if(myrank==0) then
      ! write(OUT_FWAT_LOG,*) 'This is forward simulations ...'
      call write_timestamp_log(OUT_FWAT_LOG, 'This is forward simulations ...')
      flush(OUT_FWAT_LOG)
    endif
    SIMULATION_TYPE=1
    ! save forward for adjoint simulation
    if (run_opt_num.eq.3) then
      SAVE_FORWARD=.true.
    else
      SAVE_FORWARD=.false.
    endif

    ! save STRAIN when use attenuation
    if (use_attenuation) then
      ATTENUATION = .true.
      COMPUTE_AND_STORE_STRAIN = .true.
    else
      ATTENUATION = .false.
      COMPUTE_AND_STORE_STRAIN = .false.
    endif

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

!   call prepare_timerun() !!! WK: all wavefields and kernels are reset
    call prepare_timerun_fwat(old_local_path)

    ! steps through time iterations
    call iterate_time()

    ! saves last time frame and finishes kernel calculations
    call FinalizeSpecfem() 

    !************************** measure ******************************
    if (run_opt_num.gt.1) then
      OUT_DIR=trim(OUTPUT_FILES)
      call run_preprocessing(model,evtset,ievt,simu_type,0)
    endif

    !************************** adjoint ******************************
    if (run_opt_num.gt.2) then
      if(myrank==0) then
        ! write(OUT_FWAT_LOG,*) 'This is adjoint simulations ...'
        call write_timestamp_log(OUT_FWAT_LOG, 'This is adjoint simulations ...')
        flush(OUT_FWAT_LOG)
      endif
      SIMULATION_TYPE=3
      SAVE_FORWARD=.false.
      COMPUTE_AND_STORE_STRAIN=.true.
      APPROXIMATE_HESS_KL=.true.
      COUPLE_WITH_INJECTION_TECHNIQUE = .false.
      ATTENUATION = .false.

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
      call iterate_time()
      call FinalizeSpecfem()
  
      ! clean saved abs and forward fields
      !!! WK: this is a tough way to delete adundant files, it does NOT work on 
      !       some CPUs, need to be fixed in the future. search 'system' in other
      !       files.
      !call system('rm '//trim(prname)//'save_forward_arrays.bin') 
      !call system('rm '//trim(prname)//'absorb_field.bin') 
      open(unit=1234, iostat=ier, file=trim(prname)//'save_forward_arrays.bin', status='old')
      if (ier == 0) close(1234, status='delete') 
      open(unit=1234, iostat=ier, file=trim(prname)//'absorb_field.bin', status='old')
      if (ier == 0) close(1234, status='delete') 

      !********************** Set preconditioner as hessian **************************
      call run_precond()
      !********************** summation of event kernels **************************
      if (ELASTIC_SIMULATION) then
        sum_rho_kl(:,:,:,:)=sum_rho_kl(:,:,:,:)+rho_kl(:,:,:,:) 
        if (ANISOTROPIC_KL) then
          sum_cijkl_kl(:,:,:,:,:)=sum_cijkl_kl(:,:,:,:,:)+cijkl_kl(:,:,:,:,:)      
        else
          sum_mu_kl(:,:,:,:)=sum_mu_kl(:,:,:,:)+mu_kl(:,:,:,:)      
          sum_kappa_kl(:,:,:,:)=sum_kappa_kl(:,:,:,:)+kappa_kl(:,:,:,:)      
        endif
        if (APPROXIMATE_HESS_KL) then
          sum_hess_kl(:,:,:,:)=sum_hess_kl(:,:,:,:)+hess_kl(:,:,:,:)      
        endif
      endif
    endif ! end of run_opt
    if(myrank==0)  write(OUT_FWAT_LOG,*) '----------------------------------------'
  enddo ! end loop of ievt
  do iband=1,NUM_FILTER
    if(myrank==0) then
      chi_fileid=30+iband;win_fileid=100+iband
      close(chi_fileid)
      close(win_fileid)
    endif
  enddo
  !=================================================================================
  if (run_opt_num.gt.2) then
    !!! Save summed kernels into file
    if (ELASTIC_SIMULATION) then
      rho_kl(:,:,:,:)=sum_rho_kl(:,:,:,:)
      if (ANISOTROPIC_KL) then
        cijkl_kl(:,:,:,:,:)=sum_cijkl_kl(:,:,:,:,:)     
      else
      mu_kl(:,:,:,:)=sum_mu_kl(:,:,:,:)    
      kappa_kl(:,:,:,:)=sum_kappa_kl(:,:,:,:)     
      endif
      if (APPROXIMATE_HESS_KL) then
        hess_kl(:,:,:,:)=sum_hess_kl(:,:,:,:)     
      endif
    endif
    if (SIMULATION_TYPE == 3) then
      LOCAL_PATH=trim(acqui_par%out_fwd_path(1))//'/../GRADIENT'
      call create_name_database(prname,myrank,LOCAL_PATH) 
      if (myrank==0) write(OUT_FWAT_LOG,*) 'Save gradient to '//trim(prname)
      call save_adjoint_kernels_fwat() 
    endif
  endif ! end if run_opt
  LOCAL_PATH = old_local_path
  !**************************** Post-processing  *********************************** 
  !********************* end  ******************************************************
  if(myrank==0)  write(OUT_FWAT_LOG,*) '*******************************************************'
  ! if(myrank==0)  write(OUT_FWAT_LOG,*) 'Finished simulations here!!!'
  if(myrank==0) call write_timestamp_log(OUT_FWAT_LOG, 'Finished simulations here!!!') 
  if(myrank==0) close(OUT_FWAT_LOG)
  call synchronize_all()
  if(allocated(rmass_copy)) deallocate(rmass_copy)
  if(allocated(rmassx_copy)) deallocate(rmassx_copy)
  if(allocated(rmassy_copy)) deallocate(rmassy_copy)
  if(allocated(rmassz_copy)) deallocate(rmassz_copy)
  if(allocated(rmass_acoustic_copy)) deallocate(rmass_acoustic_copy)
  if(allocated(sum_rho_kl)) deallocate(sum_rho_kl)
  if(allocated(sum_cijkl_kl)) deallocate(sum_cijkl_kl)
  if(allocated(sum_mu_kl)) deallocate(sum_mu_kl)
  if(allocated(sum_kappa_kl)) deallocate(sum_kappa_kl)
  if(allocated(sum_hess_kl)) deallocate(sum_hess_kl)
  call fwat_parameters_free()

end subroutine run_fwat1_fwd_measure_adj
