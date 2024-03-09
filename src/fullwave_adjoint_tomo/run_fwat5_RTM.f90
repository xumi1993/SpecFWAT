subroutine run_fwat5_RTM(evtsetb,evtsete,is_smooth)
  use fullwave_adjoint_tomo_par
  use fwat_input
  use tomography_par, only: USE_ALPHA_BETA_RHO,USE_ISO_KERNELS
  !!! for reading Database
  use specfem_par
  use specfem_par_elastic, only: ELASTIC_SIMULATION,ispec_is_elastic,rho_vp,rho_vs,min_resolved_period
  use specfem_par_acoustic, only: ACOUSTIC_SIMULATION,ispec_is_acoustic
  use specfem_par_poroelastic, only: POROELASTIC_SIMULATION,ispec_is_poroelastic,rho_vpI,rho_vpII,rho_vsI, &
  phistore,tortstore,rhoarraystore
  !!! for model update
  use tomography_model_iso
  use tomography_kernels_iso
  use taper3d
  use constants, only: IMAIN
  implicit none

  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob


  character(len=MAX_STRING_LEN)                   :: model 
  character(len=MAX_STRING_LEN)                   :: evtsetb,evtsete,is_smooth 
  character(len=MAX_STRING_LEN)                   :: strset
  character(len=MAX_STRING_LEN) :: output_dir,ekernel_dir_list
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  integer iset, setb, sete,ier
  !real :: t1, t2
  logical :: USE_GPU
  type(taper_cls) :: taper_3d

  model = 'M00'
  USE_GPU=.false.
  read(evtsetb(4:),'(I3)') setb 
  read(evtsete(4:),'(I3)') sete
  output_dir='RTM_kernels'
  ekernel_dir_list = trim(output_dir)//'/ekernel_dir.lst'
  tomo_par%KERNEL_TAPER = .false.
  ! output_dir='RTM_kernels'//trim(model)
  !**************** Build directories for the storage ****************************
  call world_rank(myrank)
  call read_fwat_par_file()
  if (myrank == 0) then
  open(unit=OUT_FWAT_LOG,file='output_fwat2_log_'//trim(model)//'.'//trim(evtsetb)&
                              //'-'//trim(evtsete)//'.txt')
  write(OUT_FWAT_LOG,*) 'Running reverse time migration !!!'
  write(OUT_FWAT_LOG,*) 'model,evtsetb,evtsete: ',trim(model),' ', trim(evtsetb),'-',trim(evtsete)
  write(OUT_FWAT_LOG,*) 'Is smooth: ', is_smooth
  write(OUT_FWAT_LOG,*) 'USE_SPH_SMOOTH: ',tomo_par%USE_SPH_SMOOTH
  write(OUT_FWAT_LOG,*) 'SIGMA_H, SIGMA_V: ',tomo_par%Tele_SIGMA_H, tomo_par%Tele_SIGMA_H
  write(OUT_FWAT_LOG,*) 'OPT_METHOD: ',trim(tomo_par%OPT_METHOD) 
  call system('mkdir -p RTM_kernels')
  open(unit=400,file=trim(ekernel_dir_list),iostat=ier)
    do iset=setb,sete
        write(strset,'(I0)') iset
        write(400,'(a)')'solver/'//trim(model)//'.set'//trim(strset)//'/GRADIENT'
    enddo
    close(400)
    write(OUT_FWAT_LOG,*) '*******************************************************'
  endif 
  call synchronize_all()
  ! -----------------------------------------------------------------
  ! precond and sum 
  if(myrank==0) write(OUT_FWAT_LOG,*) 'This is sum all kernels ...'
  call sum_preconditioned_kernels_fwat(ekernel_dir_list,output_dir)
  ! -----------------------------------------------------------------
  if (trim(is_smooth) == 'true') then
    !smooth
    ! Kai: Maybe we should move reding Database back to the smoothing subroutine.
    !!!  Read Database before smoothing !!! 
    ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some array sizes below
    call read_mesh_for_init()

    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),irregular_element_number(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 980')

    if (NSPEC_IRREGULAR > 0) then
      ! allocate arrays for storing the databases
      allocate(xix(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 981')
      allocate(xiy(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 982')
      allocate(xiz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 983')
      allocate(etax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 984')
      allocate(etay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 985')
      allocate(etaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 986')
      allocate(gammax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 987')
      allocate(gammay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 988')
      allocate(gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 989')
      allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 990')
    else
        ! allocate arrays for storing the databases
      allocate(xix(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 991')
      allocate(xiy(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 992')
      allocate(xiz(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 993')
      allocate(etax(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 994')
      allocate(etay(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 995')
      allocate(etaz(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 996')
      allocate(gammax(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 997')
      allocate(gammay(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 998')
      allocate(gammaz(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 999')
      allocate(jacobian(1,1,1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1000')
    endif

    ! mesh node locations
    allocate(xstore(NGLOB_AB), &
            ystore(NGLOB_AB), &
            zstore(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays for mesh nodes'

    ! material properties
    allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
            mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays for material properties'

    ! material flags
    allocate(ispec_is_acoustic(NSPEC_AB), &
            ispec_is_elastic(NSPEC_AB), &
            ispec_is_poroelastic(NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays for material flags'
    ispec_is_acoustic(:) = .false.
    ispec_is_elastic(:) = .false.
    ispec_is_poroelastic(:) = .false.

    ! reads in external mesh
    call read_mesh_databases()

    ! gets mesh dimensions
    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                              x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                              elemsize_min_glob,elemsize_max_glob, &
                              distance_min_glob,distance_max_glob)

    ! outputs infos
    if (myrank == 0) then
      print *,'mesh dimensions:'
      print *,'  Xmin and Xmax of the model = ',x_min_glob,x_max_glob
      print *,'  Ymin and Ymax of the model = ',y_min_glob,y_max_glob
      print *,'  Zmin and Zmax of the model = ',z_min_glob,z_max_glob
      print *
      print *,'  Max GLL point distance = ',distance_max_glob
      print *,'  Min GLL point distance = ',distance_min_glob
      print *,'  Max/min ratio = ',distance_max_glob/distance_min_glob
      print *
      print *,'  Max element size = ',elemsize_max_glob
      print *,'  Min element size = ',elemsize_min_glob
      print *,'  Max/min ratio = ',elemsize_max_glob/elemsize_min_glob
      print *
    endif
    if (myrank == 0 .and. IMAIN /= 6) &
      open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_mesh_resolution.txt',status='unknown')

    if (ELASTIC_SIMULATION) then
      call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                                ibool,xstore,ystore,zstore, &
                                kappastore,mustore,rho_vp,rho_vs, &
                                DT,model_speed_max,min_resolved_period, &
                                LOCAL_PATH,SAVE_MESH_FILES)

    else if (POROELASTIC_SIMULATION) then
      allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      rho_vp = 0.0_CUSTOM_REAL
      rho_vs = 0.0_CUSTOM_REAL
      call check_mesh_resolution_poro(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                      DT,model_speed_max,min_resolved_period, &
                                      phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                                      LOCAL_PATH,SAVE_MESH_FILES)
      deallocate(rho_vp,rho_vs)
    else if (ACOUSTIC_SIMULATION) then
      allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating array rho_vp'
      allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating array rho_vs'
      rho_vp = sqrt( kappastore / rhostore ) * rhostore
      rho_vs = 0.0_CUSTOM_REAL
      call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                                ibool,xstore,ystore,zstore, &
                                kappastore,mustore,rho_vp,rho_vs, &
                                DT,model_speed_max,min_resolved_period, &
                                LOCAL_PATH,SAVE_MESH_FILES)
      deallocate(rho_vp,rho_vs)
    endif
    kernel_names_comma_delimited='alpha_kernel'
    call smooth_sem_fwat(tomo_par%tele_sigma_h,tomo_par%tele_sigma_v,kernel_names_comma_delimited,output_dir,&
                  output_dir,USE_GPU, taper_3d)
    kernel_names_comma_delimited='beta_kernel'
    call smooth_sem_fwat(tomo_par%tele_sigma_h,tomo_par%tele_sigma_v,kernel_names_comma_delimited,output_dir,&
                  output_dir,USE_GPU, taper_3d)
    kernel_names_comma_delimited='rhop_kernel'
    call smooth_sem_fwat(tomo_par%tele_sigma_h,tomo_par%tele_sigma_v,kernel_names_comma_delimited,output_dir,&
                  output_dir,USE_GPU, taper_3d)
  endif
  open(unit=1234, iostat=ier, file=trim(ekernel_dir_list), status='old')
  if (ier == 0) close(1234, status='delete') 
end subroutine run_fwat5_RTM