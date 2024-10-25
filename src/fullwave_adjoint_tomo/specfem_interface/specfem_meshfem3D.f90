
subroutine meshfem3d_fwat()
  use meshfem3D_par
  use chunk_earth_mod

  implicit none

  integer :: iprocnum
  integer :: iproc_xi,iproc_eta
  integer :: ier

  ! interface parameters
  integer :: ilayer
  call world_rank(myrank)
  call world_size(sizeprocs)

  ! open main output file, only written to by process 0
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) then
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_meshfem3D.txt',status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,'Error could not open output file :',trim(OUTPUT_FILES)//'/output_meshfem3D.txt'
      stop 'Error opening output file'
    endif
  endif

  ! if meshing a chunk of the Earth, call a specific internal mesher designed specifically for that
  ! CD CD change this to have also the possibility to use a chunk without coupling
  if (MESH_A_CHUNK_OF_THE_EARTH) then

    ! user output

    !! VM VM : new way to mesh and store mesh in geocubit format
    if (myrank == 0) then  !! serial mesh and use decompose_mesh after
      call mesh_chunk_earth()
      write(*,*) 'Done creating a chunk of the earth Mesh (HEX8 elements), see directory MESH/'
    endif

    ! make sure everybody is synchronized
    call synchronize_all()

    ! stop program
    call finalize_mpi()
    stop

    !call bcast_input_param_to_all()
    !call read_mesh_parameter_file()


    !! VM VM old way to create  MESH_A_CHUNK_OF_THE_EARTH, but still some routines that will
    !! move in the new way thus not remove for now
    if (NGNOD == 8) then
      ! creates mesh in MESH/
      call earth_chunk_HEX8_Mesher(NGNOD)
      ! done with mesher
      stop 'Done creating a chunk of the earth Mesh (HEX8 elements), see directory MESH/'

    else if (NGNOD == 27) then

      ! creates mesh in MESH/
      call earth_chunk_HEX27_Mesher(NGNOD)
      ! done with mesher
      stop 'Done creating a chunk of the earth Mesh (HEX27 elements), see directory MESH/'

    else
      stop 'Bad number of nodes per hexahedron: NGNOD must be equal to 8 or 27'
    endif

    ! make sure everybody is synchronized
    call synchronize_all()
  endif
  ! read the mesh parameter file (Data/meshfem3D_files/Mesh_Par_file)
  ! nullify(subregions,material_properties)
  ! call read_mesh_parameter_file()

  ! read the parameter file (DATA/Par_file)
  ! BROADCAST_AFTER_READ = .true.
  ! call read_parameter_file(myrank,BROADCAST_AFTER_READ)
  ! get interface data from external file to count the spectral elements along Z
  call get_interfaces_mesh_count_fwat()

  ! compute total number of spectral elements in vertical direction
  NER = sum(ner_layer)

  ! checks if regions and vertical layers from interfaces file match
  if (maxval(subregions(:,6)) /= NER) then
    print *,'Error invalid total number of element layers in vertical direction!'
    print *,'from interface file, total layers = ',NER
    print *,'should be equal to maximum layer NZ_END specified in regions:', maxval(subregions(:,6))
    stop 'Error invalid total number of vertical layers'
  endif

  ! checks irregular grid entries
  if (.not. USE_REGULAR_MESH) then
    if (maxval(ner_doublings(:)) == NER) then
      print *,'Error invalid doubling layer NZ index too close to surface layer ',NER
      print *,'Please decrease maximum doubling layer index NZ_DOUBLING'
      stop 'Error invalid doubling layer index too close to surface layer'
    endif
  endif


  ! compute other parameters based upon values read
  call compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
                          NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                          NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
                          NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
                          NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                          NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB, &
                          USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

  ! dynamic allocation of mesh arrays
  call initialize_mesh_arrays()
  allocate(rns(0:2*NER),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1352')
  if (ier /= 0) stop 'Error allocating array rns'

  allocate(xgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1353')
  if (ier /= 0) stop 'Error allocating array xgrid'
  allocate(ygrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1354')
  if (ier /= 0) stop 'Error allocating array ygrid'
  allocate(zgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1355')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  allocate(addressing(0:NPROC_XI-1,0:NPROC_ETA-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1356')
  if (ier /= 0) stop 'Error allocating array addressing'
  allocate(iproc_xi_slice(0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1357')
  if (ier /= 0) stop 'Error allocating array iproc_xi_slice'
  allocate(iproc_eta_slice(0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1358')
  if (ier /= 0) stop 'Error allocating array iproc_eta_slice'

  ! clear arrays
  xgrid(:,:,:) = 0.d0
  ygrid(:,:,:) = 0.d0
  zgrid(:,:,:) = 0.d0

  iproc_xi_slice(:) = 0
  iproc_eta_slice(:) = 0

  ! create global slice addressing for solver
  do iproc_eta = 0,NPROC_ETA-1
    do iproc_xi = 0,NPROC_XI-1
      iprocnum = iproc_eta * NPROC_XI + iproc_xi
      iproc_xi_slice(iprocnum) = iproc_xi
      iproc_eta_slice(iprocnum) = iproc_eta
      addressing(iproc_xi,iproc_eta) = iprocnum
    enddo
  enddo

    ! check that the constants.h file is correct
  if (NGNOD /= 8) call exit_MPI(myrank,'volume elements should have 8 control nodes in our internal mesher')
  if (NGNOD2D /= 4) call exit_MPI(myrank,'surface elements should have 4 control nodes in our internal mesher')

  ! check that reals are either 4 or 8 bytes
  if (CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) call exit_MPI(myrank,'wrong size of CUSTOM_REAL for reals')

  ! check that number of slices is at least 1 in each direction
  if (NPROC_XI < 1) call exit_MPI(myrank,'NPROC_XI must be greater than 1')
  if (NPROC_ETA < 1) call exit_MPI(myrank,'NPROC_ETA must be greater than 1')

  ! check that mesh can be cut into the right number of slices
  if (mod(NEX_XI,NPROC_XI) /= 0) call exit_MPI(myrank,'NEX_XI must be a multiple of NPROC_XI for a regular mesh')
  if (mod(NEX_ETA,NPROC_ETA) /= 0) call exit_MPI(myrank,'NEX_ETA must be a multiple of NPROC_ETA for a regular mesh')

  ! also check that mesh can be coarsened in depth twice (block size multiple of 8)
  ! i.e. check that NEX is divisible by 8 and that NEX_PER_PROC is divisible by 8
  ! This is not required for a regular mesh
  if (.not. USE_REGULAR_MESH) then
    if (mod(NEX_XI,8) /= 0) call exit_MPI(myrank,'NEX_XI must be a multiple of 8')
    if (mod(NEX_ETA,8) /= 0) call exit_MPI(myrank,'NEX_ETA must be a multiple of 8')

    if (mod(NEX_PER_PROC_XI,8) /= 0) call exit_MPI(myrank,'NEX_PER_PROC_XI must be a multiple of 8')
    if (mod(NEX_PER_PROC_ETA,8) /= 0) call exit_MPI(myrank,'NEX_PER_PROC_ETA must be a multiple of 8')

    if (mod(NEX_PER_PROC_XI, 2**NDOUBLINGS * 2) /= 0 ) &
      call exit_MPI(myrank,'NEX_PER_PROC_XI must be a multiple of 2 * 2**NDOUBLINGS')
    if (mod(NEX_PER_PROC_ETA, 2**NDOUBLINGS * 2) /= 0 ) &
      call exit_MPI(myrank,'NEX_PER_PROC_ETA must be a multiple of 2 * 2**NDOUBLINGS')
  endif

  ! get addressing for this process
  iproc_xi_current = iproc_xi_slice(myrank)
  iproc_eta_current = iproc_eta_slice(myrank)

  ! number of elements in each slice
  npx_element_steps = 2*NEX_PER_PROC_XI
  npy_element_steps = 2*NEX_PER_PROC_ETA
  ner_layer(:) = 2 * ner_layer(:)

  ! creates mesh interfaces
  call create_interfaces_mesh()

  ! creates mesh element points
  call create_meshfem_mesh()

  ! make sure everybody is synchronized
  if (myrank == 0) close(IMAIN)
  call synchronize_all()
end subroutine meshfem3d_fwat

subroutine generate_database_fwat()
  use adios_manager_mod
  use generate_databases_par

  implicit none

  ! open main output file, only written to by process 0
  call world_rank(myrank)
  call world_size(sizeprocs)
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
  open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_generate_databases.txt',status='unknown')

  ! read the parameter file
  ! call read_parameters() !! MJ: not needed, already read outside of this subroutine

  ! reads topography and bathymetry file
  call read_topography()
  !! MJ: ignore ADIOS_FOR_DATABASES for now, always use the non-ADIOS version
  ! Initialize ADIOS I/O
  ! if (ADIOS_ENABLED) then
  !   call adios_setup()
  ! endif

  !! MJ: ignore ADIOS_FOR_DATABASES for now, always use the non-ADIOS version
  ! reads Databases files
  ! if (ADIOS_FOR_DATABASES) then
  !   call read_partition_files_adios()
  ! else
  !   call read_partition_files()
  ! endif 

  call initialize_partition_arrays()
  call read_partition_files()
  ! external mesh creation
  call setup_mesh_fwat()

  ! finalize mesher
  call finalize_databases()
  ! call save_header_file(NSPEC_AB,NGLOB_AB,NPROC, &
  !                         ATTENUATION,ANISOTROPY,NSTEP,DT,STACEY_INSTEAD_OF_FREE_SURFACE, &
  !                         SIMULATION_TYPE,max_memory_size,nfaces_surface_glob_ext_mesh)

  ! if (ADIOS_ENABLED) then
  !   call adios_cleanup()
  ! endif
  if (myrank == 0) close(IMAIN)
  
end subroutine generate_database_fwat

subroutine setup_mesh_fwat()
  use generate_databases_par
  use create_regions_mesh_ext_par
  implicit none

  ! compute maximum number of points
  npointot = NSPEC_AB * NGLLX * NGLLY * NGLLZ

  ! use dynamic allocation to allocate memory for arrays
  if(allocated(ibool)) deallocate(ibool)
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 605')
  if (ier /= 0) stop 'error allocating array ibool'
  if(allocated(xstore)) deallocate(xstore)
  allocate(xstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 606')
  if (ier /= 0) stop 'error allocating array xstore'
  if(allocated(ystore)) deallocate(ystore)
  allocate(ystore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 607')
  if (ier /= 0) stop 'error allocating array ystore'
  if(allocated(zstore)) deallocate(zstore)
  allocate(zstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 608')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  call memory_eval_mesher(myrank,NSPEC_AB,npointot,nnodes_ext_mesh, &
                      nelmnts_ext_mesh,nmat_ext_mesh,num_interfaces_ext_mesh, &
                      max_interface_size_ext_mesh,nspec2D_xmin,nspec2D_xmax, &
                      nspec2D_ymin,nspec2D_ymax,nspec2D_bottom,nspec2D_top, &
                      max_memory_size_request)
  max_memory_size = max_memory_size_request

  ! make sure everybody is synchronized
  call synchronize_all()
  
  call crm_ext_deallocate_arrays()
  ! call create_regions_mesh()
  call create_regions_mesh_fwat()

  ! print min and max of topography included
  min_elevation = HUGEVAL
  max_elevation = -HUGEVAL
  do iface = 1,nspec2D_top_ext
    do icorner = 1,NGNOD2D
      inode = nodes_ibelm_top(icorner,iface)
      if (nodes_coords_ext_mesh(3,inode) < min_elevation) then
        min_elevation = nodes_coords_ext_mesh(3,inode)
      endif
      if (nodes_coords_ext_mesh(3,inode) > max_elevation) then
        max_elevation = nodes_coords_ext_mesh(3,inode)
      endif
    enddo
  enddo
end subroutine setup_mesh_fwat

subroutine initialize_mesh_arrays()
  use meshfem3D_par
  use create_meshfem_par
  implicit none

  if (allocated(rns)) deallocate(rns)
  if (allocated(xgrid)) deallocate(xgrid)
  if (allocated(ygrid)) deallocate(ygrid)
  if (allocated(zgrid)) deallocate(zgrid)
  if (allocated(addressing)) deallocate(addressing)
  if (allocated(iproc_xi_slice)) deallocate(iproc_xi_slice)
  if (allocated(iproc_eta_slice)) deallocate(iproc_eta_slice)
  if (allocated(iboun)) deallocate(iboun)
  if (allocated(iMPIcut_xi)) deallocate(iMPIcut_xi)
  if (allocated(iMPIcut_eta)) deallocate(iMPIcut_eta)
  if (allocated(ispec_material_id)) deallocate(ispec_material_id)
  
end subroutine initialize_mesh_arrays

subroutine initialize_partition_arrays()
  use generate_databases_par
  use create_regions_mesh_ext_par
  implicit none

  if (allocated(nodes_coords_ext_mesh)) deallocate(nodes_coords_ext_mesh)
  if (allocated(materials_ext_mesh)) deallocate(materials_ext_mesh)
  if (allocated(undef_mat_prop)) deallocate(undef_mat_prop)
  if (allocated(elmnts_ext_mesh)) deallocate(elmnts_ext_mesh)
  if (allocated(mat_ext_mesh)) deallocate(mat_ext_mesh)
  if (allocated(ibelm_xmin)) deallocate(ibelm_xmin)
  if (allocated(nodes_ibelm_xmin)) deallocate(nodes_ibelm_xmin)
  if (allocated(ibelm_xmax)) deallocate(ibelm_xmax)
  if (allocated(nodes_ibelm_xmax)) deallocate(nodes_ibelm_xmax)
  if (allocated(ibelm_ymin)) deallocate(ibelm_ymin)
  if (allocated(nodes_ibelm_ymin)) deallocate(nodes_ibelm_ymin)
  if (allocated(ibelm_ymax)) deallocate(ibelm_ymax)
  if (allocated(nodes_ibelm_ymax)) deallocate(nodes_ibelm_ymax)
  if (allocated(ibelm_bottom)) deallocate(ibelm_bottom)
  if (allocated(nodes_ibelm_bottom)) deallocate(nodes_ibelm_bottom)
  if (allocated(ibelm_top)) deallocate(ibelm_top)
  if (allocated(nodes_ibelm_top)) deallocate(nodes_ibelm_top)
  if (allocated(CPML_to_spec)) deallocate(CPML_to_spec)
  if (allocated(CPML_regions)) deallocate(CPML_regions)
  if (allocated(is_CPML)) deallocate(is_CPML)
  if (allocated(my_neighbors_ext_mesh)) deallocate(my_neighbors_ext_mesh)
  if (allocated(my_nelmnts_neighbors_ext_mesh)) deallocate(my_nelmnts_neighbors_ext_mesh)
  if (allocated(my_interfaces_ext_mesh)) deallocate(my_interfaces_ext_mesh)
  if (allocated(ibool_interfaces_ext_mesh)) deallocate(ibool_interfaces_ext_mesh)
  if (allocated(nibool_interfaces_ext_mesh)) deallocate(nibool_interfaces_ext_mesh)
  if (allocated(ibelm_moho)) deallocate(ibelm_moho)
  if (allocated(nodes_ibelm_moho)) deallocate(nodes_ibelm_moho)
end subroutine initialize_partition_arrays

subroutine crm_ext_deallocate_arrays()
  use generate_databases_par
  use create_regions_mesh_ext_par
  implicit none

  if (allocated(xelm)) deallocate(xelm)
  if (allocated(yelm)) deallocate(yelm)
  if (allocated(zelm)) deallocate(zelm)
  if (allocated(qkappa_attenuation_store)) deallocate(qkappa_attenuation_store)
  if (allocated(qmu_attenuation_store)) deallocate(qmu_attenuation_store)
  if (allocated(xigll)) deallocate(xigll)
  if (allocated(yigll)) deallocate(yigll)
  if (allocated(zigll)) deallocate(zigll)
  if (allocated(wxgll)) deallocate(wxgll)
  if (allocated(wygll)) deallocate(wygll)
  if (allocated(wzgll)) deallocate(wzgll)
  if (allocated(shape3D)) deallocate(shape3D)
  if (allocated(dershape3D)) deallocate(dershape3D)
  if (allocated(shape2D_x)) deallocate(shape2D_x)
  if (allocated(shape2D_y)) deallocate(shape2D_y)
  if (allocated(shape2D_bottom)) deallocate(shape2D_bottom)
  if (allocated(shape2D_top)) deallocate(shape2D_top)
  if (allocated(dershape2D_x)) deallocate(dershape2D_x)
  if (allocated(dershape2D_y)) deallocate(dershape2D_y)
  if (allocated(dershape2D_bottom)) deallocate(dershape2D_bottom)
  if (allocated(dershape2D_top)) deallocate(dershape2D_top)
  if (allocated(wgllwgll_xy)) deallocate(wgllwgll_xy)
  if (allocated(wgllwgll_xz)) deallocate(wgllwgll_xz)
  if (allocated(wgllwgll_yz)) deallocate(wgllwgll_yz)
  if (allocated(rho_vp)) deallocate(rho_vp)
  if (allocated(rho_vs)) deallocate(rho_vs)
  if (allocated(rhostore)) deallocate(rhostore)
  if (allocated(kappastore)) deallocate(kappastore)
  if (allocated(mustore)) deallocate(mustore)
  if (allocated(rhoarraystore)) deallocate(rhoarraystore)
  if (allocated(kappaarraystore)) deallocate(kappaarraystore)
  if (allocated(etastore)) deallocate(etastore)
  if (allocated(tortstore)) deallocate(tortstore)
  if (allocated(phistore)) deallocate(phistore)
  if (allocated(rho_vpI)) deallocate(rho_vpI)
  if (allocated(rho_vpII)) deallocate(rho_vpII)
  if (allocated(rho_vsI)) deallocate(rho_vsI)
  if (allocated(permstore)) deallocate(permstore)
  if (allocated(xixstore)) deallocate(xixstore)
  if (allocated(xiystore)) deallocate(xiystore)
  if (allocated(xizstore)) deallocate(xizstore)
  if (allocated(etaxstore)) deallocate(etaxstore)
  if (allocated(etaystore)) deallocate(etaystore)
  if (allocated(etazstore)) deallocate(etazstore)
  if (allocated(gammaxstore)) deallocate(gammaxstore)
  if (allocated(gammaystore)) deallocate(gammaystore)
  if (allocated(gammazstore)) deallocate(gammazstore)
  if (allocated(jacobianstore)) deallocate(jacobianstore)
  if (allocated(abs_boundary_ispec)) deallocate(abs_boundary_ispec)
  if (allocated(abs_boundary_ijk)) deallocate(abs_boundary_ijk)
  if (allocated(abs_boundary_jacobian2Dw)) deallocate(abs_boundary_jacobian2Dw)
  if (allocated(abs_boundary_normal)) deallocate(abs_boundary_normal)
  if (allocated(free_surface_ispec)) deallocate(free_surface_ispec)
  if (allocated(free_surface_ijk)) deallocate(free_surface_ijk)
  if (allocated(free_surface_jacobian2Dw)) deallocate(free_surface_jacobian2Dw)
  if (allocated(free_surface_normal)) deallocate(free_surface_normal)
  if (allocated(c11store)) deallocate(c11store)
  if (allocated(c12store)) deallocate(c12store)
  if (allocated(c13store)) deallocate(c13store)
  if (allocated(c14store)) deallocate(c14store)
  if (allocated(c15store)) deallocate(c15store)
  if (allocated(c16store)) deallocate(c16store)
  if (allocated(c22store)) deallocate(c22store)
  if (allocated(c23store)) deallocate(c23store)
  if (allocated(c24store)) deallocate(c24store)
  if (allocated(c25store)) deallocate(c25store)
  if (allocated(c26store)) deallocate(c26store)
  if (allocated(c33store)) deallocate(c33store)
  if (allocated(c34store)) deallocate(c34store)
  if (allocated(c35store)) deallocate(c35store)
  if (allocated(c36store)) deallocate(c36store)
  if (allocated(c44store)) deallocate(c44store)
  if (allocated(c45store)) deallocate(c45store)
  if (allocated(c46store)) deallocate(c46store)
  if (allocated(c55store)) deallocate(c55store)
  if (allocated(c56store)) deallocate(c56store)
  if (allocated(c66store)) deallocate(c66store)
  if (allocated(ispec_is_acoustic)) deallocate(ispec_is_acoustic)
  if (allocated(ispec_is_elastic)) deallocate(ispec_is_elastic)
  if (allocated(ispec_is_poroelastic)) deallocate(ispec_is_poroelastic)
  if (allocated(ibelm_moho_bot)) deallocate(ibelm_moho_bot)
  if (allocated(ibelm_moho_top)) deallocate(ibelm_moho_top)
  if (allocated(normal_moho_top)) deallocate(normal_moho_top)
  if (allocated(normal_moho_bot)) deallocate(normal_moho_bot)
  if (allocated(ijk_moho_bot)) deallocate(ijk_moho_bot)
  if (allocated(ijk_moho_top)) deallocate(ijk_moho_top)
  if (allocated(is_moho_top)) deallocate(is_moho_top)
  if (allocated(is_moho_bot)) deallocate(is_moho_bot)
  if (allocated(ispec_is_inner)) deallocate(ispec_is_inner)
  if (allocated(phase_ispec_inner_acoustic)) deallocate(phase_ispec_inner_acoustic)
  if (allocated(phase_ispec_inner_elastic)) deallocate(phase_ispec_inner_elastic)
  if (allocated(phase_ispec_inner_poroelastic)) deallocate(phase_ispec_inner_poroelastic)
  if (allocated(ispec_is_surface_external_mesh)) deallocate(ispec_is_surface_external_mesh)
  if (allocated(iglob_is_surface_external_mesh)) deallocate(iglob_is_surface_external_mesh)
  if (allocated(itopo_bathy)) deallocate(itopo_bathy)
  if (allocated(irregular_element_number)) deallocate(irregular_element_number)
  if (allocated(num_elem_colors_acoustic)) deallocate(num_elem_colors_acoustic)
  if (allocated(num_elem_colors_elastic)) deallocate(num_elem_colors_elastic)
  if (allocated(rmass_ocean_load)) deallocate(rmass_ocean_load)
end subroutine crm_ext_deallocate_arrays

subroutine get_interfaces_mesh_count_fwat()

  use meshfem3D_par, only: myrank,INTERFACES_FILE, &
    number_of_interfaces,max_npx_interface,max_npy_interface, &
    number_of_layers,ner_layer

  use constants, only: IMAIN,IIN,MF_IN_DATA_FILES,DONT_IGNORE_JUNK,IUTM2LONGLAT,MAX_STRING_LEN

  implicit none

  ! local parameters
  integer :: npx_interface,npy_interface
  integer :: interface_current,ilayer
  integer :: ier

  ! dummy values
  double precision :: orig_x_dummy,orig_y_dummy
  double precision :: spacing_x_dummy,spacing_y_dummy
  character(len=MAX_STRING_LEN) :: dummy_file
  logical :: SUPPRESS_UTM_PROJECTION_DUMMY

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Reading interface data from file ',trim(MF_IN_DATA_FILES)//trim(INTERFACES_FILE)
    call flush_IMAIN()
  endif

  ! opens interfaces file
  open(unit=IIN,file=trim(MF_IN_DATA_FILES)//trim(INTERFACES_FILE),status='old',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening interface file: ',trim(MF_IN_DATA_FILES)//trim(INTERFACES_FILE)
    stop 'Error opening interface file'
  endif

  ! initializes maximum count
  max_npx_interface  = -1
  max_npy_interface  = -1

  ! read number of interfaces
  call read_value_integer_mesh(IIN,DONT_IGNORE_JUNK,number_of_interfaces,'NINTERFACES', ier)
  if (ier /= 0) stop 'Error reading interface parameter for NINTERFACES'

  if (number_of_interfaces < 1) stop 'Error not enough interfaces (minimum is 1, for topography)'

  ! loop on all the interfaces
  do interface_current = 1,number_of_interfaces
    call read_interface_parameters(IIN,SUPPRESS_UTM_PROJECTION_DUMMY,dummy_file, &
                                   npx_interface,npy_interface, &
                                   orig_x_dummy,orig_y_dummy, &
                                   spacing_x_dummy,spacing_y_dummy,ier)
    if (ier /= 0) then
      print *,'Error reading interface parameters: interface ',interface_current
      stop 'Error reading interface parameters for interfaces'
    endif

    max_npx_interface = max(npx_interface,max_npx_interface)
    max_npy_interface = max(npy_interface,max_npy_interface)

    if ((max_npx_interface < 2) .or. (max_npy_interface < 2)) then
      print *,'Error interface ',interface_current,': has not enough interface points (minimum is 2x2)'
      stop 'Error not enough interface points (minimum is 2x2)'
    endif
  enddo

  if (myrank == 0) then
    write(IMAIN,*) 'maximum interface points x/y = ',max_npx_interface,max_npy_interface
    call flush_IMAIN()
  endif

  ! define number of layers
  number_of_layers = number_of_interfaces

  if (allocated(ner_layer)) deallocate(ner_layer)
  allocate(ner_layer(number_of_layers),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1283')
  if (ier /= 0) stop 'Error allocating array ner_layer'

  ! loop on all the layers
  do ilayer = 1,number_of_layers

    ! read number of spectral elements in vertical direction in this layer
    call read_value_integer_mesh(IIN,DONT_IGNORE_JUNK,ner_layer(ilayer),'NER_LAYER', ier)
    if (ier /= 0) stop 'Error reading interface parameter for NER_LAYER'

    ! checks
    if (ner_layer(ilayer) < 1) then
      print *,'Error invalid layering number ',ner_layer(ilayer),' for layer ',ilayer
      print *,'Please use a minimum element layer number of 1 in interface file'
      stop 'Error not enough spectral elements along Z in layer (minimum is 1)'
    endif
  enddo

  close(IIN)

  end subroutine get_interfaces_mesh_count_fwat