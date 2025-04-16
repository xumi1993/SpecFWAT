module generate_databases_subs
  use manager_adios
  use generate_databases_par

  implicit none
  
  include 'version.fh'

contains
  subroutine generate_databases_fwat(is_get_model)
    logical, intent(in) :: is_get_model

    call world_rank(myrank)
    call world_size(sizeprocs)
    if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_generate_databases.txt',status='unknown')
    ! reads topography and bathymetry file
    call read_topography()

    ! Initialize ADIOS I/O
    if (ADIOS_ENABLED) then
      call initialize_adios()
    endif

    ! reads Databases files
    if (ADIOS_FOR_DATABASES) then
      ! call read_partition_files_adios()
      call exit_MPI(0, 'FWAT: ADIOS_FOR_DATABASES not implemented yet')
    else if (HDF5_ENABLED) then
      ! call read_partition_files_hdf5()
      call exit_MPI(0, 'FWAT: HDF5_ENABLED not implemented yet')
    else
      call initialize_partition_arrays()
      call read_partition_files()
    endif

    call setup_mesh_fwat(is_get_model)

    ! finalize mesher
    call finalize_databases()
  
    if (myrank == 0) close(IMAIN)

  end subroutine generate_databases_fwat

  subroutine setup_mesh_fwat(is_get_model)

  ! local parameters
  integer :: iface,icorner,inode,ier
  logical, intent(in) :: is_get_model

  ! compute maximum number of points
  npointot = NSPEC_AB * NGLLX * NGLLY * NGLLZ

  ! use dynamic allocation to allocate memory for arrays
  ! local to global indices array
  if(allocated(ibool)) deallocate(ibool)
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 605')
  if (ier /= 0) stop 'error allocating array ibool'
  ibool(:,:,:,:) = 0

  ! node coordinates defined on local level
  if(allocated(xstore)) deallocate(xstore)
  allocate(xstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 606')
  if (ier /= 0) stop 'error allocating array xstore'
  xstore(:,:,:,:) = 0.d0
  if(allocated(ystore)) deallocate(ystore)
  allocate(ystore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 607')
  if (ier /= 0) stop 'error allocating array ystore'
  ystore(:,:,:,:) = 0.d0
  if(allocated(zstore)) deallocate(zstore)
  allocate(zstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 608')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  zstore(:,:,:,:) = 0.d0

  ! estimates memory requirement
  call memory_eval_mesher(NSPEC_AB,npointot,nnodes_ext_mesh, &
                          nelmnts_ext_mesh,nmat_ext_mesh,num_interfaces_ext_mesh, &
                          max_interface_size_ext_mesh,nspec2D_xmin,nspec2D_xmax, &
                          nspec2D_ymin,nspec2D_ymax,nspec2D_bottom,nspec2D_top, &
                          max_memory_size_request)
  max_memory_size = max_memory_size_request

  ! make sure everybody is synchronized
  call synchronize_all()

  call crm_ext_deallocate_arrays()
  if (is_get_model) then
    call create_regions_mesh()
  else
    call create_regions_mesh_coord()
  endif

  end subroutine setup_mesh_fwat

  subroutine create_regions_mesh_coord()
    use constants, only: myrank,PARALLEL_FAULT,IMAIN

    use shared_parameters, only: &
      ACOUSTIC_SIMULATION, ELASTIC_SIMULATION, POROELASTIC_SIMULATION, &
      STACEY_ABSORBING_CONDITIONS,SAVE_MESH_FILES,PML_CONDITIONS, &
      ANISOTROPY,APPROXIMATE_OCEAN_LOAD,OLSEN_ATTENUATION_RATIO, &
      ATTENUATION,USE_OLSEN_ATTENUATION, &
      SAVE_MOHO_MESH,ATTENUATION_f0_REFERENCE, &
      LOCAL_PATH,LTS_MODE

    use generate_databases_par, only: &
        nnodes_ext_mesh,nelmnts_ext_mesh, &
        nodes_coords_ext_mesh, elmnts_ext_mesh, &
        max_memory_size,num_interfaces_ext_mesh, &
        nibool_interfaces_ext_mesh, &
        nspec2D_xmin, nspec2D_xmax, &
        nspec2D_ymin, nspec2D_ymax, &
        NSPEC2D_BOTTOM, NSPEC2D_TOP, &
        ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top, &
        nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
        nodes_ibelm_bottom,nodes_ibelm_top, &
        nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho

    ! global index array
    use generate_databases_par, only: nspec => NSPEC_AB, nglob => NGLOB_AB, ibool, xstore, ystore, zstore

    use create_regions_mesh_ext_par

    use fault_generate_databases, only: fault_read_input,fault_setup, &
                            fault_save_arrays,fault_save_arrays_test, &
                            nnodes_coords_open,nodes_coords_open,ANY_FAULT_IN_THIS_PROC, &
                            ANY_FAULT

    !! setup wavefield discontinuity interface
    use shared_parameters, only: IS_WAVEFIELD_DISCONTINUITY
    use wavefield_discontinuity_db, only: &
                                setup_boundary_wavefield_discontinuity, &
                                read_partition_files_wavefield_discontinuity

    implicit none

    ! local parameters
    ! memory size needed by the solver
    double precision :: memory_size
    real(kind=CUSTOM_REAL) :: model_speed_max,min_resolved_period
    ! timing
    double precision, external :: wtime
    double precision :: time_start
    logical, parameter :: DEBUG_TIMING = .false.

    ! get MPI starting time
    if (DEBUG_TIMING) time_start = wtime()

    ! initializes arrays
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...allocating arrays '
      call flush_IMAIN()
    endif
    call crm_ext_allocate_arrays(nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                nspec2D_bottom,nspec2D_top, &
                                nodes_coords_ext_mesh,nnodes_ext_mesh,elmnts_ext_mesh,nelmnts_ext_mesh)

    ! if faults exist this reads nodes_coords_open
    call fault_read_input(prname)

    ! user output
    ! call print_timing()

    ! fills location and weights for Gauss-Lobatto-Legendre points, shape and derivations,
    ! returns jacobianstore,xixstore,...gammazstore
    ! and GLL-point locations in xstore,ystore,zstore
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...setting up jacobian '
      call flush_IMAIN()
    endif
    if (ANY_FAULT_IN_THIS_PROC) then
    ! compute jacobians with fault open and *store needed for ibool.
      call crm_ext_setup_jacobian(nodes_coords_open, nnodes_coords_open, elmnts_ext_mesh, nelmnts_ext_mesh)
    else ! with no fault
      call crm_ext_setup_jacobian(nodes_coords_ext_mesh, nnodes_ext_mesh, elmnts_ext_mesh, nelmnts_ext_mesh)
    endif

    ! user output
    ! call print_timing()

    ! creates ibool index array for projection from local to global points
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...indexing global points'
      call flush_IMAIN()
    endif
    if (ANY_FAULT_IN_THIS_PROC) then
      call crm_ext_setup_indexing(nnodes_coords_open, nodes_coords_open)
    else ! with no fault
      call crm_ext_setup_indexing(nnodes_ext_mesh, nodes_coords_ext_mesh)
    endif

    ! user output
    ! call print_timing()

    if (ANY_FAULT) then
      ! recalculate *store with faults closed
      call synchronize_all()
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) '  ...resetting up jacobian in fault domains'
        call flush_IMAIN()
      endif
      if (ANY_FAULT_IN_THIS_PROC) then
        call crm_ext_setup_jacobian(nodes_coords_ext_mesh, nnodes_ext_mesh, elmnts_ext_mesh, nelmnts_ext_mesh)
      endif
      ! at this point (xyz)store_unique are still open
      if (.not. PARALLEL_FAULT) then
        call fault_setup(ibool,nnodes_ext_mesh,nodes_coords_ext_mesh, &
                        xstore,ystore,zstore,nspec,nglob)
      endif
      ! this closes (xyz)store_unique

      ! user output
      ! call print_timing()
    endif


    ! sets up MPI interfaces between partitions
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...preparing MPI interfaces '
      call flush_IMAIN()
    endif
    call get_MPI_interface(nglob_unique,nspec,ibool)

    ! setting up parallel fault
    if (PARALLEL_FAULT .and. ANY_FAULT) then
      call synchronize_all()
      !at this point (xyz)store_unique are still open
      call fault_setup(ibool,nnodes_ext_mesh,nodes_coords_ext_mesh, &
                      xstore,ystore,zstore,nspec,nglob)
    ! this closes (xyz)store_unique
    endif

    ! user output
    ! call print_timing()

    ! sets up absorbing/free surface boundaries
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...setting up absorbing boundaries'
      call flush_IMAIN()
    endif
    call get_absorbing_boundary(nspec,ibool, &
                                nodes_coords_ext_mesh,nnodes_ext_mesh, &
                                ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                                nodes_ibelm_bottom,nodes_ibelm_top, &
                                nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                nspec2D_bottom,nspec2D_top)

    ! user output
    ! call print_timing()

    ! sets up mesh surface
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...setting up mesh surface'
      call flush_IMAIN()
    endif
    call crm_setup_mesh_surface()

    ! user output
    ! call print_timing()

    ! sets up up Moho surface
    if (SAVE_MOHO_MESH) then
      call synchronize_all()
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) '  ...setting up Moho surface'
        call flush_IMAIN()
      endif
      call crm_setup_moho(nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
                          nodes_coords_ext_mesh,nnodes_ext_mesh)

      ! user output
      ! call print_timing()
    endif

  end subroutine create_regions_mesh_coord

  subroutine crm_ext_deallocate_arrays()
    use create_regions_mesh_ext_par

    if (allocated(xstore_unique)) deallocate(xstore_unique)
    if (allocated(ystore_unique)) deallocate(ystore_unique)
    if (allocated(zstore_unique)) deallocate(zstore_unique)
    if (allocated(neighbors_adjncy)) deallocate(neighbors_adjncy)
    if (allocated(neighbors_xadj)) deallocate(neighbors_xadj)
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

  subroutine initialize_partition_arrays()
    use create_regions_mesh_ext_par

    if (allocated(mat_prop)) deallocate(mat_prop)
    if (allocated(nodes_coords_ext_mesh)) deallocate(nodes_coords_ext_mesh)
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

end module generate_databases_subs