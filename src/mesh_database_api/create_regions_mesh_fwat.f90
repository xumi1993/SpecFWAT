!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine create_regions_mesh_fwat()

! create the different regions of the mesh

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
  call print_timing()

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
  call print_timing()

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
  call print_timing()

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
    call print_timing()
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
  call print_timing()

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
  call print_timing()

  ! sets up mesh surface
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...setting up mesh surface'
    call flush_IMAIN()
  endif
  call crm_setup_mesh_surface()

  ! user output
  call print_timing()

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
    call print_timing()
  endif

  ! sets material velocities
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...determining velocity model'
    call flush_IMAIN()
  endif
  call get_model_fwat()

  ! user output
  call print_timing()

  ! sets up acoustic-elastic-poroelastic coupling surfaces
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...detecting acoustic-elastic-poroelastic surfaces '
    call flush_IMAIN()
  endif
  call get_coupling_surfaces(nspec,ibool)

  ! user output
  call print_timing()

  ! locates inner and outer elements
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...element inner/outer separation '
    call flush_IMAIN()
  endif
  call crm_setup_inner_outer_elemnts()

  ! user output
  call print_timing()

  ! LTS mode
  if (LTS_MODE) then
    call synchronize_all()
    if ( myrank == 0) then
      write(IMAIN,*) '  ...determining LTS arrays'
      call flush_IMAIN()
    endif
    ! sets up elements based on current element size and velocity model
    call lts_generate_databases()

    ! user output
    call print_timing()
  endif

  ! colors mesh if requested
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...element mesh coloring '
    call flush_IMAIN()
  endif
  call setup_color_perm(nspec,nglob,ibool,ANISOTROPY,SAVE_MESH_FILES)

  ! user output
  call print_timing()

  ! overwrites material parameters from external binary files
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...external binary models '
    call flush_IMAIN()
  endif
  call get_model_binaries(nspec,LOCAL_PATH)

  ! user output
  call print_timing()

  ! calculates damping profiles and auxiliary coefficients on all C-PML points
  if (PML_CONDITIONS) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...creating C-PML damping profiles '
      call flush_IMAIN()
    endif
    call pml_set_local_dampingcoeff(xstore_unique,ystore_unique,zstore_unique)

    ! user output
    call print_timing()
  endif

  ! creates mass matrix
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...creating mass matrix '
    call flush_IMAIN()
  endif
  call create_mass_matrices(nglob_unique)

  ! user output
  call print_timing()

  ! sets up mesh adjacency for point searches
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...setting up mesh adjacency '
    call flush_IMAIN()
  endif
  call setup_mesh_adjacency()

  ! user output
  call print_timing()

  ! setup wavefield discontinuity interface
  if (IS_WAVEFIELD_DISCONTINUITY) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...setting up wavefield discontinuity boundary '
      call flush_IMAIN()
    endif
    call read_partition_files_wavefield_discontinuity()
    call setup_boundary_wavefield_discontinuity()
  endif

  ! saves the binary mesh files
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...saving mesh databases'
    call flush_IMAIN()
  endif
  call save_arrays_solver_mesh()

  ! user output
  call print_timing()

  ! saves faults
  if (ANY_FAULT) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...saving fault databases'
      call flush_IMAIN()
    endif
    ! for debugging
    !call fault_save_arrays_test(prname)

    ! saves fault arrays
    call fault_save_arrays(prname)

    ! user output
    call print_timing()
  endif

  ! saves moho surface
  if (SAVE_MOHO_MESH) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...saving Moho surfaces'
      call flush_IMAIN()
    endif
    call crm_save_moho()

    ! user output
    call print_timing()
  endif
  call synchronize_all()

  ! computes the approximate amount of memory needed to run the solver
  call memory_eval(nspec,nglob_unique,maxval(nibool_interfaces_ext_mesh),num_interfaces_ext_mesh, &
                   APPROXIMATE_OCEAN_LOAD,memory_size)

  call max_all_dp(memory_size, max_memory_size)

  ! checks the mesh, stability and resolved period
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...checking mesh resolution'
    call flush_IMAIN()
  endif
  call check_mesh_resolution(nspec,nglob_unique, &
                             ibool,xstore_unique,ystore_unique,zstore_unique, &
                             ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
                             kappastore,mustore,rhostore, &
                             phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                             -1.0d0,model_speed_max,min_resolved_period)

  ! user output
  call print_timing()

  ! saves binary mesh files for attenuation
  if (ATTENUATION) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...saving attenuation databases'
      call flush_IMAIN()
    endif
    call get_attenuation_model(nspec,USE_OLSEN_ATTENUATION,OLSEN_ATTENUATION_RATIO, &
                               mustore,rho_vs,kappastore,rho_vp, &
                               qkappa_attenuation_store,qmu_attenuation_store, &
                               ispec_is_elastic,min_resolved_period,prname,ATTENUATION_f0_REFERENCE)

    ! user output
    call print_timing()
  endif

  ! synchronizes processes, making sure everybody has finished
  call synchronize_all()

  ! cleanup
  call pml_cleanup_arrays()

  ! model arrays
  deallocate(xixstore,xiystore,xizstore, &
             etaxstore,etaystore,etazstore, &
             gammaxstore,gammaystore,gammazstore,jacobianstore)
  deallocate(qkappa_attenuation_store,qmu_attenuation_store)
  deallocate(kappastore,mustore,rhostore,rho_vp,rho_vs)
  deallocate(rho_vpI,rho_vpII,rho_vsI)
  deallocate(rhoarraystore,kappaarraystore,etastore,phistore,tortstore,permstore)

  if (.not. SAVE_MOHO_MESH) deallocate(xstore_unique,ystore_unique,zstore_unique)

  if (ACOUSTIC_SIMULATION) deallocate(rmass_acoustic)
  if (ELASTIC_SIMULATION) deallocate(rmass)
  if (POROELASTIC_SIMULATION) deallocate(rmass_solid_poroelastic,rmass_fluid_poroelastic)

  if (STACEY_ABSORBING_CONDITIONS) then
    if (ELASTIC_SIMULATION) deallocate(rmassx,rmassy,rmassz)
    if (ACOUSTIC_SIMULATION) deallocate(rmassz_acoustic)
  endif

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'mesh regions done'
    call flush_IMAIN()
  endif

contains

  ! helper routine to output timing info
  subroutine print_timing()
    implicit none
    double precision :: tCPU
    ! user output
    if (DEBUG_TIMING) then
      if (myrank == 0) then
        tCPU = wtime() - time_start
        write(IMAIN,*) "     time in seconds = ", tCPU
        call flush_IMAIN()
      endif
      time_start = wtime()
    endif
  end subroutine print_timing

  end subroutine create_regions_mesh_fwat
