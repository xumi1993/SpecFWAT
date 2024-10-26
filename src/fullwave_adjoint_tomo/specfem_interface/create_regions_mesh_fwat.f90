!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

  use generate_databases_par, only: nspec => NSPEC_AB,nglob => NGLOB_AB, &
      ibool,xstore,ystore,zstore, &
      npointot,myrank,LOCAL_PATH, &
      nnodes_ext_mesh,nelmnts_ext_mesh, &
      nodes_coords_ext_mesh, elmnts_ext_mesh, &
      max_memory_size,num_interfaces_ext_mesh, max_interface_size_ext_mesh, &
      my_neighbors_ext_mesh, my_nelmnts_neighbors_ext_mesh, &
      my_interfaces_ext_mesh, &
      ibool_interfaces_ext_mesh, nibool_interfaces_ext_mesh, &
      STACEY_ABSORBING_CONDITIONS, nspec2D_xmin, nspec2D_xmax, &
      nspec2D_ymin, nspec2D_ymax, &
      NSPEC2D_BOTTOM, NSPEC2D_TOP, &
      ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top, &
      nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
      nodes_ibelm_bottom,nodes_ibelm_top, &
      SAVE_MESH_FILES,PML_CONDITIONS, &
      ANISOTROPY,NPROC,APPROXIMATE_OCEAN_LOAD,OLSEN_ATTENUATION_RATIO, &
      ATTENUATION,USE_OLSEN_ATTENUATION, &
      nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
      ADIOS_FOR_MESH,IMAIN,SAVE_MOHO_MESH,ATTENUATION_f0_REFERENCE

  use create_regions_mesh_ext_par
  use fault_generate_databases, only: fault_read_input,fault_setup, &
                          fault_save_arrays,fault_save_arrays_test, &
                          nnodes_coords_open,nodes_coords_open,ANY_FAULT_IN_THIS_PROC, &
                          ANY_FAULT, PARALLEL_FAULT

  implicit none

! local parameters
! memory size needed by the solver
  double precision :: memory_size
  real(kind=CUSTOM_REAL) :: model_speed_max,min_resolved_period

! initializes arrays
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...allocating arrays '
    call flush_IMAIN()
  endif
  call crm_ext_allocate_arrays(nspec,LOCAL_PATH,myrank, &
                               nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                               nspec2D_bottom,nspec2D_top,ANISOTROPY, &
                               nodes_coords_ext_mesh,nnodes_ext_mesh,elmnts_ext_mesh,nelmnts_ext_mesh,ANY_FAULT_IN_THIS_PROC)

 ! if faults exist this reads nodes_coords_open
  call fault_read_input(prname,myrank)

! fills location and weights for Gauss-Lobatto-Legendre points, shape and derivations,
! returns jacobianstore,xixstore,...gammazstore
! and GLL-point locations in xstore,ystore,zstore
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...setting up jacobian '
    call flush_IMAIN()
  endif
  if (ANY_FAULT_IN_THIS_PROC) then
   ! compute jacobians with fault open and *store needed for ibool.
    call crm_ext_setup_jacobian(myrank,xstore,ystore,zstore,nspec, &
                                nodes_coords_open, nnodes_coords_open,elmnts_ext_mesh,nelmnts_ext_mesh)
  else ! with no fault
    call crm_ext_setup_jacobian(myrank,xstore,ystore,zstore,nspec, &
                                nodes_coords_ext_mesh,nnodes_ext_mesh,elmnts_ext_mesh,nelmnts_ext_mesh)
  endif


! creates ibool index array for projection from local to global points
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...indexing global points'
    call flush_IMAIN()
  endif
  if (ANY_FAULT_IN_THIS_PROC) then
    call crm_ext_setup_indexing(ibool, &
                                xstore,ystore,zstore,nspec,nglob,npointot, &
                                nnodes_coords_open,nodes_coords_open,myrank)
  else ! with no fault
    call crm_ext_setup_indexing(ibool, &
                                xstore,ystore,zstore,nspec,nglob,npointot, &
                                nnodes_ext_mesh,nodes_coords_ext_mesh,myrank)
  endif

  if (ANY_FAULT) then
    ! recalculate *store with faults closed
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ... resetting up jacobian in fault domains'
      call flush_IMAIN()
    endif
    if (ANY_FAULT_IN_THIS_PROC) then
      call crm_ext_setup_jacobian(myrank,xstore,ystore,zstore,nspec, &
                                  nodes_coords_ext_mesh,nnodes_ext_mesh, &
                                  elmnts_ext_mesh,nelmnts_ext_mesh)
    endif
    ! at this point (xyz)store_dummy are still open
    if (.not. PARALLEL_FAULT) then
      call fault_setup (ibool,nnodes_ext_mesh,nodes_coords_ext_mesh, &
                        xstore,ystore,zstore,nspec,nglob,myrank)
    endif
    ! this closes (xyz)store_dummy
  endif


! sets up MPI interfaces between partitions
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...preparing MPI interfaces '
    call flush_IMAIN()
  endif
  call get_MPI(myrank,nglob_dummy,nspec,ibool, &
               nelmnts_ext_mesh,elmnts_ext_mesh, &
               my_nelmnts_neighbors_ext_mesh, my_interfaces_ext_mesh, &
               ibool_interfaces_ext_mesh, &
               nibool_interfaces_ext_mesh, &
               num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
               my_neighbors_ext_mesh)

! setting up parallel fault
  if (PARALLEL_FAULT .and. ANY_FAULT) then
    call synchronize_all()
    !at this point (xyz)store_dummy are still open
    call fault_setup (ibool,nnodes_ext_mesh,nodes_coords_ext_mesh, &
                      xstore,ystore,zstore,nspec,nglob,myrank)
   ! this closes (xyz)store_dummy
  endif

! sets up absorbing/free surface boundaries
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...setting up absorbing boundaries'
    call flush_IMAIN()
  endif
  call get_absorbing_boundary(myrank,nspec,ibool, &
                              nodes_coords_ext_mesh,nnodes_ext_mesh, &
                              ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                              nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                              nodes_ibelm_bottom,nodes_ibelm_top, &
                              nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                              nspec2D_bottom,nspec2D_top)

! sets up mesh surface
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...setting up mesh surface'
    call flush_IMAIN()
  endif
  call crm_setup_mesh_surface()

! sets up up Moho surface
  if (SAVE_MOHO_MESH) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ...setting up Moho surface'
      call flush_IMAIN()
    endif
    call crm_setup_moho(myrank,nspec, &
                        nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,ibool )
  endif

! sets material velocities
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...determining velocity model'
    call flush_IMAIN()
  endif
  call get_model(myrank)

! sets up acoustic-elastic-poroelastic coupling surfaces
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...detecting acoustic-elastic-poroelastic surfaces '
    call flush_IMAIN()
  endif
  call get_coupling_surfaces(myrank, &
                             nspec,ibool,NPROC, &
                             nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                             num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                             my_neighbors_ext_mesh)

! locates inner and outer elements
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...element inner/outer separation '
    call flush_IMAIN()
  endif
  call crm_setup_inner_outer_elemnts(myrank,nspec, &
                                     num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                                     nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                     ibool,SAVE_MESH_FILES)

! colors mesh if requested
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...element mesh coloring '
    call flush_IMAIN()
  endif
  call setup_color_perm(myrank,nspec,nglob,ibool,ANISOTROPY,SAVE_MESH_FILES)

! overwrites material parameters from external binary files
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...external binary models '
    call flush_IMAIN()
  endif
  call get_model_binaries(myrank,nspec,LOCAL_PATH)

! calculates damping profiles and auxiliary coefficients on all C-PML points
  if (PML_CONDITIONS) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ...creating C-PML damping profiles '
      call flush_IMAIN()
    endif
    call pml_set_local_dampingcoeff(myrank,xstore_dummy,ystore_dummy,zstore_dummy)
  endif

! creates mass matrix
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...creating mass matrix '
    call flush_IMAIN()
  endif
  call create_mass_matrices(nglob_dummy,nspec,ibool,PML_CONDITIONS,STACEY_ABSORBING_CONDITIONS)

! saves the binary mesh files
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...saving databases'
    call flush_IMAIN()
  endif
  !call create_name_database(prname,myrank,LOCAL_PATH)
  if (ADIOS_FOR_MESH) then
    call save_arrays_solver_ext_mesh_adios(nspec,nglob,APPROXIMATE_OCEAN_LOAD, &
                                           ibool,num_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           max_interface_size_ext_mesh, &
                                           ibool_interfaces_ext_mesh,SAVE_MESH_FILES,ANISOTROPY)
  else
    call save_arrays_solver_ext_mesh(nspec,nglob_dummy,APPROXIMATE_OCEAN_LOAD,ibool, &
                                     num_interfaces_ext_mesh,my_neighbors_ext_mesh,nibool_interfaces_ext_mesh, &
                                     max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                                     SAVE_MESH_FILES,ANISOTROPY)
  endif

! saves faults
  if (ANY_FAULT) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ...saving fault databases'
      call flush_IMAIN()
    endif
    !  call fault_save_arrays_test(prname)  ! for debugging
    call fault_save_arrays(prname)
  endif

! saves moho surface
  if (SAVE_MOHO_MESH) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ...saving Moho surfaces'
      call flush_IMAIN()
    endif
    call crm_save_moho()
  endif

! computes the approximate amount of memory needed to run the solver
  call synchronize_all()
  call memory_eval(nspec,nglob_dummy,maxval(nibool_interfaces_ext_mesh),num_interfaces_ext_mesh, &
                   APPROXIMATE_OCEAN_LOAD,memory_size)
  call max_all_dp(memory_size, max_memory_size)

! checks the mesh, stability and resolved period
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...checking mesh resolution'
    call flush_IMAIN()
  endif

  if (POROELASTIC_SIMULATION) then
    !chris: check for poro: At the moment cpI & cpII are for eta=0
    call check_mesh_resolution_poro(myrank,nspec,nglob_dummy,ibool, &
                                    xstore_dummy,ystore_dummy,zstore_dummy, &
                                    -1.0d0, model_speed_max,min_resolved_period, &
                                    phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                                    LOCAL_PATH,SAVE_MESH_FILES )
  else
    call check_mesh_resolution(myrank,nspec,nglob_dummy, &
                               ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
                               kappastore,mustore,rho_vp,rho_vs, &
                               -1.0d0,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)
  endif

! saves binary mesh files for attenuation
  if (ATTENUATION) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ...saving attenuation databases'
      call flush_IMAIN()
    endif
    call get_attenuation_model(myrank,nspec,USE_OLSEN_ATTENUATION,OLSEN_ATTENUATION_RATIO, &
                               mustore,rho_vs,kappastore,rho_vp,qkappa_attenuation_store,qmu_attenuation_store, &
                               ispec_is_elastic,min_resolved_period,prname,ATTENUATION_f0_REFERENCE)
  endif

  ! cleanup
  deallocate(xixstore,xiystore,xizstore, &
             etaxstore,etaystore,etazstore, &
             gammaxstore,gammaystore,gammazstore)
  deallocate(jacobianstore)
  deallocate(qkappa_attenuation_store,qmu_attenuation_store)
  deallocate(kappastore,mustore,rho_vp,rho_vs)
  deallocate(rho_vpI,rho_vpII,rho_vsI)
  deallocate(rhoarraystore,kappaarraystore,etastore,phistore,tortstore,permstore)

  if (.not. SAVE_MOHO_MESH) then
    deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
  endif

  if (ACOUSTIC_SIMULATION) then
    deallocate(rmass_acoustic)
  endif

  if (ELASTIC_SIMULATION) then
    deallocate(rmass)
  endif

  if (POROELASTIC_SIMULATION) then
    deallocate(rmass_solid_poroelastic,rmass_fluid_poroelastic)
  endif

  if (STACEY_ABSORBING_CONDITIONS) then
     if (ELASTIC_SIMULATION) then
       deallocate(rmassx,rmassy,rmassz)
     endif
     if (ACOUSTIC_SIMULATION) then
       deallocate(rmassz_acoustic)
     endif
  endif

end subroutine create_regions_mesh_fwat



  subroutine get_model_binaries_fwat(myrank,nspec,LOCAL_PATH)

! reads in material parameters from external binary files

  use generate_databases_par, only: IMAIN, IMODEL, IMODEL_GLL,IMODEL_IPATI,IMODEL_IPATI_WATER,&
                                    IMODEL_SEP, ADIOS_FOR_MESH, IMODEL_USER_EXTERNAL

  use create_regions_mesh_ext_par

  use model_ipati_adios_mod, only: model_ipati_adios
  use fullwave_adjoint_tomo_par, only: GRID_PATH

  implicit none

  ! number of spectral elements in each block
  integer :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! external GLL models
  ! variables for importing models from files in SPECFEM format, e.g.,  proc000000_vp.bin etc.
  ! can be used for importing updated model in iterative inversions

  ! note: we read in these binary files after mesh coloring, since mesh coloring is permuting arrays.
  !          here, the ordering in **_vp.bin etc. can be permuted as they are outputted when saving mesh files

  select case (IMODEL)
  case (IMODEL_USER_EXTERNAL)
    call model_grid(myrank,nspec,LOCAL_PATH, GRID_PATH)
  case (IMODEL_GLL)
    ! note:
    ! import the model from files in SPECFEM format
    ! note that those those files should be saved in LOCAL_PATH
    if (ADIOS_FOR_MESH) then
      call model_gll_adios(myrank,nspec,LOCAL_PATH)
    else
      call model_gll(myrank,nspec,LOCAL_PATH)
    endif

  case (IMODEL_IPATI)
    ! import the model from modified files in SPECFEM format
    if (ADIOS_FOR_MESH) then
      call model_ipati_adios(myrank,nspec,LOCAL_PATH)
    else
      call model_ipati(myrank,nspec,LOCAL_PATH)
    endif

  case (IMODEL_IPATI_WATER)
    ! import the model from modified files in SPECFEM format
    if (ADIOS_FOR_MESH) then
      call model_ipati_water_adios(myrank,nspec,LOCAL_PATH)
    else
      call model_ipati_water(myrank,nspec,LOCAL_PATH)
    endif

  case default
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '     no external binary model used '
    endif
  end select

  end subroutine get_model_binaries_fwat


  subroutine model_grid(myrank,nspec,LOCAL_PATH,filename)

  use constants
  use utils
  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IMAIN,MAX_STRING_LEN,ATTENUATION

  use create_regions_mesh_ext_par, only: rhostore,kappastore,mustore,rho_vp,rho_vs,qkappa_attenuation_store,qmu_attenuation_store
  use projection_on_FD_grid_fwat

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: vp_read,vs_read,rho_read, qmu_read, qkappa_read
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: vp_gll, vs_gll, rho_gll
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: x, y, z
  integer :: ier
  character(len=*), intent(in) :: filename

  ! read h5file
  if (myrank == 0) then
    write(IMAIN,*) '     reading in: ',trim(filename)
    call h5read(filename, '/vp', vp_read)
    vp_read = transpose_3(vp_read)
    call h5read(filename, '/vs', vs_read)
    vs_read = transpose_3(vs_read)
    call h5read(filename, '/rho', rho_read)
    rho_read = transpose_3(rho_read)
    if (ATTENUATION) then
      call h5read(filename, '/qmu', qmu_read)
      qmu_read = transpose_3(qmu_read)
      call h5read(filename, '/qkappa', qkappa_read)
      qkappa_read = transpose_3(qkappa_read)
    endif
    call h5read(filename, '/x', x)
    call h5read(filename, '/y', y)
    call h5read(filename, '/z', z)
    ox_fd_proj = x(1); hx_fd_proj = x(2) - x(1); nx_fd_proj = size(x)
    oy_fd_proj = y(1); hy_fd_proj = y(2) - y(1); ny_fd_proj = size(y)
    oz_fd_proj = z(1); hz_fd_proj = z(2) - z(1); nz_fd_proj = size(z)
  endif
  call synchronize_all()
  call bcast_all_singlei(nx_fd_proj)
  call bcast_all_singlei(ny_fd_proj)
  call bcast_all_singlei(nz_fd_proj)
  call bcast_all_singlecr(ox_fd_proj)
  call bcast_all_singlecr(hx_fd_proj)
  call bcast_all_singlecr(oy_fd_proj)
  call bcast_all_singlecr(hy_fd_proj)
  call bcast_all_singlecr(oz_fd_proj)
  call bcast_all_singlecr(hz_fd_proj)
  if (myrank /= 0) then
    allocate(vp_read(nx_fd_proj, ny_fd_proj, nz_fd_proj))
    allocate(vs_read(nx_fd_proj, ny_fd_proj, nz_fd_proj))
    allocate(rho_read(nx_fd_proj, ny_fd_proj, nz_fd_proj))
  endif
  call bcast_all_cr(vp_read, nx_fd_proj*ny_fd_proj*nz_fd_proj)
  call bcast_all_cr(vs_read, nx_fd_proj*ny_fd_proj*nz_fd_proj)
  call bcast_all_cr(rho_read, nx_fd_proj*ny_fd_proj*nz_fd_proj)
  if (ATTENUATION) then
    if (myrank /= 0) then
      allocate(qmu_read(nx_fd_proj, ny_fd_proj, nz_fd_proj))
      allocate(qkappa_read(nx_fd_proj, ny_fd_proj, nz_fd_proj))
    endif
    call bcast_all_cr(qmu_read, nx_fd_proj*ny_fd_proj*nz_fd_proj)
    call bcast_all_cr(qkappa_read, nx_fd_proj*ny_fd_proj*nz_fd_proj)
  endif

  allocate(vp_gll(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(vs_gll(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho_gll(NGLLX,NGLLY,NGLLZ,nspec))

  call Project_model_FD_grid2SEM(vp_gll, vp_read, myrank)
  call Project_model_FD_grid2SEM(vs_gll, vs_read, myrank)
  call Project_model_FD_grid2SEM(rho_gll, rho_read, myrank)

  rhostore(:,:,:,:) = rho_gll(:,:,:,:)
  kappastore(:,:,:,:) = rhostore(:,:,:,:) * ( vp_gll(:,:,:,:) * vp_gll(:,:,:,:) &
                                              - FOUR_THIRDS * vs_gll(:,:,:,:) * vs_gll(:,:,:,:) )
  mustore(:,:,:,:) = rhostore(:,:,:,:) * vs_gll(:,:,:,:) * vs_gll(:,:,:,:)

  rho_vp(:,:,:,:) = rhostore(:,:,:,:) * vp_gll(:,:,:,:)
  rho_vs(:,:,:,:) = rhostore(:,:,:,:) * vs_gll(:,:,:,:)

  if (ATTENUATION) then
    call Project_model_FD_grid2SEM(qmu_attenuation_store, qmu_read, myrank)
    call Project_model_FD_grid2SEM(qkappa_attenuation_store, qkappa_read, myrank)
  endif

  end subroutine
