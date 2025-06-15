module meshfem3D_subs
  use meshfem_par
  use chunk_earth_mod
  use create_meshfem_par
  use config, only : local_path_fwat

  implicit none

  include 'version.fh'

contains

  subroutine meshfem3D_fwat(filename)
    ! local parameters
    integer :: iprocnum
    integer :: iproc_xi,iproc_eta
    integer :: ier
    integer(kind=8) :: nglob_total
    character(len=3) :: str_unit
    character(len=*), intent(in) :: filename

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
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'creating chunk of the Earth mesh'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

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

    call read_mesh_parameter_file_fwat(filename)
    LOCAL_PATH = local_path_fwat
    
    ! get interface data from external file to count the spectral elements along Z
    if (allocated(ner_layer)) deallocate(ner_layer)
    call get_interfaces_mesh_count()

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
                            NGLOB_AB, &
                            USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

    ! check that the code is running with the requested nb of processes
    if (sizeprocs /= NPROC) then
      if (myrank == 0) then
        write(IMAIN,*) 'Error: number of processors supposed to run on: ',NPROC
        write(IMAIN,*) 'Error: number of MPI processors actually run on: ',sizeprocs
        print *
        print *, 'Error meshfem3D: number of processors supposed to run on: ',NPROC
        print *, 'Error meshfem3D: number of MPI processors actually run on: ',sizeprocs
        print *
      endif
      call exit_MPI(myrank,'wrong number of MPI processes')
    endif
    call synchronize_all()

    call initialize_mesh_arrays()
    ! dynamic allocation of mesh arrays
    allocate(xgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1353')
    if (ier /= 0) stop 'Error allocating array xgrid'
    xgrid(:,:,:) = 0.d0

    allocate(ygrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1354')
    if (ier /= 0) stop 'Error allocating array ygrid'
    ygrid(:,:,:) = 0.d0

    allocate(zgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1355')
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
    zgrid(:,:,:) = 0.d0

    allocate(addressing(0:NPROC_XI-1,0:NPROC_ETA-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1356')
    if (ier /= 0) stop 'Error allocating array addressing'
    addressing(:,:) = 0

    allocate(iproc_xi_slice(0:NPROC-1), &
            iproc_eta_slice(0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1358')
    if (ier /= 0) stop 'Error allocating array iproc_eta_slice'
    ! clear arrays
    iproc_xi_slice(:) = 0
    iproc_eta_slice(:) = 0

    ! create global slice addressing for solver
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Creating global slice addressing'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    do iproc_eta = 0,NPROC_ETA-1
      do iproc_xi = 0,NPROC_XI-1
        iprocnum = iproc_eta * NPROC_XI + iproc_xi
        iproc_xi_slice(iprocnum) = iproc_xi
        iproc_eta_slice(iprocnum) = iproc_eta
        addressing(iproc_xi,iproc_eta) = iprocnum
      enddo
    enddo

    if (myrank == 0) then
      write(IMAIN,*) 'Spatial distribution of slice numbers:'
      do iproc_eta = NPROC_ETA-1, 0, -1
        do iproc_xi = 0, NPROC_XI-1, 1
          write(IMAIN,'(i5)',advance='no') addressing(iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(a1)',advance='yes') ' '
      enddo
      call flush_IMAIN()
    endif

    if (myrank == 0) then
      write(IMAIN,*) 'This is process ',myrank
      write(IMAIN,*) 'There are ',sizeprocs,' MPI processes'
      write(IMAIN,*) 'Processes are numbered from 0 to ',sizeprocs-1
      write(IMAIN,*)
      write(IMAIN,*) 'There are ',NEX_XI,' elements along xi'
      write(IMAIN,*) 'There are ',NEX_ETA,' elements along eta'
      write(IMAIN,*) 'There are ',NER,' elements along Z'
      write(IMAIN,*)
      do ilayer = 1,number_of_layers
        write(IMAIN,*) 'There are ',ner_layer(ilayer),' spectral elements along Z in layer ',ilayer
      enddo
      write(IMAIN,*)
      write(IMAIN,*) 'There are ',NPROC_XI,' slices along xi'
      write(IMAIN,*) 'There are ',NPROC_ETA,' slices along eta'
      write(IMAIN,*) 'There is a total of ',NPROC,' slices'

      write(IMAIN,*)
      write(IMAIN,*) 'Shape functions defined by NGNOD = ',NGNOD,' control nodes'
      write(IMAIN,*) 'Surface shape functions defined by NGNOD2D = ',NGNOD2D,' control nodes'
      write(IMAIN,*) 'Beware! Curvature (i.e. HEX27 elements) is not handled by our internal mesher'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! check that the constants.h file is correct
    if (NGNOD /= 8 .and. NGNOD /= 27) &
      call exit_MPI(myrank,'Error must have set NGNOD == 8 or NGNOD == 27')
    if (NGNOD2D /= 4 .and. NGNOD2D /= 9) &
      call exit_MPI(myrank,'Error must have set NGNOD2D == 4 or NGNOD2D == 9')

    if (NGLLX_M == 2 .and. NGLLY_M == 2 .and. NGLLZ_M == 2) then
      if (NGNOD /= 8) &
        call exit_MPI(myrank,'With NGLLX_M == 2, volume elements should have NGNOD == 8 control nodes in our internal mesher')
      if (NGNOD2D /= 4) &
        call exit_MPI(myrank,'With NGLLX_M == 2, surface elements should have NGNOD2D == 4 control nodes in our internal mesher')
    endif

    if (NGNOD == 27 .and. (NGLLX_M < 3 .or. NGLLY_M < 3 .or. NGLLZ_M < 3)) &
      call exit_MPI(myrank,'NGNOD = 27 control nodes needs at least NGLLX_M == NGLLY_M == NGLLZ_M >= 3 in our internal mesher')
    if (NGNOD2D == 9 .and. (NGLLX_M < 3 .or. NGLLY_M < 3 .or. NGLLZ_M < 3)) &
      call exit_MPI(myrank,'NGNOD2D = 9 control nodes needs at least NGLLX_M == NGLLY_M == NGLLZ_M >= 3 in our internal mesher')

  !  if (.not. USE_REGULAR_MESH .and. NGLLX_M >= 3) &
  !    call exit_MPI(myrank,'NGLLX_M == NGLLY_M == NGLLZ_M >= 3 only supported with USE_REGULAR_MESH = .true. at the moment')

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

    if (myrank == 0) then
      write(IMAIN,*) 'region selected:'
      write(IMAIN,*)
      write(IMAIN,*) 'latitude min = ',LATITUDE_MIN
      write(IMAIN,*) 'latitude max = ',LATITUDE_MAX
      write(IMAIN,*)
      write(IMAIN,*) 'longitude min = ',LONGITUDE_MIN
      write(IMAIN,*) 'longitude max = ',LONGITUDE_MAX
      write(IMAIN,*)
      if (SUPPRESS_UTM_PROJECTION) then
        write(IMAIN,*) 'this is given directly as UTM'
      else
        write(IMAIN,*) 'this is mapped to UTM in region ',UTM_PROJECTION_ZONE
      endif
      write(IMAIN,*)
      write(IMAIN,*) 'UTM X min = ',UTM_X_MIN
      write(IMAIN,*) 'UTM X max = ',UTM_X_MAX
      write(IMAIN,*)
      write(IMAIN,*) 'UTM Y min = ',UTM_Y_MIN
      write(IMAIN,*) 'UTM Y max = ',UTM_Y_MAX
      write(IMAIN,*)
      write(IMAIN,*) 'UTM size of model along X is ',(UTM_X_MAX-UTM_X_MIN)/1000.,' km'
      write(IMAIN,*) 'UTM size of model along Y is ',(UTM_Y_MAX-UTM_Y_MIN)/1000.,' km'
      write(IMAIN,*)
      write(IMAIN,*) 'Bottom of the mesh is at a depth of ',dabs(Z_DEPTH_BLOCK)/1000.,' km'
      write(IMAIN,*)
      write(IMAIN,*)
      if (SUPPRESS_UTM_PROJECTION) then
        write(IMAIN,*) 'suppressing UTM projection'
      else
        write(IMAIN,*) 'using UTM projection in region ',UTM_PROJECTION_ZONE
      endif
      if (PML_CONDITIONS) then
        if (SUPPRESS_UTM_PROJECTION) then
          ! no UTM, thickness given in m
          str_unit = '(m)'
        else
          ! UTM, thickness given in degree
          str_unit = 'deg'
        endif
        write(IMAIN,*)
        write(IMAIN,*) 'PML thickness in X direction = ',sngl(THICKNESS_OF_X_PML),str_unit
        write(IMAIN,*) 'PML thickness in Y direction = ',sngl(THICKNESS_OF_Y_PML),str_unit
        write(IMAIN,*) 'PML thickness in Z direction = ',sngl(THICKNESS_OF_Z_PML),str_unit
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! get addressing for this process
    iproc_xi_current = iproc_xi_slice(myrank)
    iproc_eta_current = iproc_eta_slice(myrank)

    ! number of elements in each slice
    npx_element_steps = 2*NEX_PER_PROC_XI
    npy_element_steps = 2*NEX_PER_PROC_ETA
    ner_layer(:) = 2 * ner_layer(:)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '**************************'
      write(IMAIN,*) 'Creating interfaces'
      write(IMAIN,*) '**************************'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! creates mesh interfaces
    call create_interfaces_mesh()

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '**************************'
      write(IMAIN,*) 'Creating mesh in the model'
      write(IMAIN,*) '**************************'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! creates mesh element points
    call create_meshfem_mesh()

    ! make sure everybody is synchronized
    call synchronize_all()

    ! to avoid overflow for large meshes
    nglob_total = int( dble(NGLOB_AB)*dble(NPROC) ,kind=8)

    !--- print number of points and elements in the mesh
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Repartition of elements:'
      write(IMAIN,*) '-----------------------'
      write(IMAIN,*)
      write(IMAIN,*) 'total number of elements in mesh slice 0: ',NSPEC_AB
      write(IMAIN,*) 'total number of points in mesh slice 0: ',NGLOB_AB
      write(IMAIN,*)
      write(IMAIN,*) 'total number of elements in entire mesh: ',NSPEC_AB*NPROC
      write(IMAIN,*) 'approximate total number of points in entire mesh (with duplicates on MPI edges): ',nglob_total
      write(IMAIN,*) 'approximate total number of DOFs in entire mesh (with duplicates on MPI edges): ',nglob_total*NDIM
      write(IMAIN,*)
      ! write information about precision used for floating-point operations
      if (CUSTOM_REAL == SIZE_REAL) then
        write(IMAIN,*) 'using single precision for the calculations'
      else
        write(IMAIN,*) 'using double precision for the calculations'
      endif
      write(IMAIN,*)
      write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
      write(IMAIN,*)
      call flush_IMAIN()
    endif   ! end of section executed by main process only

    if (myrank == 0) close(IMAIN)
  end subroutine meshfem3D_fwat

  subroutine read_mesh_parameter_file_fwat(filename)
    use meshfem_par, only: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
      UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
      NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE, &
      LOCAL_PATH,SUPPRESS_UTM_PROJECTION, &
      INTERFACES_FILE,CAVITY_FILE,NSUBREGIONS,subregions, &
      NMATERIALS,material_properties,material_properties_undef, &
      CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES, &
      USE_REGULAR_MESH,NDOUBLINGS,ner_doublings, &
      THICKNESS_OF_X_PML,THICKNESS_OF_Y_PML,THICKNESS_OF_Z_PML, &
      myrank,sizeprocs,NUMBER_OF_MATERIAL_PROPERTIES, &
      SAVE_MESH_AS_CUBIT

    use constants, only: IIN,MF_IN_DATA_FILES,MAX_STRING_LEN,IMAIN, &
      IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC, &
      ILONGLAT2UTM,IGNORE_JUNK,DONT_IGNORE_JUNK

    implicit none

    integer :: NEX_MAX
    double precision :: UTM_MAX
    double precision :: DEPTH_BLOCK_KM

    integer :: ier

    integer :: ix_beg_region,ix_end_region,iy_beg_region,iy_end_region
    integer :: iz_beg_region,iz_end_region,imaterial_number

    double precision :: rho,vp,vs,Q_Kappa,Q_mu,anisotropy_flag

    integer :: ireg,imat,ndef,nundef
    integer :: mat_id,domain_id
    integer :: ner_value,ner_max,idoub
    logical :: found
    character(len=*), intent(in) :: filename


    if (myrank == 0) then
      write(IMAIN,*) 'Reading mesh parameters from file ',trim(filename)
      call flush_IMAIN()
    endif

    ! note: please be careful that the order of reading in parameters here
    !       must match the order of appearance in the Mesh_Par_file
    !
    ! open parameter file Mesh_Par_file
    call open_parameter_file_mesh(filename)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  input parameters...'
      call flush_IMAIN()
    endif

    call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LATITUDE_MIN, 'LATITUDE_MIN', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter LATITUDE_MIN'
    call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LATITUDE_MAX, 'LATITUDE_MAX', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter LATITUDE_MAX'
    call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LONGITUDE_MIN, 'LONGITUDE_MIN', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter LONGITUDE_MIN'
    call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LONGITUDE_MAX, 'LONGITUDE_MAX', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter LONGITUDE_MAX'
    call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,DEPTH_BLOCK_KM, 'DEPTH_BLOCK_KM', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter DEPTH_BLOCK_KM'
    call read_value_integer_mesh(IIN,IGNORE_JUNK,UTM_PROJECTION_ZONE, 'UTM_PROJECTION_ZONE', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter UTM_PROJECTION_ZONE'
    call read_value_logical_mesh(IIN,IGNORE_JUNK,SUPPRESS_UTM_PROJECTION, 'SUPPRESS_UTM_PROJECTION', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter SUPPRESS_UTM_PROJECTION'
    call read_value_string_mesh(IIN,IGNORE_JUNK,INTERFACES_FILE, 'INTERFACES_FILE', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter INTERFACES_FILE'
    call read_value_string_mesh(IIN,IGNORE_JUNK,CAVITY_FILE, 'CAVITY_FILE', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter CAVITY_FILE'

    call read_value_integer_mesh(IIN,IGNORE_JUNK,NEX_XI, 'NEX_XI', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter NEX_XI'
    call read_value_integer_mesh(IIN,IGNORE_JUNK,NEX_ETA, 'NEX_ETA', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter NEX_ETA'
    call read_value_integer_mesh(IIN,IGNORE_JUNK,NPROC_XI, 'NPROC_XI', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter NPROC_XI'
    call read_value_integer_mesh(IIN,IGNORE_JUNK,NPROC_ETA, 'NPROC_ETA', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter NPROC_ETA'

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  doubling layers...'
      call flush_IMAIN()
    endif

    call read_value_logical_mesh(IIN,IGNORE_JUNK,USE_REGULAR_MESH, 'USE_REGULAR_MESH', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter USE_REGULAR_MESH'
    call read_value_integer_mesh(IIN,IGNORE_JUNK,NDOUBLINGS, 'NDOUBLINGS', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter NDOUBLINGS'

    ! checks value
    if (USE_REGULAR_MESH) then
      ! sets NDOUBLINGS to zero for regular grids
      NDOUBLINGS = 0
    else
      ! irregular grid with doubling layer
      if (NDOUBLINGS < 1) stop 'Error parameter NDOUBLINGS is invalid! value must be at least 1 for irregular grids'
    endif

    ! allocate doubling array
    if (allocated(ner_doublings)) deallocate(ner_doublings)
    allocate(ner_doublings(NDOUBLINGS),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1318')
    if (ier /= 0) stop 'Error allocating ner_doublings array'
    ner_doublings(:) = 0

    do idoub = 1,NDOUBLINGS
      call read_value_doubling_integer_mesh(IIN,DONT_IGNORE_JUNK, ner_value, 'NZ_DOUBLING', ier)
      if (ier /= 0) then
        print *,'Error reading doubling entry for doubling layer: ',idoub
        print *,'Please check NDOUBLINGS value and corresponding lines with NZ_DOUBLING entries'
        stop 'Error reading Mesh parameter NZ_DOUBLING'
      endif
      ner_doublings(idoub) = ner_value
    enddo
    ! jump over unused doubling entries to reach lines with visualization parameters below
    call read_value_doubling_skip_mesh(IIN,DONT_IGNORE_JUNK,'NZ_DOUBLING', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter after NDOUBLINGS'

    ! only fix 2 doubling layer entries
    !call read_value_integer_mesh(IIN,IGNORE_JUNK, ner_doublings(1), 'NZ_DOUBLING_1', ier)
    !if (ier /= 0) stop 'Error reading Mesh parameter NZ_DOUBLING_1'
    !call read_value_integer_mesh(IIN,IGNORE_JUNK, ner_doublings(2), 'NZ_DOUBLING_2', ier)
    !if (ier /= 0) stop 'Error reading Mesh parameter NZ_DOUBLING_2'

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  visualization...'
      call flush_IMAIN()
    endif

    ! visualization file output
    call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_ABAQUS_FILES, 'CREATE_ABAQUS_FILES', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter CREATE_ABAQUS_FILES'
    call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_DX_FILES, 'CREATE_DX_FILES', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter CREATE_DX_FILES'
    call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_VTK_FILES, 'CREATE_VTK_FILES', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter CREATE_VTK_FILES'

    ! save as CUBIT-mesh
    call read_value_logical_mesh(IIN,IGNORE_JUNK,SAVE_MESH_AS_CUBIT, 'SAVE_MESH_AS_CUBIT', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter SAVE_MESH_AS_CUBIT'

    ! file in which we store the databases
    call read_value_string_mesh(IIN,IGNORE_JUNK,LOCAL_PATH, 'LOCAL_PATH', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter LOCAL_PATH'

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  CPML...'
      call flush_IMAIN()
    endif

    ! CPML thickness
    call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,THICKNESS_OF_X_PML, 'THICKNESS_OF_X_PML', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter THICKNESS_OF_X_PML'
    call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,THICKNESS_OF_Y_PML, 'THICKNESS_OF_Y_PML', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter THICKNESS_OF_Y_PML'
    call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,THICKNESS_OF_Z_PML, 'THICKNESS_OF_Z_PML', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter THICKNESS_OF_Z_PML'

    if (THICKNESS_OF_X_PML < 0.0) stop 'Error invalid negative value THICKNESS_OF_X_PML'
    if (THICKNESS_OF_Y_PML < 0.0) stop 'Error invalid negative value THICKNESS_OF_Y_PML'
    if (THICKNESS_OF_Z_PML < 0.0) stop 'Error invalid negative value THICKNESS_OF_Z_PML'

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  domain materials...'
      call flush_IMAIN()
    endif

    ! read number of materials
    call read_value_integer_mesh(IIN,IGNORE_JUNK,NMATERIALS, 'NMATERIALS', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter NMATERIALS'

    ! read materials properties
    if (allocated(material_properties)) deallocate(material_properties)
    allocate(material_properties(NMATERIALS,NUMBER_OF_MATERIAL_PROPERTIES),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1319')
    if (ier /= 0) stop 'Error allocation of material_properties'
    material_properties(:,:) = 0.d0

    if (allocated(material_properties_undef)) deallocate(material_properties_undef)
    allocate(material_properties_undef(NMATERIALS,3),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1319b')
    if (ier /= 0) stop 'Error allocation of material_properties_undef'
    material_properties_undef(:,:) = ""

    do imat = 1,NMATERIALS
      call read_material_parameters(IIN,material_properties,material_properties_undef,imat,NMATERIALS,ier)
      if (ier /= 0) then
        print *,'Error reading material ',imat,' out of ',NMATERIALS
        stop 'Error reading materials in Mesh_Par_file'
      endif
      ! user output
      domain_id = material_properties(imat,7)
      mat_id = material_properties(imat,8)
      if (myrank == 0) then
        select case(domain_id)
        case(IDOMAIN_ACOUSTIC)
          write(IMAIN,*) '    material ',mat_id,' acoustic'
        case(IDOMAIN_ELASTIC)
          write(IMAIN,*) '    material ',mat_id,' elastic'
        case (IDOMAIN_POROELASTIC)
          write(IMAIN,*) '    material ',mat_id,' poroelastic'
        end select
        call flush_IMAIN()
      endif
    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  domain regions...'
      call flush_IMAIN()
    endif

    ! read number of subregions
    call read_value_integer_mesh(IIN,IGNORE_JUNK,NSUBREGIONS, 'NSUBREGIONS', ier)
    if (ier /= 0) stop 'Error reading Mesh parameter NSUBREGIONS'

    ! read subregions properties
    if (allocated(subregions)) deallocate(subregions)
    allocate(subregions(NSUBREGIONS,7),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1320')
    if (ier /= 0) stop 'Error allocation of subregions'
    subregions(:,:) = 0
    do ireg = 1,NSUBREGIONS
      call read_region_parameters(IIN,ix_beg_region,ix_end_region,iy_beg_region,iy_end_region, &
                                  iz_beg_region,iz_end_region,imaterial_number,ier)
      if (ier /= 0) then
        print *,'Error reading region ',ireg,' out of ',NSUBREGIONS
        stop 'Error reading regions in Mesh_Par_file'
      endif
      ! check for negative values: if ix or iy == -1, it means we take full range
      if (ix_end_region == -1) ix_end_region = NEX_XI
      if (iy_end_region == -1) iy_end_region = NEX_ETA

      ! stores region
      subregions(ireg,1) = ix_beg_region
      subregions(ireg,2) = ix_end_region
      subregions(ireg,3) = iy_beg_region
      subregions(ireg,4) = iy_end_region
      subregions(ireg,5) = iz_beg_region
      subregions(ireg,6) = iz_end_region
      subregions(ireg,7) = imaterial_number

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '    region ',ireg,' with material ',imaterial_number
        write(IMAIN,*) '      nex_xi  begin/end = ',ix_beg_region,ix_end_region
        write(IMAIN,*) '      nex_eta begin/end = ',iy_beg_region,iy_end_region
        write(IMAIN,*) '      nz      begin/end = ',iz_beg_region,iz_end_region
        call flush_IMAIN()
      endif
    enddo

    ! close parameter file
    call close_parameter_file_mesh()

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  reading Mesh_Par_file done successfully'
      write(IMAIN,*)
      write(IMAIN,*) '  checking mesh setup...'
      call flush_IMAIN()
    endif

    ! convert model size to UTM coordinates and depth of mesh to meters
    call utm_geo(LONGITUDE_MIN,LATITUDE_MIN,UTM_X_MIN,UTM_Y_MIN,ILONGLAT2UTM)
    call utm_geo(LONGITUDE_MAX,LATITUDE_MAX,UTM_X_MAX,UTM_Y_MAX,ILONGLAT2UTM)

    Z_DEPTH_BLOCK = - dabs(DEPTH_BLOCK_KM) * 1000.d0

    ! check that parameters computed are consistent
    if (UTM_X_MIN >= UTM_X_MAX) stop 'horizontal dimension of UTM block incorrect'
    if (UTM_Y_MIN >= UTM_Y_MAX) stop 'vertical dimension of UTM block incorrect'

    ! set time step and radial distribution of elements
    ! right distribution is determined based upon maximum value of NEX
    NEX_MAX = max(NEX_XI,NEX_ETA)
    UTM_MAX = max(UTM_Y_MAX-UTM_Y_MIN, UTM_X_MAX-UTM_X_MIN)/1000.0 ! in KM

    !------------------------------------
    ! Mesh_Par_file parameter check
    !------------------------------------

    ! counts defined/undefined materials
    ndef = 0
    nundef = 0
    do imat = 1,NMATERIALS
      mat_id = material_properties(imat,8)
      ! checks material id (must be positive or negative)
      if (mat_id == 0) stop 'Error incorrect material ID 0 found'
      ! counters
      if (mat_id > 0) ndef = ndef + 1
      if (mat_id < 0) nundef = nundef + 1
    enddo

    ! checks given material properties
    do imat = 1,NMATERIALS
      ! material properties
      rho = material_properties(imat,1)
      vp = material_properties(imat,2)
      vs = material_properties(imat,3)
      Q_Kappa = material_properties(imat,4)
      Q_mu = material_properties(imat,5)
      anisotropy_flag = material_properties(imat,6)
      domain_id = material_properties(imat,7)
      mat_id = material_properties(imat,8)

      ! checks material parameters
      if (rho <= 0.d0 .or. vp <= 0.d0 .or. vs < 0.d0) stop 'Error material with negative value of velocity or density found'
      if (Q_Kappa < 0.d0) stop 'Error material with negative Q_Kappa value found'
      if (Q_mu < 0.d0) stop 'Error material with negative Q_mu value found'

      ! checks domain id (1 = acoustic / 2 = elastic / 3 = poroelastic)
      select case (domain_id)
      case (IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC)
        continue
      case default
        stop 'Error invalid domain ID in materials'
      end select

      ! checks that material id doesn't exceed bounds of defined/undefined materials
      if (mat_id > 0 .and. mat_id > ndef) stop 'Error material ID bigger than total number of defined materials'
      if (mat_id < 0 .and. abs(mat_id) > nundef) stop 'Error negative material ID exceeds total number of undefined materials'
    enddo

    ! checks subregion ranges
    ner_max = 0
    do ireg = 1,NSUBREGIONS
      ix_beg_region = subregions(ireg,1)
      ix_end_region = subregions(ireg,2)
      iy_beg_region = subregions(ireg,3)
      iy_end_region = subregions(ireg,4)
      iz_beg_region = subregions(ireg,5)
      iz_end_region = subregions(ireg,6)
      imaterial_number = subregions(ireg,7)

      ! xi-direction range
      if (ix_beg_region < 1) stop 'XI coordinate of region negative!'
      if (ix_end_region > NEX_XI) stop 'XI coordinate of region too high!'

      ! eta-direction range
      if (iy_beg_region < 1) stop 'ETA coordinate of region negative!'
      if (iy_end_region > NEX_ETA) stop 'ETA coordinate of region too high!'

      ! depth range
      if (iz_beg_region < 1 .or. iz_end_region < 1) stop 'Z coordinate of region negative!'
      if (iz_end_region > ner_max) ner_max = iz_end_region

      if (imaterial_number == 0) stop 'Material ID of region zero!'
      ! searches material in given material_properties
      found = .false.
      do imat = 1,NMATERIALS
        mat_id = material_properties(imat,8)
        if (imaterial_number == mat_id) then
          found = .true.
          exit
        endif
      enddo
      if (.not. found) then
        print *,'Error: material id ',imaterial_number,' given in region ',ireg,' not found in materials section'
        stop 'Material ID of region not matching any given material'
      endif
    enddo
    ! checks if full range is covers
    if (minval(subregions(:,1)) > 1) stop 'Error minimum XI coordinate of regions must start at 1'
    if (maxval(subregions(:,2)) < NEX_XI) stop 'Error maximum XI coordinate of regions must end at NEX_XI'

    if (minval(subregions(:,3)) > 1) stop 'Error minimum ETA coordinate of regions must start at 1'
    if (maxval(subregions(:,4)) < NEX_ETA) stop 'Error maximum ETA coordinate of regions must end at NEX_ETA'

    if (minval(subregions(:,5)) > 1) stop 'Error minimum Z coordinate of regions must start at 1'
    ! maximum Z-layers must match with interface which is not read in here yet...
    !if (maxval(subregions(:,6)) < NER) stop 'Error maximum Z coordinate of regions must end at NER'

    ! checks doubling layer selection
    if (.not. USE_REGULAR_MESH) then
      do idoub = 1,NDOUBLINGS
        if (ner_doublings(idoub) < 1) then
          print *,'Error doubling layer ',idoub,' has invalid NZ value ',ner_doublings(idoub)
          stop 'Error doubling layer value must be at least 1'
        else if (ner_doublings(idoub) > ner_max) then
          print *,'Error doubling layer ',idoub,' with NZ value ',ner_doublings(idoub),'exceeds regions layers ',ner_max
          stop 'Error invalid doubling layer number, NZ exceeds regions NZ_END specification'
        endif
      enddo
    endif

    ! sorts doubling layer entries
    ! note: creating mesh layers starts from bottom, i.e. doubling builds layers from bottom to top,
    !       and a smaller NZ entry means the layer is closer to the bottom
    if (NDOUBLINGS > 1) then
      ! sorts with decreasing order, e.g. (14,10) instead of (10,14)
      ! switches doubling layer entries to have first doubling number with biggest and last entry with smallest value
      call bubble_sort_decreasing(NDOUBLINGS,ner_doublings)

      ! checks entries
      do idoub = 1,NDOUBLINGS-1
        ! checks if there are duplicate entries
        if (ner_doublings(idoub) == ner_doublings(idoub+1)) then
          print *,'Error doubling layer indices are invalid: ',ner_doublings(:)
          stop 'Error doubling layer entries contain duplicate values, please use unique layering'
        endif

        ! checks that entries are at least 2 layers apart
        ! (we alternate between regular and doubling subregions)
        if (ner_doublings(idoub) - ner_doublings(idoub+1) < 2 ) then
          print *,'Error doubling layer indices are too close: ',ner_doublings(idoub),' and ',ner_doublings(idoub+1)
          stop 'Error doubling layer entries must be at least 2 layers apart, please use larger layer spacings'
        endif
      enddo
    endif

    ! checks number of processes
    if (sizeprocs == 1 .and. (NPROC_XI /= 1 .or. NPROC_ETA /= 1)) &
      stop 'Error: must have NPROC_XI = NPROC_ETA = 1 for a serial run'

    ! checks CUBIT output
    if (SAVE_MESH_AS_CUBIT .and. NPROC_XI /= 1 .and. NPROC_ETA /= 1) &
      stop 'Error: SAVE_MESH_AS_CUBIT must have NPROC_XI = NPROC_ETA = 1'

    ! make sure everybody is synchronized
    call synchronize_all()

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  all okay'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  end subroutine read_mesh_parameter_file_fwat

  subroutine initialize_mesh_arrays()

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
    ! if (allocated(flag_sediments)) deallocate(flag_sediments)
    ! if (allocated(not_fully_in_bedrock)) deallocate(not_fully_in_bedrock)
    if (allocated(ibelm_xmin)) deallocate(ibelm_xmin)
    if (allocated(ibelm_xmax)) deallocate(ibelm_xmax)
    if (allocated(ibelm_ymin)) deallocate(ibelm_ymin)
    if (allocated(ibelm_ymax)) deallocate(ibelm_ymax)
    if (allocated(ibelm_bottom)) deallocate(ibelm_bottom)
    if (allocated(ibelm_top)) deallocate(ibelm_top)
    
  end subroutine initialize_mesh_arrays


end module meshfem3D_subs