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


  subroutine setup_sources_receivers_fwat(ievt)

  use specfem_par

  implicit none

  integer :: ievt

  ! checks if anything to do
  if (.not. IO_compute_task) return

  ! setup for point search
  call setup_point_search_arrays()

  ! sets number of timestep parameters
  call setup_timesteps()

  ! locates sources and determines simulation start time t0
  call setup_sources_fwat(ievt)

  ! reads in stations file and locates receivers
  call setup_receivers_fwat(ievt)

  ! pre-compute source arrays
  call setup_sources_precompute_arrays()

  ! pre-compute receiver interpolation factors
  call setup_receivers_precompute_intp()

  ! write source and receiver VTK files for Paraview
  ! call setup_sources_receivers_VTKfile()

  ! user output
  ! if (myrank == 0) then
  !   write(IMAIN,*)
  !   write(IMAIN,*) 'Total number of samples for seismograms = ',NSTEP
  !   write(IMAIN,*)
  !   if (NSOURCES > 1) then
  !     write(IMAIN,*) 'Using ',NSOURCES,' point sources'
  !     write(IMAIN,*)
  !   endif
  !   call flush_IMAIN()
  ! endif

  ! setup for coupling boundary points
  call setup_coupling_boundary_points()

  ! ! frees memory
  ! if (.not. INVERSE_FWI_FULL_PROBLEM) then
    ! deallocate(xyz_midpoints)
    ! deallocate(anchor_iax,anchor_iay,anchor_iaz)
  !   ! no need to keep mesh adjacency after point searches
    ! deallocate(neighbors_xadj,neighbors_adjncy)
  ! endif
  if (.not. DO_BRUTE_FORCE_POINT_SEARCH) call setup_free_kdtree()

  ! synchronizes processes
  call synchronize_all()

  end subroutine setup_sources_receivers_fwat

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_fwat(ievt)

  use specfem_par

  implicit none

  ! local parameters
  integer :: isource,ier, ievt

  ! initializes onset time (depends on hdur and source time function)
  t0 = 0.d0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'sources: ',NSOURCES
    ! finite fault simulations ignore CMT and force sources
    if (HAS_FINITE_FAULT_SOURCE) &
      write(IMAIN,*) '         finite fault simulation'
    call flush_IMAIN()
  endif

  ! allocate arrays for source
  allocate(islice_selected_source(NSOURCES), &
           ispec_selected_source(NSOURCES), &
           Mxx(NSOURCES), &
           Myy(NSOURCES), &
           Mzz(NSOURCES), &
           Mxy(NSOURCES), &
           Mxz(NSOURCES), &
           Myz(NSOURCES), &
           xi_source(NSOURCES), &
           eta_source(NSOURCES), &
           gamma_source(NSOURCES), &
           tshift_src(NSOURCES), &
           hdur(NSOURCES), &
           hdur_Gaussian(NSOURCES), &
           utm_x_source(NSOURCES), &
           utm_y_source(NSOURCES), &
           nu_source(NDIM,NDIM,NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2049')
  if (ier /= 0) stop 'error allocating arrays for sources'

  ! initializes arrays
  islice_selected_source(:) = -1
  ispec_selected_source(:) = 0
  Mxx(:) = 0.d0; Myy(:) = 0.d0; Mzz(:) = 0.d0
  Mxy(:) = 0.d0; Mxz(:) = 0.d0; Myz(:) = 0.d0
  xi_source(:) = 0.d0; eta_source(:) = 0.d0; gamma_source(:) = 0.d0
  tshift_src(:) = 0.d0; hdur(:) = 0.d0; hdur_Gaussian(:) = 0.d0
  utm_x_source(:) = 0.d0; utm_y_source(:) = 0.d0
  nu_source(:,:,:) = 0.d0

  if (USE_FORCE_POINT_SOURCE) then
    allocate(force_stf(NSOURCES), &
             factor_force_source(NSOURCES), &
             comp_dir_vect_source_E(NSOURCES), &
             comp_dir_vect_source_N(NSOURCES), &
             comp_dir_vect_source_Z_UP(NSOURCES),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2053')
  else
    allocate(force_stf(1), &
             factor_force_source(1), &
             comp_dir_vect_source_E(1), &
             comp_dir_vect_source_N(1), &
             comp_dir_vect_source_Z_UP(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2057')
  endif
  if (ier /= 0) stop 'error allocating arrays for force point sources'
  force_stf(:) = 0
  factor_force_source(:) = 0.d0
  comp_dir_vect_source_E(:) = 0.d0
  comp_dir_vect_source_N(:) = 0.d0
  comp_dir_vect_source_Z_UP(:) = 0.d0

  ! sets the size of user_source_time_function array
  if (USE_EXTERNAL_SOURCE_FILE) then
    NSTEP_STF = NSTEP
    NSOURCES_STF = NSOURCES
  else
    ! We don't need the array user_source_time_function : use a small dummy array
    NSTEP_STF = 1
    NSOURCES_STF = 1
  endif
  ! allocate array that contains the user defined source time function
  allocate(user_source_time_function(NSTEP_STF, NSOURCES_STF),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2058')
  if (ier /= 0) stop 'error allocating arrays for user sources time function'
  user_source_time_function(:,:) = 0.0_CUSTOM_REAL

  ! fused wavefield simulations
  call get_run_number_of_the_source()

  ! checks if anything left to do
  if (INVERSE_FWI_FULL_PROBLEM) then
    ! fwi will determine acquisition receivers later
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'running inverse FWI full problem, will determine sources later...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    ! all done
    return
  endif

  ! locate sources in the mesh
  call locate_source_fwat(ievt)

  ! count number of sources located in this slice
  nsources_local = 0
  do isource = 1, NSOURCES
    if (myrank == islice_selected_source(isource)) nsources_local = nsources_local + 1
  enddo

  ! checks if source is in an acoustic element and exactly on the free surface because pressure is zero there
  call setup_sources_check_acoustic()

  ! determines onset time
  call setup_stf_constants()

  ! prints source time functions to output files
  if (PRINT_SOURCE_TIME_FUNCTION) call print_stf_file()

  end subroutine setup_sources_fwat
!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers_fwat(ievt)

  use specfem_par
  use specfem_par_acoustic
  use input_params, only: fpar => fwat_par_global

  implicit none

  ! local parameters
  integer :: nrec_simulation,nrec_filtered, ievt
  integer :: nrec_tot_found
  integer :: irec,isource,ier
  character(len=MAX_STRING_LEN) :: rec_filename,filtered_rec_filename
  character(len=MAX_STRING_LEN) :: path_to_add

  ! adjoint sources/receivers
  integer :: icomp,itime,nadj_files_found,nadj_files_found_tot
  real(kind=CUSTOM_REAL) :: junk
  character(len=3),dimension(NDIM) :: comp
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: adj_source_file

  ! initializes
  nrec = 0
  nrec_filtered = 0
  nrec_local = 0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'receivers:'
    call flush_IMAIN()
  endif

  ! gets number of stations
  if (INVERSE_FWI_FULL_PROBLEM) then
    ! stations will be determined later...
    nrec = 0
    nrec_local = 0
  else
    ! reads in station file
    if (SIMULATION_TYPE == 1) then
      ! rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS'
      rec_filename = fpar%acqui%station_file(ievt)
      ! filtered_rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_FILTERED'
      filtered_rec_filename = trim(fpar%acqui%station_file(ievt))//'_FILTERED'
    else
      ! rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_ADJOINT'
      rec_filename = trim(fpar%acqui%station_file(ievt))//'_ADJOINT'
      ! filtered_rec_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_ADJOINT_FILTERED'
      filtered_rec_filename = trim(fpar%acqui%station_file(ievt))//'_ADJOINT_FILTERED'
    endif

    ! see if we are running several independent runs in parallel
    ! if so, add the right directory for that run
    ! (group numbers start at zero, but directory names start at run0001, thus we add one)
    ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      rec_filename = path_to_add(1:len_trim(path_to_add))//rec_filename(1:len_trim(rec_filename))
      filtered_rec_filename = path_to_add(1:len_trim(path_to_add))//filtered_rec_filename(1:len_trim(filtered_rec_filename))
    endif

    call station_filter(rec_filename,filtered_rec_filename,nrec_filtered)

    ! sets actual number of stations
    nrec = nrec_filtered

    if (nrec < 1) call exit_MPI(myrank,'need at least one receiver')
  endif

  call synchronize_all()

  ! user output
  if (myrank == 0) then
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      write(IMAIN,*) '  Total number of receivers = ', nrec
    else
      write(IMAIN,*) '  Total number of adjoint sources = ', nrec
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocate memory for receiver arrays, i.e. stations given in STATIONS file
  allocate(islice_selected_rec(nrec), &
           ispec_selected_rec(nrec), &
           xi_receiver(nrec), &
           eta_receiver(nrec), &
           gamma_receiver(nrec), &
           nu_rec(NDIM,NDIM,nrec), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2067')
  if (ier /= 0) stop 'error allocating arrays for receivers'
  ! initializes arrays
  islice_selected_rec(:) = -1
  ispec_selected_rec(:) = 0
  xi_receiver(:) = 0.d0; eta_receiver(:) = 0.d0; gamma_receiver(:) = 0.d0
  nu_rec(:,:,:) = 0.0d0

  allocate(station_name(nrec), &
           network_name(nrec), &
           stlat(nrec), &
           stlon(nrec), &
           stele(nrec), &
           stbur(nrec),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating receiver arrays')
  stlat(:) = 0.d0; stlon(:) = 0.d0; stele(:) = 0.d0; stbur(:) = 0.d0
  station_name(:) = ""; network_name(:) = ""

  ! checks if anything left to do
  if (INVERSE_FWI_FULL_PROBLEM) then
    ! fwi will determine acquisition receivers later
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'running inverse FWI full problem, will determine receivers later...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    ! all done
    return
  endif

  ! locate receivers in the mesh
  call locate_receivers(filtered_rec_filename,utm_x_source(1),utm_y_source(1))

  ! count number of receivers located in this slice
  nrec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! number of receivers are given by stations
    ! in STATIONS (forward runs) or STATIONS_ADJOINT (kernel runs) file
    nrec_simulation = nrec
    do irec = 1,nrec
      if (myrank == islice_selected_rec(irec)) nrec_local = nrec_local + 1
    enddo
  else
    ! adjoint simulation:
    ! station locations (in STATIONS_ADJOINT file) become adjoint sources
    ! and source locations (in CMTSOLUTION file) become adjoint "receivers"
    nrec_simulation = NSOURCES
    do isource = 1, NSOURCES
      if (myrank == islice_selected_source(isource)) nrec_local = nrec_local + 1
    enddo
  endif

  ! check that the sum of the number of receivers in each slice is nrec
  call sum_all_i(nrec_local,nrec_tot_found)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'found a total of ',nrec_tot_found,' receivers in all the slices'
    write(IMAIN,*)
  endif

  ! checks
  if (myrank == 0) then
    if (nrec_tot_found /= nrec_simulation) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    endif
  endif
  call synchronize_all()

  ! checks if acoustic receiver is exactly on the free surface because pressure is zero there
  call setup_receivers_check_acoustic()

  ! counter for adjoint receiver stations in local slice, used to allocate adjoint source arrays
  nadj_rec_local = 0

  ! counts receivers for adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! ADJOINT simulations
    ! gets number of local adjoint sources, i.e. located in this slice (nadj_rec_local)
    if (.not. SU_FORMAT) then
      ! prepares channel names
      do icomp = 1,NDIM
        call write_channel_name(icomp,comp(icomp))
      enddo

      ! temporary counter to check if any files are found at all
      nadj_files_found = 0

      do irec = 1,nrec
        ! adjoint receiver station
        if (myrank == islice_selected_rec(irec)) then
          ! checks that the source slice number is okay
          if (islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROC-1) &
            call exit_MPI(myrank,'something is wrong with the source slice number in adjoint simulation')

          ! updates counter
          nadj_rec_local = nadj_rec_local + 1

          ! checks **net**.**sta**.**MX**.adj files for correct number of time steps
          if (READ_ADJSRC_ASDF) then
            ! ASDF format
            call check_adjoint_sources_asdf(irec,nadj_files_found)
          else
            ! ASCII format
            adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))

            ! loops over file components E/N/Z
            do icomp = 1,NDIM
              ! trace name
              ! format: **net**.**sta**.**comp**.adj
              !         example: DB.X01.HXX.adj
              filename = trim(adj_source_file) // '.' // comp(icomp) // '.adj'

              ! adds path to SEM/ folder relative to OUTPUT_FILES/
              filename = OUTPUT_FILES(1:len_trim(OUTPUT_FILES))// '/../SEM/' // trim(filename)

              ! reads in adjoint trace
              open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
              if (ier == 0) then
                ! checks length of file
                itime = 0
                do while (ier == 0)
                  read(IIN,*,iostat=ier) junk,junk
                  if (ier == 0) itime = itime + 1
                enddo

                ! checks length
                if (itime /= NSTEP) then
                  print *,'adjoint source error: ',trim(filename),' has length',itime,' but should be',NSTEP
                  call exit_MPI(myrank, &
                    'file '//trim(filename)//' has wrong length, please check your adjoint sources and your simulation duration')
                endif

                ! updates counter for found files
                nadj_files_found = nadj_files_found + 1
              endif
              ! closes file
              close(IIN)
            enddo
          endif
        endif
      enddo

    else
      ! SU_FORMAT file
      ! adjoint sources given in single SU_FORMAT file;
      ! skip counting, because only one file per component per proc in SU_FORMAT
      nadj_rec_local = nrec_local
      nadj_files_found = nrec_local
    endif !if (.not. SU_FORMAT)

    ! checks if any adjoint source files found at all
    call sum_all_i(nadj_files_found,nadj_files_found_tot)
    if (myrank == 0) then
      ! user output
      write(IMAIN,*) '    ',nadj_files_found_tot,' adjoint component trace files found in all slices'
      write(IMAIN,*)
      call flush_IMAIN()

      ! main process checks if any adjoint files found
      if (nadj_files_found_tot == 0) then
        print *,'Error no adjoint traces found: ',nadj_files_found_tot
        print *,'in directory : ',OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//'/../SEM/'
        if (.not. SU_FORMAT .and. .not. READ_ADJSRC_ASDF) then
          print *,'with endings : ', '**.'//comp(1)//'.adj',' ','**.'//comp(2)//'.adj',' ','**.'//comp(3)//'.adj'
        endif
        print *
        call exit_MPI(myrank,'no adjoint traces found, please check adjoint sources in directory SEM/')
      endif
    endif

    ! initializes adjoint sources
    allocate(source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2073')
    if (ier /= 0) stop 'error allocating array source_adjoint'
    source_adjoint(:,:,:) = 0.0_CUSTOM_REAL

    ! note:
    ! computes adjoint sources in chunks/blocks during time iterations.
    ! we moved it to compute_add_sources_viscoelastic.f90 & compute_add_sources_acoustic.f90,
    ! because we may need to read in adjoint sources block by block
  endif

  end subroutine setup_receivers_fwat


!
