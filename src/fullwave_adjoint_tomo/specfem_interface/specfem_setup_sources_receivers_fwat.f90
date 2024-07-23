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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine setup_sources_receivers_fwat(source_fname,station_fname)

  use specfem_par

  use kdtree_search, only: kdtree_delete,kdtree_nodes_location,kdtree_nodes_index

  implicit none
  character(len= MAX_STRING_LEN)              :: source_fname
  character(len= MAX_STRING_LEN)              :: station_fname

  ! builds search tree
  if (.not. DO_BRUTE_FORCE_POINT_SEARCH) call setup_search_kdtree()

  ! sets number of timestep parameters
  call setup_timesteps()

  ! locates sources and determines simulation start time t0
  call setup_sources_fwat(source_fname)

  ! reads in stations file and locates receivers
  call setup_receivers_fwat(station_fname)

  ! pre-compute source arrays
  call setup_sources_precompute_arrays()

  ! pre-compute receiver interpolation factors
  call setup_receivers_precompute_intp()

  ! write source and receiver VTK files for Paraview
  call setup_sources_receivers_VTKfile()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of samples for seismograms = ',NSTEP
    write(IMAIN,*)
    write(IMAIN,*) 'found a total of ',nrec_tot_found,' receivers in all the slices'
    write(IMAIN,*)
    if (NSOURCES > 1) then
      write(IMAIN,*) 'Using ',NSOURCES,' point sources'
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! frees tree memory
  if (.not. DO_BRUTE_FORCE_POINT_SEARCH) then
    ! deletes tree arrays
    deallocate(kdtree_nodes_location)
    deallocate(kdtree_nodes_index)
    ! deletes search tree nodes
    call kdtree_delete()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine setup_sources_receivers_fwat

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_fwat(source_fname)

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  implicit none

  double precision :: min_tshift_src_original
  integer :: isource,ier
  character(len=MAX_STRING_LEN) :: SOURCE_FILE,path_to_add
  character(len= MAX_STRING_LEN)              :: source_fname

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'sources:',NSOURCES
    call flush_IMAIN()
  endif

  ! allocate arrays for source
  allocate(islice_selected_source(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2033')
  allocate(ispec_selected_source(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2034')
  allocate(Mxx(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2035')
  allocate(Myy(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2036')
  allocate(Mzz(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2037')
  allocate(Mxy(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2038')
  allocate(Mxz(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2039')
  allocate(Myz(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2040')
  allocate(xi_source(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2041')
  allocate(eta_source(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2042')
  allocate(gamma_source(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2043')
  allocate(tshift_src(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2044')
  allocate(hdur(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2045')
  allocate(hdur_Gaussian(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2046')
  allocate(utm_x_source(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2047')
  allocate(utm_y_source(NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2048')
  allocate(nu_source(NDIM,NDIM,NSOURCES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2049')
  if (ier /= 0) stop 'error allocating arrays for sources'

  if (USE_FORCE_POINT_SOURCE) then
    allocate(factor_force_source(NSOURCES),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2050')
    allocate(comp_dir_vect_source_E(NSOURCES),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2051')
    allocate(comp_dir_vect_source_N(NSOURCES),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2052')
    allocate(comp_dir_vect_source_Z_UP(NSOURCES),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2053')
  else
    allocate(factor_force_source(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2054')
    allocate(comp_dir_vect_source_E(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2055')
    allocate(comp_dir_vect_source_N(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2056')
    allocate(comp_dir_vect_source_Z_UP(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2057')
  endif
  if (ier /= 0) stop 'error allocating arrays for force point sources'

  !! allocate the array contains the user defined source time function
  allocate(user_source_time_function(NSTEP_STF, NSOURCES_STF),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2058')
  if (ier /= 0) stop 'error allocating arrays for user sources time function'

  if (USE_FORCE_POINT_SOURCE) then
    SOURCE_FILE = trim(source_fname)!IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'FORCESOLUTION'
  else
    SOURCE_FILE = trim(source_fname)!IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'CMTSOLUTION'
  endif

  ! see if we are running several independent runs in parallel
  ! if so, add the right directory for that run
  ! (group numbers start at zero, but directory names start at run0001, thus we add one)
  ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    SOURCE_FILE = path_to_add(1:len_trim(path_to_add))//SOURCE_FILE(1:len_trim(SOURCE_FILE))
  endif

  ! locate sources in the mesh
  call locate_source(SOURCE_FILE,tshift_src,min_tshift_src_original,utm_x_source,utm_y_source, &
                     hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                     islice_selected_source,ispec_selected_source, &
                     factor_force_source,comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP, &
                     xi_source,eta_source,gamma_source,nu_source,user_source_time_function)

  ! determines onset time
  call setup_stf_constants(hdur,hdur_Gaussian,tshift_src,min_tshift_src_original, &
                           islice_selected_source,ispec_selected_source,t0)

  ! count number of sources located in this slice
  nsources_local = 0
  do isource = 1, NSOURCES
    if (myrank == islice_selected_source(isource)) nsources_local = nsources_local + 1
  enddo

  ! checks if source is in an acoustic element and exactly on the free surface because pressure is zero there
  call setup_sources_check_acoustic()

  ! prints source time functions to output files
  if (PRINT_SOURCE_TIME_FUNCTION) call print_stf_file()

  call get_run_number_of_the_source()

  end subroutine setup_sources_fwat
!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers_fwat(station_fname)

  use specfem_par
  use specfem_par_acoustic

  implicit none

  ! local parameters
  integer :: nrec_simulation
  integer :: irec,isource,ier
  character(len=MAX_STRING_LEN) :: path_to_add
  character(len= MAX_STRING_LEN)              :: station_fname

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'receivers:'
    call flush_IMAIN()
  endif

  ! reads in station file
  if (SIMULATION_TYPE == 1) then
    rec_filename = trim(station_fname)!IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS'
    filtered_rec_filename = trim(station_fname)//'_FILTERED'!IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_FILTERED'
  else
    rec_filename = trim(station_fname)//'_ADJOINT'!IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_ADJOINT'
    filtered_rec_filename = trim(station_fname)//'_ADJOINT_FILTERED'!IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'STATIONS_ADJOINT_FILTERED'
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

  call station_filter(rec_filename,filtered_rec_filename,nrec)

  if (nrec < 1) call exit_MPI(myrank,'need at least one receiver')
  call synchronize_all()

  if (myrank == 0) then
    write(IMAIN,*)
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      write(IMAIN,*) 'Total number of receivers = ', nrec
    else
      write(IMAIN,*) 'Total number of adjoint sources = ', nrec
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocate memory for receiver arrays, i.e. stations given in STATIONS file
  allocate(islice_selected_rec(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2060')
  allocate(ispec_selected_rec(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2061')
  allocate(xi_receiver(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2062')
  allocate(eta_receiver(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2063')
  allocate(gamma_receiver(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2064')
  allocate(station_name(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2065')
  allocate(network_name(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2066')
  allocate(nu(NDIM,NDIM,nrec), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2067')
  if (ier /= 0) stop 'error allocating arrays for receivers'

  ! locate receivers in the mesh
  call locate_receivers(filtered_rec_filename,nrec,islice_selected_rec,ispec_selected_rec, &
                        xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
                        utm_x_source(1),utm_y_source(1))

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
  if (myrank == 0) then
    if (nrec_tot_found /= nrec_simulation) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    endif
  endif
  call synchronize_all()

! checks if acoustic receiver is exactly on the free surface because pressure is zero there
  call setup_receivers_check_acoustic()

end subroutine setup_receivers_fwat


subroutine read_receiver_file(filename, stla, stlo, stel)
  use specfem_par, only : nrec, IIN, MAX_STRING_LEN, CUSTOM_REAL, myrank
  use fwat_utils, only : append

  implicit none

  character(len=*) :: filename
  character(len=MAX_STRING_LEN) :: dummystring, station_name, network_name
  real(kind=CUSTOM_REAL) :: junk_lat, junk_lon, junk_ele, junk_bur
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: stla, stlo, stel 
  integer :: ier
  integer :: i, this_nrec

  if (myrank == 0) then

    open(unit=IIN, file=trim(filename), status='old', iostat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'No file '//trim(filename)//', exit')
    
    do 
      read(IIN,*,iostat=ier) dummystring
      if (ier /= 0) exit
      if (len_trim(dummystring) > 0) then
        read(IIN,*,iostat=ier) station_name, network_name, junk_lat, junk_lon, junk_ele, junk_bur
        call append(stla, junk_lat)
        call append(stlo, junk_lon)
        call append(stel, junk_ele)
      endif
    enddo
    close (IIN)
  endif
  this_nrec = size(stla)
  call bcast_all_singlei(this_nrec)
  if (myrank /= 0) then
    allocate(stla(this_nrec), stlo(this_nrec), stel(this_nrec))
  endif
  call bcast_all_cr(stla, this_nrec)
  call bcast_all_cr(stlo, this_nrec)
  call bcast_all_cr(stel, this_nrec)

end subroutine read_receiver_file