program xmeasure_adj
  use measure_adj_mod
  use specfem_par, only: OUTPUT_FILES,MAX_STRING_LEN, SPECFEM_T0 => T0 
  use shared_input_parameters, only: NSTEP, SPECFEM_DT => DT 

  character(len=512)                               :: datafile,synfile,synfile_phydisp,file_prefix0,bandname
  double precision, dimension(NDIM)                :: data,syn,syn_phydisp,adj_syn_all
  integer                                          :: npt1,npt2,npt3,npts,nn
  double precision                                 :: t0,t01,t02,t03,dt,dt1,dt2,dt3,tt,dtt
  double precision                                 :: tstart,tend,fend0, fstart0
  double precision, dimension(100)                 :: tstart_arr,tend_arr
  integer                                          :: yr,jda,ho,mi
  integer, dimension(100)                          :: out_imeas

  double precision                                 :: sec,dist,az,baz,slat,slon
  character(len=10)                                :: net,sta,chan_dat,chan
  double precision, dimension(NCHI,100)            :: window_chi
  double precision, dimension(100)                 :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn
  integer                                          :: ios, j, ipair, npairs, nwin, num_meas
  double precision                                 :: T_surfacewaves
  
  !************* use origianl code for benchmark ***********************
  !call measure_adj() ! need MESUREMENT.PAR and MEASUREMENT.WINDOWS
  
  !******************The following are commented out in measure_adj_fwat************************************
  ! uses weights to balance love and rayleigh measurements
  ! we do a normalization of P_SV, P_SH, Love, Rayleigh with the number of measurement picks
  !if( DO_WEIGHTING ) call setup_weighting(chan)

  ! input file: MEASUREMENT.WINDOWS
  open(11,file='MEASUREMENT.WINDOWS',status='old',iostat=ios)
  if (ios /= 0) stop 'Error opening input file: MEASUREMENT WINDOWS'

  read(11,*,iostat=ios) npairs
  if (ios /= 0) stop 'Error reading number of pairs of data/syn'
  print *, 'reading in the data and synthetics...'

  ! output files
  open(12,file='window_index',status='unknown',iostat=ios)
  open(13,file='window_chi',status='unknown',iostat=ios)
  call read_parameter_file(0, .false.)
  call read_par_file(fstart0,fend0,tt,dtt,nn,chan)
  SPECFEM_DT = dtt
  NSTEP = nn
  SPECFEM_T0 = tt

  nwin = 0;
  do ipair = 1, npairs

    data(:) = 0.0
    syn(:)  = 0.0
    ! reads in file names for data and synthetics
    read(11,'(a)',iostat=ios) datafile
    if (ios /= 0) stop 'Error reading windows file'
    read(11,'(a)',iostat=ios) synfile
    if (ios /= 0) stop 'Error reading windows file'
    if (USE_PHYSICAL_DISPERSION) then
            synfile_phydisp=trim(synfile)//'.phydisp'
    endif

    ! read data and syn (in double precision)
    ! LQY: data is read last to be stored in memory for header retrieval later

    call drsac1(datafile,data,npt1,t01,dt1)
    call drsac1(synfile,syn,npt2,t02,dt2)
    if (USE_PHYSICAL_DISPERSION) then
            call drsac1(synfile_phydisp,syn_phydisp,npt3,t03,dt3)
    endif

    if (DISPLAY_DETAILS) then
       print *
       print *, 'data: ',trim(datafile)
       print *, '  min/max: ',sngl(minval(data(:))),sngl(maxval(data(:)))

       print *, 'syn:   ',trim(synfile)
       print *, '  min/max: ',sngl(minval(syn(:))),sngl(maxval(syn(:)))
    endif

    ! check if npts, dt, t0 matches for data and synthetics
    ! no interpolation is necessary at any point
    if (max(npt1,npt2) > NDIM) &
         stop 'Error: Too many number of points in data or syn'
    npts = min(npt1,npt2)

    if (abs(dt1-dt2) > TOL) stop 'Error: check if dt match'
    dt = dt1

    if (abs(t01-t02) > dt)  stop 'Check if t0 match'
    t0 = t01

    if (DISPLAY_DETAILS) print *,'  time, dt, npts :',sngl(t01), sngl(dt), npts


    ! apply bandpass filter to data and synthetics with saclib, if desired
    ! http://www.iris.washington.edu/pipermail/sac-help/2008-March/000376.html
    ! Access to the kidate, xapiir, and getfil is not simple and not
    ! supported under the current state of the SAC code base.

    ! if(RUN_BANDPASS) then
    !    call bandpass(data,npts,dt,fstart0,fend0)
    !    call bandpass(syn,npts,dt,fstart0,fend0)
    !    if (USE_PHYSICAL_DISPERSION) then
    !            call bandpass(syn_phydisp,npts,dt,fstart0,fend0)
    !    endif
    ! endif

    ! find out station/network/comp names,etc from synthetics
    call get_sacfile_header(trim(synfile),yr,jda,ho,mi,sec,net,sta, &
                          chan_dat,dist,az,baz,slat,slon)

    ! theoretical surface wave arrival time
    T_surfacewaves = dist / surface_vel

    ! reads number of measurement windows
    read(11,*,iostat=ios) num_meas
    if (ios /= 0) stop 'Error reading num_meas'

    do j = 1, num_meas
      ! reads in start and end time of the measurement window
      read(11,*,iostat=ios) tstart, tend
      if (ios /= 0) stop 'Error reading tstart and tend'
      tstart_arr(j)=tstart
      tend_arr(j)=tend
    enddo ! end j

    !   ! checks start and end times of window compared to trace lengths
    !   tstart = max(tstart,t0)
    !   tend = min(tend, t0+(npts-1)*dt)
    !   nlen = floor((tend-tstart)/dt) + 1  ! dummy, replaced later in mt_measure()

    !   ! write values to output file
    !   nwin = nwin + 1       ! overall window counter
    !   write(12,'(a3,a8,a5,a5,3i5,2f12.3)') net,sta,chan_syn,chan_dat,nwin,ipair,j,tstart,tend
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      OUT_DIR=OUTPUT_FILES
      bandname='T006_T030'
      call measure_adj_fwat(data,syn,tstart_arr(1),tend_arr(1),t0,dt,npts,net,sta,&
                            chan_dat,window_chi(:,1),tr_chi(1),am_chi(1),&
                            T_pmax_dat(1),T_pmax_syn(1),adj_syn_all,file_prefix0,out_imeas(1),bandname)
    do j = 1, num_meas
      write(13,'(a14,a8,a3,a5,i4,i4,2e14.6,20e14.6,2e14.6,2f14.6)') &
               file_prefix0,sta,net,chan_dat,j,out_imeas(j),&
               tstart_arr(j),tend_arr(j),window_chi(:,j),tr_chi(j),am_chi(j),T_pmax_dat(j),T_pmax_syn(j)
    enddo
  enddo ! end ipair



  close(12)
  close(13)
 
end program xmeasure_adj


  subroutine usage()

  print *, 'Usage: xmeasure_adj datafile synfile bandname tb te'
  stop

  end subroutine usage
