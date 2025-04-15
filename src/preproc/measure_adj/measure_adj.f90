module measure_adj_mod
  use ma_variables
  use ma_constants
  use ascii_rw       ! dwascii()
  use ma_sub2        ! fft(), fftinv()
  use ma_sub         ! mt_measure(), mt_adj()
  use ma_weighting
  use config, only: worldrank

  use specfem_par, only: OUTPUT_FILES,MAX_STRING_LEN, SPECFEM_T0 => T0,&
                         CUSTOM_REAL, NSTEP, SPECFEM_DT => DT 
  implicit none

contains

  subroutine measure_adj_fwat(data_in,syn_in,tstart,tend,shortp,longp,net_in,sta_in,&
                              chan_dat_in,window_chi,tr_chi,am_chi,&
                              T_pmax_dat,T_pmax_syn,adj_syn_all,file_prefix0,out_imeas)
  !  main program that calls the subroutines to make measurements within input time windows
  !  and then compute the corresponding adjoint sources

  ! input parameter:
  !  1. imeas = 1, normalized waveform difference. Adjoint source is constructed from the data
  !                      only, with the form -d(t)/ || d(t) || 2
  !  2. imeas = 2, waveform difference, s(t) - d(t).
  !  3. imeas = 3, cross-correlation traveltime difference for a (banana-doughtnut) sensitivity ker-
  !                       nel. The measurement between data and synthetics is not used in constructing the adjoint
  !                       source.
  !  4. imeas = 4, amplitude difference for a (banana-doughtnut) sensitivity kernel. The measure-
  !                       ment between data and synthetics is not used in constructing the adjoint source.
  !  5. imeas = 5, cross-correlation traveltime difference for an event kernel. The measurement
  !                       between data and synthetics is used in constructing the adjoint source.
  !  6. imeas = 6, amplitude difference for an event kernel. The measurement between data and
  !                        synthetics is used in constructing the adjoint source.
  !  7. imeas = 7, multitaper traveltime difference for an event kernel. The measurement between
  !                       data and synthetics is used in constructing the adjoint source. See multitaper_notes.pdf.
  !  8. imeas = 8, multitaper amplitude difference for an event kernel. The measurement between
  !                       data and synthetics is used in constructing the adjoint source. See multitaper_notes.pdf.

  implicit none

  character(len=150) :: datafile,file_prefix,file_prefix0,file_prefix2,measure_file_prefix,adj_file_prefix
  integer :: num_meas, j, npts, nn !, ios, npt3
  double precision, dimension(:), intent(in) :: data_in, syn_in
  character(len=*), intent(in) :: net_in,sta_in,chan_dat_in
  double precision, dimension(:), intent(inout) :: window_chi
  double precision, intent(in) :: shortp,longp, tstart, tend

  double precision, dimension(NDIM_MA) :: data, syn, syn_phydisp, adj_syn_all, &
                        tr_adj_src, am_adj_src, recon_cc_all, syn_dtw_cc, syn_dtw_mt
  double precision :: t0, dt, tt, dtt, df
  double precision :: fend0, fstart0, fend, fstart

  character(len=10) :: net,sta,chan_dat,chan,cmp,chan_syn
  double precision :: tshift, sigma_dt_cc, dlnA, sigma_dlnA_cc, sigma_dt, sigma_dlnA
  double precision :: all_chi, tr_chi, am_chi, cc_max, T_pmax_dat, T_pmax_syn
  !double precision :: tshift_f1f2, cc_max_f1f2
  double precision, dimension(NPT) :: dtau_w, dlnA_w, err_dt, err_dlnA, syn_dtw, data_dtw,syn_dtw_phydisp
  complex(kind=8), dimension(NPT) :: trans_mtm
  integer :: nlen, i_left, i_pmax_dat, i_pmax_syn, i_right, i_right0, istart, &
        ipair, nwin, itmax !, npairs
  logical :: use_trace
  !double precision :: trbdndw, a
  !integer :: iord, passes
  integer :: ipick_type
  double precision :: T_surfacewaves

  integer :: out_imeas
  character(len=MAX_STRING_LEN)                     :: bandname

  !********* PROGRAM STARTS HERE *********************
  !! read in MEASUREMENT.PAR (see ma_sub.f90 and write_par_file.pl)
  !! most parameters are global (see ma_variables.f90)
  ! call read_par_file(fstart0,fend0,tt,dtt,nn,chan)
  ! read(bandname(2:4),'(i3.3)') shortp
  ! read(bandname(7:9),'(i3.3)') longp
  fstart0=1.0/longp
  fend0=1.0/shortp
  TSHORT=shortp
  TLONG=longp

  t0 = -dble(SPECFEM_T0)
  dt = dble(SPECFEM_DT)
  npts = NSTEP
  data(1:npts) = data_in
  syn(1:npts) = syn_in
  net = trim(net_in)
  sta = trim(sta_in)
  chan_dat = trim(chan_dat_in)

  nwin = 0; all_chi=0.
  do ipair = 1, 1

    adj_syn_all(:) = 0.0
    recon_cc_all(:) = 0.0

  !  ! apply bandpass filter to data and synthetics with saclib, if desired
  !  ! http://www.iris.washington.edu/pipermail/sac-help/2008-March/000376.html
  !  ! Access to the kidate, xapiir, and getfil is not simple and not
  !  ! supported under the current state of the SAC code base.

    ! if(RUN_BANDPASS) then
    !    call bandpass(data,npts,dt,fstart0,fend0)
    !    call bandpass(syn,npts,dt,fstart0,fend0)
    !    if (USE_PHYSICAL_DISPERSION) then
    !            call bandpass(syn_phydisp,npts,dt,fstart0,fend0)
    !    endif
    ! endif
  !if (DISPLAY_DETAILS) then  
  !endif

    ! theoretical surface wave arrival time
    !T_surfacewaves = dist / surface_vel

    cmp = chan_dat(3:3)
    chan_syn = chan_dat

    ! example: OUT/PAS.CI.BHZ
    file_prefix0 = trim(sta)//'.'//trim(net)//'.'//trim(chan_syn)
    file_prefix2 = trim(OUT_DIR)//'/'//trim(file_prefix0)
    ! print *
    ! print *,  trim(file_prefix2), ' --- '

    ! note: MT measurement could revert to CC, but still keep the MT suffix
    !write(adj_file_prefix,'(a,i2.2)') trim(file_prefix2)//'.iker', imeas0
    adj_file_prefix = trim(net)//'.'//trim(sta)//'.'//trim(chan_syn)

    ! reads number of measurement windows
    !!!!!! WK added for noise/tele !!!!!!!
    num_meas=1 
    !!!!!!!!!!!!!!!!!!!!!!!

    !read(11,*,iostat=ios) num_meas
    !if (ios /= 0) stop 'Error reading num_meas'

    do j = 1, num_meas
      ! reads in start and end time of the measurement window
      !read(11,*,iostat=ios) tstart, tend
      !if (ios /= 0) stop 'Error reading tstart and tend'

      ! checks start and end times of window compared to trace lengths
      !tstart = max(tstart,t0)
      !tend = min(tend, t0+(npts-1)*dt)
      nlen = floor((tend-tstart)/dt) + 1  ! dummy, replaced later in mt_measure()

      ! write values to output file
      nwin = nwin + 1       ! overall window counter

      ! add taper type to file prefix: OUT/PAS.CI.BHZ.01.mtm
      write(file_prefix,'(a,i2.2)') trim(file_prefix2)//'.', j

      if (is_mtm == 1) then
        measure_file_prefix = trim(file_prefix) // '.mtm'  ! multitaper taper
      else if (is_mtm == 2) then
        measure_file_prefix = trim(file_prefix) // '.ctp'  ! cosine taper
      else
        measure_file_prefix = trim(file_prefix) // '.btp'  ! boxcar taper
      endif

      !print *
      !print *, ' Measurement window No.', j, ' ... '

      ! initialize the measurements

      window_chi(:) = 0.

      ! compute integrated waveform difference, normalized by duration of the record
      ! NOTE: (1) this is for the FULL record, not the windowed record
      !       (2) for comparison with waveform_chi, we include the 0.5 factor
      !       (3) we might want to include dt as an integration factor (also for waveform_chi),
      !           but the ratio (d-s)^2 / d^2 avoids the need for dt, nstep, or length of record
      window_chi(17) = 0.5 * sum( data**2 )
      window_chi(18) = 0.5 * sum( syn**2 )
      window_chi(19) = 0.5 * sum( (data-syn)**2 )
      window_chi(20) = npts*dt

      ! make measurements
      ! also compute reconstructed synthetics for CC (and MT, if specified) measurements
      T_pmax_dat = 0.
      T_pmax_syn = 0.
      call mt_measure(datafile,measure_file_prefix,data,syn,syn_phydisp,t0,dt,npts,tstart,tend,&
            istart,data_dtw,syn_dtw,syn_dtw_phydisp,nlen,tshift,sigma_dt_cc,dlnA,sigma_dlnA_cc,cc_max,syn_dtw_cc,&
            i_pmax_dat,i_pmax_syn,i_right0,trans_mtm,dtau_w,dlnA_w,sigma_dt,sigma_dlnA,syn_dtw_mt,err_dt,err_dlnA)
      i_right = i_right0
      i_left = 1  ! LQY: is it feasible that i_left is not 1? mt_adj() inherently assumes it.

      ! adjust frequency ranges for MT measurements
      ! fstart is constrained by NCYCLE_IN_WINDOW/tlen, fend constrained by i_right
      if (is_mtm == 1) then
        ! period of the max power of the synthetic record
        T_pmax_dat = (dt*NPT) / dble(i_pmax_dat)
        T_pmax_syn = (dt*NPT) / dble(i_pmax_syn)

        fstart = fstart0  ; fend = fend0
        call mt_measure_select(nlen,tshift,i_pmax_syn,dtau_w,err_dt, &
                            dt,i_left,i_right,fstart,fend,use_trace)
        !print *, '     Tlong/Tshort (input) :', sngl(1/fstart0), sngl(1/fend0)
        !print *, '     Tlong/Tshort (adjusted)  :', sngl(1/fstart), sngl(1/fend)
        !print *, '     period of max data/syn power    :', sngl(T_pmax_dat), sngl(T_pmax_syn)

        ! if MT measurement window is rejected by mt_measure_select, then use a CC measurement
        if(.not. use_trace) then
          !stop 'Check why this MT measurement was rejected'
          ! print *, '   reverting from MT measurement to CC measurement...'
          imeas = imeas0 - 2
          is_mtm = 3  ! LQY: WHY not is_mtm = 2?
          call mt_measure(datafile,measure_file_prefix,data,syn,syn_phydisp,t0,dt,npts,tstart,tend,&
                istart,data_dtw,syn_dtw,syn_dtw_phydisp,nlen,&
                tshift,sigma_dt_cc,dlnA,sigma_dlnA_cc,cc_max,syn_dtw_cc,&
                i_pmax_dat,i_pmax_syn,i_right,trans_mtm,dtau_w,dlnA_w,sigma_dt,sigma_dlnA,syn_dtw_mt)
        ! else
        !   ! print *, '     using this MTM. '
        endif
      endif

      out_imeas=imeas
      ! check that the CC measurements are within the specified input range
      if (imeas >= 5) call cc_measure_select(tshift,dlnA,cc_max)

      ! ! write frequency limits to file
      ! if (OUTPUT_MEASUREMENT_FILES) then
      !   df = 1./(dt*NPT)
      !   open(71,file=trim(measure_file_prefix)//'.freq_limits')
      !   write(71,'(6f18.8)') fstart0, fend0, df, i_right0*df, fstart, fend
      !   close(71)
      ! endif

      ! compute adjoint sources and misfit function values and also the CC-reconstructed records
      if (COMPUTE_ADJOINT_SOURCE) then
         !print *, '   Generating adjoint source and chi value for imeas = ', imeas

        ! banana-doughnut kernel (needs only synthetic trace)
        ! LQY: what is this section intended to do?
        ! reset imeas == 3 for adjoint sources without time shift and uncertainty scaling
        ! (pure cross-correlation adjoint source for banana-doughnuts)
        !if(imeas == 5 .and. trim(datafile) == trim(synfile) ) then
        !   print *,'cross-correlation measurement:'
        !   print *,'  only synthetic file: ',trim(synfile)
        !   print *,'    without traveltime difference/uncertainty'
        !   print *
        !   imeas = 3
        !endif

        tr_chi = 0.0 ; am_chi = 0.0    ! must be initialized
        call mt_adj(istart,data_dtw,syn_dtw,syn_dtw_phydisp,nlen,dt,tshift,dlnA,sigma_dt_cc,sigma_dlnA_cc,&
             dtau_w,dlnA_w,err_dt,err_dlnA,sigma_dt,sigma_dlnA,i_left,i_right,&
             window_chi,tr_adj_src,tr_chi,am_adj_src,am_chi)

        ! KEY: write misfit function values to file (two for each window)
        ! Here are the 20 columns of the vector window_chi
        !  1: MT-TT chi,    2: MT-dlnA chi,    3: XC-TT chi,    4: XC-dlnA chi
        !  5: MT-TT meas,   6: MT-dlnA meas,   7: XC-TT meas,   8: XC-dlnA meas
        !  9: MT-TT error, 10: MT-dlnA error, 11: XC-TT error, 12: XC-dlnA error
        ! WINDOW     : 13: data power, 14: syn power, 15: (data-syn) power, 16: window duration
        ! FULL RECORD: 17: data power, 18: syn power, 19: (data-syn) power, 20: record duration
        ! Example of a reduced file: awk '{print $2,$3,$4,$5,$6,$31,$32}' window_chi > window_chi_sub
        !write(13,'(a14,a8,a3,a5,i4,i4,2e14.6,20e14.6,2e14.6,2f14.6)') &
        !   file_prefix0,sta,net,chan_syn,j,imeas,&
        !   tstart,tend,window_chi(:),tr_chi,am_chi,T_pmax_dat,T_pmax_syn
        !print *, '     tr_chi = ', sngl(tr_chi), '  am_chi = ', sngl(am_chi)

        ! uses weighting to balance love / rayleigh measurements
        if( DO_WEIGHTING ) then
           ipick_type = 0
           if( tend <= T_surfacewaves ) then
              ! body wave picks
              if( cmp(1:1) == "Z" ) ipick_type = P_SV_V
              if( cmp(1:1) == "R" ) ipick_type = P_SV_R
              if( cmp(1:1) == "T" ) ipick_type = SH_T
           else
              ! surface wave picks
              if( cmp(1:1) == "Z" ) ipick_type = Rayleigh_V
              if( cmp(1:1) == "R" ) ipick_type = Rayleigh_R
              if( cmp(1:1) == "T" ) ipick_type = Love_T
           endif

          ! LQY: shouldn't chi values be changed accordingly?????
          ! No total chi value is calculated ...

          ! weights by phase types
          select case(ipick_type)
            case( P_SV_V )
              tr_adj_src(:) = tr_adj_src(:) * num_P_SV_V
              tr_chi = tr_chi * num_P_SV_V
            case( P_SV_R )
              tr_adj_src(:) = tr_adj_src(:) * num_P_SV_R
              tr_chi = tr_chi * num_P_SV_R
            case( SH_T )
              tr_adj_src(:) = tr_adj_src(:) * num_SH_T
              tr_chi = tr_chi * num_SH_T
            case( Rayleigh_V )
              tr_adj_src(:) = tr_adj_src(:) * num_Rayleigh_V
              tr_chi = tr_chi * num_Rayleigh_V
            case( Rayleigh_R )
              tr_adj_src(:) = tr_adj_src(:) * num_Rayleigh_R
              tr_chi = tr_chi * num_Rayleigh_R
            case( Love_T )
              tr_adj_src(:) = tr_adj_src(:) * num_Love_T
              tr_chi = tr_chi * num_Love_T
            case default
              stop 'error ipick_type unknown'
          end select
       endif

        ! combine adjoint sources from different measurement windows
       if (mod(imeas,2)==1) then
          adj_syn_all(:) = adj_syn_all(:) + tr_adj_src(:)   ! imeas = 1,3,5,7
          all_chi = all_chi + tr_chi
       else
          adj_syn_all(:) = adj_syn_all(:) + am_adj_src(:)   ! imeas = 2,4,6,8
          all_chi = all_chi + am_chi
       endif

        ! combine CC-reconstructed records
       if (imeas >= 7) then
          recon_cc_all(istart:istart+nlen-1) = recon_cc_all(istart:istart+nlen-1) + syn_dtw_mt(1:nlen)
       else
          recon_cc_all(istart:istart+nlen-1) = recon_cc_all(istart:istart+nlen-1) + syn_dtw_cc(1:nlen)
       endif

     endif ! COMPUTE_ADJOINT_SOURCE

      ! CHT: (re-)set to multitaper parameters, if originally specified
      if (is_mtm0 == 1) then
         imeas = imeas0
         is_mtm = is_mtm0
      endif

   enddo ! nmeas

    !----------------------------
    ! write out the adjoint source for the trace (STA.NI.CMP) by combining contribution from all wins

    if (COMPUTE_ADJOINT_SOURCE) then

    ! write out the CC-reconstructed data from synthetics
      !  if (OUTPUT_MEASUREMENT_FILES) &
      !       call dwsac1(trim(file_prefix2)//'.recon.sac',recon_cc_all,npts,t0,dt)

      ! OPTIONAL: A conservative choice is to filter the adjoint source,
      !   since higher frequencies could enter from the tapering operations.
      ! Note: time_window in mt_adj.f90 tapers the windows.

      ! note also:
      ! measurements are done on filtered synthetics F(s) and filtered data F(d), such that DeltaT
      ! is given for filtered data & synthetics.
      ! then kernels,
      ! i.e. for a traveltime measurement: DeltaT = 1/N * int  F(d/dt s) F(ds)
      ! should contain this filter as well.
      !
      ! when we construct the adjoint source here,it is initially a filtered version
      ! as well F(s_adj) since we use/depend on filtered synthetics F(s).
      ! however, for kernel simulations, we do run with a reconstructed forward wavefield,
      ! which is unfiltered (only filter there is by source half-time), but we want to convolve
      !  K = int F*(s_adj) F(s)
      ! using the same (bandpass) filter F() as used for filtereing data & synthetics in the meausurements
      ! We can write the kernel expression as K = int F*{F* (s_adj)}  s
      ! thus we should apply the filter F() twice on the adjoint source
      !
      ! why is this important? the filter, like bandpassing, is usually acausal, that is, it can
      ! introduce a slight phase-shift to the data. but, phase-shifts is what we are interested in
      ! and invert for. so, filtering might affect our inversions...

      ! we do use a bandpass filter here again on the adjoint source. this is slightly different
      ! to the transfer function filter in SAC used initially to filter data & synthetics.
      ! but this seems to be the best and fairly easy what we can do here...
      call bandpass(adj_syn_all,npts,dt,fstart0,fend0) ! sac butterworth filter

      ! cut and interpolate to match time-stepping for SEM
      ! NOTE: This can leave a non-zero value to start the record,
      !       which is NOT GOOD for the SEM simulation.
      call interpolate_syn(adj_syn_all,t0,dt,npts,-SPECFEM_T0,SPECFEM_DT,NSTEP)

      ! Taper the start of the adjoint source, since cutting the record
      ! may have left a non-zero value to start the record,
      ! which is not good for the SEM simulation.
      itmax = int(TSHORT/SPECFEM_DT)
      call taper_start(adj_syn_all,NSTEP,itmax)

      ! output the adjoint source (or ray density) as ASCII or SAC format
      !print *, 'writing adjoint source to file for the full seismogram'
      if( DO_RAY_DENSITY_SOURCE ) then
        call dwascii(trim(adj_file_prefix)//'.density.adj',adj_syn_all,NSTEP,-SPECFEM_T0,SPECFEM_DT)
      ! else
      !   call dwsac1(trim(OUTPUT_FILES)//'/../SEM/'//trim(adj_file_prefix)// &
      !               '.adj.sac'//'.'//trim(bandname),adj_syn_all,NSTEP,-SPECFEM_T0,SPECFEM_DT)
        ! LQY add the sum of chi values (total misfit), weights included if DO_WEIGHTING on
      endif

    endif

  enddo ! npairs

  !close(11)  ! read: MEASUREMENT.WINDOWS
  !close(12)  ! write: window_index
  !close(13)  ! write: window_chi

 end subroutine measure_adj_fwat


subroutine measure_adj()

  !  main program that calls the subroutines to make measurements within input time windows
  !  and then compute the corresponding adjoint sources

  ! input parameter:
  !  1. imeas = 1, normalized waveform difference. Adjoint source is constructed from the data
  !                      only, with the form -d(t)/ || d(t) || 2
  !  2. imeas = 2, waveform difference, s(t) - d(t).
  !  3. imeas = 3, cross-correlation traveltime difference for a (banana-doughtnut) sensitivity ker-
  !                       nel. The measurement between data and synthetics is not used in constructing the adjoint
  !                       source.
  !  4. imeas = 4, amplitude difference for a (banana-doughtnut) sensitivity kernel. The measure-
  !                       ment between data and synthetics is not used in constructing the adjoint source.
  !  5. imeas = 5, cross-correlation traveltime difference for an event kernel. The measurement
  !                       between data and synthetics is used in constructing the adjoint source.
  !  6. imeas = 6, amplitude difference for an event kernel. The measurement between data and
  !                        synthetics is used in constructing the adjoint source.
  !  7. imeas = 7, multitaper traveltime difference for an event kernel. The measurement between
  !                       data and synthetics is used in constructing the adjoint source. See multitaper_notes.pdf.
  !  8. imeas = 8, multitaper amplitude difference for an event kernel. The measurement between
  !                       data and synthetics is used in constructing the adjoint source. See multitaper_notes.pdf.

  use ma_variables
  use ma_constants
  use ascii_rw       ! dwascii()
  use ma_sub2        ! fft(), fftinv()
  use ma_sub         ! mt_measure(), mt_adj()
  use ma_weighting

  implicit none

  character(len=150) :: datafile,synfile,synfile_phydisp,file_prefix,file_prefix0,file_prefix2,measure_file_prefix,adj_file_prefix
  integer :: num_meas, j, ios, npt1, npt2,npt3, npts, nn
  double precision, dimension(NDIM_MA) :: data, syn, syn_phydisp, adj_syn_all, &
                        tr_adj_src, am_adj_src, recon_cc_all, syn_dtw_cc, syn_dtw_mt
  double precision :: t01, dt1, t02, dt2, t03, dt3, t0, dt, tstart, tend, tt, dtt, df
  double precision, dimension(NCHI) :: window_chi
  double precision :: fend0, fstart0, fend, fstart

  ! sac header information
  integer :: yr,jda,ho,mi
  double precision :: sec,dist,az,baz,slat,slon
  character(len=10) :: net,sta,chan_dat,chan,cmp,chan_syn
  double precision :: tshift, sigma_dt_cc, dlnA, sigma_dlnA_cc, sigma_dt, sigma_dlnA
  double precision :: all_chi, tr_chi, am_chi, cc_max, T_pmax_dat, T_pmax_syn
  !double precision :: tshift_f1f2, cc_max_f1f2
  double precision, dimension(NPT) :: dtau_w, dlnA_w, err_dt, err_dlnA, syn_dtw, data_dtw,syn_dtw_phydisp
  complex(kind=8), dimension(NPT) :: trans_mtm
  integer :: nlen, i_left, i_pmax_dat, i_pmax_syn, i_right, i_right0, istart, &
        ipair, npairs, nwin, itmax
  logical :: use_trace
  !double precision :: trbdndw, a
  !integer :: iord, passes
  integer :: ipick_type
  double precision :: T_surfacewaves

  !********* PROGRAM STARTS HERE *********************
  ! read in MEASUREMENT.PAR (see ma_sub.f90 and write_par_file.pl)
  ! most parameters are global (see ma_variables.f90)
  call read_par_file(fstart0,fend0,tt,dtt,nn,chan)

  ! uses weights to balance love and rayleigh measurements
  ! we do a normalization of P_SV, P_SH, Love, Rayleigh with the number of measurement picks
  if( DO_WEIGHTING ) call setup_weighting(chan)

  ! input file: MEASUREMENT.WINDOWS
  open(11,file='MEASUREMENT.WINDOWS',status='old',iostat=ios)
  if (ios /= 0) stop 'Error opening input file: MEASUREMENT WINDOWS'

  read(11,*,iostat=ios) npairs
  if (ios /= 0) stop 'Error reading number of pairs of data/syn'
  print *, 'reading in the data and synthetics...'

  ! output files
  open(12,file='window_index',status='unknown',iostat=ios)
  open(13,file='window_chi',status='unknown',iostat=ios)

  nwin = 0; all_chi=0.
  do ipair = 1, npairs

    data(:) = 0.0
    syn(:)  = 0.0

    adj_syn_all(:) = 0.0
    recon_cc_all(:) = 0.0

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
    if (max(npt1,npt2) > NDIM_MA) &
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

    if(RUN_BANDPASS) then
       call bandpass(data,npts,dt,fstart0,fend0)
       call bandpass(syn,npts,dt,fstart0,fend0)
       if (USE_PHYSICAL_DISPERSION) then
               call bandpass(syn_phydisp,npts,dt,fstart0,fend0)
       endif
    endif

    ! find out station/network/comp names,etc from synthetics
    call get_sacfile_header(trim(synfile),yr,jda,ho,mi,sec,net,sta, &
                          chan_dat,dist,az,baz,slat,slon)

    ! theoretical surface wave arrival time
    T_surfacewaves = dist / surface_vel

    ! synthetics always have the form BH_ or LH_, but the data may not (HH_, LH_, BL_, etc).
    cmp = chan_dat(3:3)
    chan_syn = trim(chan)//trim(cmp)

    ! example: OUT/PAS.CI.BHZ
    file_prefix0 = trim(sta)//'.'//trim(net)//'.'//trim(chan_syn)
    file_prefix2 = trim(OUT_DIR)//'/'//trim(file_prefix0)
    print *
    print *,  trim(file_prefix2), ' --- '

    ! note: MT measurement could revert to CC, but still keep the MT suffix
    write(adj_file_prefix,'(a,i2.2)') trim(file_prefix2)//'.iker', imeas0

    ! reads number of measurement windows
    read(11,*,iostat=ios) num_meas
    if (ios /= 0) stop 'Error reading num_meas'

    do j = 1, num_meas
      ! reads in start and end time of the measurement window
      read(11,*,iostat=ios) tstart, tend
      if (ios /= 0) stop 'Error reading tstart and tend'

      ! checks start and end times of window compared to trace lengths
      tstart = max(tstart,t0)
      tend = min(tend, t0+(npts-1)*dt)
      nlen = floor((tend-tstart)/dt) + 1  ! dummy, replaced later in mt_measure()

      ! write values to output file
      nwin = nwin + 1       ! overall window counter
      write(12,'(a3,a8,a5,a5,3i5,2f12.3)') net,sta,chan_syn,chan_dat,nwin,ipair,j,tstart,tend

      ! add taper type to file prefix: OUT/PAS.CI.BHZ.01.mtm
      write(file_prefix,'(a,i2.2)') trim(file_prefix2)//'.', j

      if (is_mtm == 1) then
        measure_file_prefix = trim(file_prefix) // '.mtm'  ! multitaper taper
      else if (is_mtm == 2) then
        measure_file_prefix = trim(file_prefix) // '.ctp'  ! cosine taper
      else
        measure_file_prefix = trim(file_prefix) // '.btp'  ! boxcar taper
      endif

      print *
      print *, ' Measurement window No.', j, ' ... '

      ! initialize the measurements
      window_chi(:) = 0.

      ! compute integrated waveform difference, normalized by duration of the record
      ! NOTE: (1) this is for the FULL record, not the windowed record
      !       (2) for comparison with waveform_chi, we include the 0.5 factor
      !       (3) we might want to include dt as an integration factor (also for waveform_chi),
      !           but the ratio (d-s)^2 / d^2 avoids the need for dt, nstep, or length of record
      window_chi(17) = 0.5 * sum( data**2 )
      window_chi(18) = 0.5 * sum( syn**2 )
      window_chi(19) = 0.5 * sum( (data-syn)**2 )
      window_chi(20) = npts*dt

      ! make measurements
      ! also compute reconstructed synthetics for CC (and MT, if specified) measurements
      call mt_measure(datafile,measure_file_prefix,data,syn,syn_phydisp,t0,dt,npts,tstart,tend,&
            istart,data_dtw,syn_dtw,syn_dtw_phydisp,nlen,tshift,sigma_dt_cc,dlnA,sigma_dlnA_cc,cc_max,syn_dtw_cc,&
            i_pmax_dat,i_pmax_syn,i_right0,trans_mtm,dtau_w,dlnA_w,sigma_dt,sigma_dlnA,syn_dtw_mt,err_dt,err_dlnA)
      i_right = i_right0
      i_left = 1  ! LQY: is it feasible that i_left is not 1? mt_adj() inherently assumes it.

      ! period of the max power of the synthetic record
      T_pmax_dat = (dt*NPT) / dble(i_pmax_dat)
      T_pmax_syn = (dt*NPT) / dble(i_pmax_syn)

      ! adjust frequency ranges for MT measurements
      ! fstart is constrained by NCYCLE_IN_WINDOW/tlen, fend constrained by i_right
      if (is_mtm == 1) then
         fstart = fstart0  ; fend = fend0
         call mt_measure_select(nlen,tshift,i_pmax_syn,dtau_w,err_dt, &
                              dt,i_left,i_right,fstart,fend,use_trace)
         print *, '     Tlong/Tshort (input) :', sngl(1/fstart0), sngl(1/fend0)
         print *, '     Tlong/Tshort (adjusted)  :', sngl(1/fstart), sngl(1/fend)
         print *, '     period of max data/syn power    :', sngl(T_pmax_dat), sngl(T_pmax_syn)

         ! if MT measurement window is rejected by mt_measure_select, then use a CC measurement
         if(.not. use_trace) then
            !stop 'Check why this MT measurement was rejected'
            print *, '   reverting from MT measurement to CC measurement...'
            imeas = imeas0 - 2
            is_mtm = 3  ! LQY: WHY not is_mtm = 2?
            call mt_measure(datafile,measure_file_prefix,data,syn,syn_phydisp,t0,dt,npts,tstart,tend,&
                  istart,data_dtw,syn_dtw,syn_dtw_phydisp,nlen,&
                  tshift,sigma_dt_cc,dlnA,sigma_dlnA_cc,cc_max,syn_dtw_cc,&
                  i_pmax_dat,i_pmax_syn,i_right,trans_mtm,dtau_w,dlnA_w,sigma_dt,sigma_dlnA,syn_dtw_mt)
         else
            print *, '     using this MTM. '
         endif
      endif

      ! check that the CC measurements are within the specified input range
      if (imeas >= 5) call cc_measure_select(tshift,dlnA,cc_max)

      ! write frequency limits to file
      if (OUTPUT_MEASUREMENT_FILES) then
        df = 1./(dt*NPT)
        open(71,file=trim(measure_file_prefix)//'.freq_limits')
        write(71,'(6f18.8)') fstart0, fend0, df, i_right0*df, fstart, fend
        close(71)
      endif

      ! compute adjoint sources and misfit function values and also the CC-reconstructed records
      if (COMPUTE_ADJOINT_SOURCE) then
         print *, '   Generating adjoint source and chi value for imeas = ', imeas

        ! banana-doughnut kernel (needs only synthetic trace)
        ! LQY: what is this section intended to do?
        ! reset imeas == 3 for adjoint sources without time shift and uncertainty scaling
        ! (pure cross-correlation adjoint source for banana-doughnuts)
        if(imeas == 5 .and. trim(datafile) == trim(synfile) ) then
           print *,'cross-correlation measurement:'
           print *,'  only synthetic file: ',trim(synfile)
           print *,'    without traveltime difference/uncertainty'
           print *
           imeas = 3
        endif

        tr_chi = 0.0 ; am_chi = 0.0    ! must be initialized
        call mt_adj(istart,data_dtw,syn_dtw,syn_dtw_phydisp,nlen,dt,tshift,dlnA,sigma_dt_cc,sigma_dlnA_cc,&
             dtau_w,dlnA_w,err_dt,err_dlnA,sigma_dt,sigma_dlnA,i_left,i_right,&
             window_chi,tr_adj_src,tr_chi,am_adj_src,am_chi)

        ! KEY: write misfit function values to file (two for each window)
        ! Here are the 20 columns of the vector window_chi
        !  1: MT-TT chi,    2: MT-dlnA chi,    3: XC-TT chi,    4: XC-dlnA chi
        !  5: MT-TT meas,   6: MT-dlnA meas,   7: XC-TT meas,   8: XC-dlnA meas
        !  9: MT-TT error, 10: MT-dlnA error, 11: XC-TT error, 12: XC-dlnA error
        ! WINDOW     : 13: data power, 14: syn power, 15: (data-syn) power, 16: window duration
        ! FULL RECORD: 17: data power, 18: syn power, 19: (data-syn) power, 20: record duration
        ! Example of a reduced file: awk '{print $2,$3,$4,$5,$6,$31,$32}' window_chi > window_chi_sub
        write(13,'(a14,a8,a3,a5,i4,i4,2e14.6,20e14.6,2e14.6,2f14.6)') &
           file_prefix0,sta,net,chan_syn,j,imeas,&
           tstart,tend,window_chi(:),tr_chi,am_chi,T_pmax_dat,T_pmax_syn
        print *, '     tr_chi = ', sngl(tr_chi), '  am_chi = ', sngl(am_chi)

        ! uses weighting to balance love / rayleigh measurements
        if( DO_WEIGHTING ) then
           ipick_type = 0
           if( tend <= T_surfacewaves ) then
              ! body wave picks
              if( cmp(1:1) == "Z" ) ipick_type = P_SV_V
              if( cmp(1:1) == "R" ) ipick_type = P_SV_R
              if( cmp(1:1) == "T" ) ipick_type = SH_T
           else
              ! surface wave picks
              if( cmp(1:1) == "Z" ) ipick_type = Rayleigh_V
              if( cmp(1:1) == "R" ) ipick_type = Rayleigh_R
              if( cmp(1:1) == "T" ) ipick_type = Love_T
           endif

          ! LQY: shouldn't chi values be changed accordingly?????
          ! No total chi value is calculated ...

          ! weights by phase types
          select case(ipick_type)
            case( P_SV_V )
              tr_adj_src(:) = tr_adj_src(:) * num_P_SV_V
              tr_chi = tr_chi * num_P_SV_V
            case( P_SV_R )
              tr_adj_src(:) = tr_adj_src(:) * num_P_SV_R
              tr_chi = tr_chi * num_P_SV_R
            case( SH_T )
              tr_adj_src(:) = tr_adj_src(:) * num_SH_T
              tr_chi = tr_chi * num_SH_T
            case( Rayleigh_V )
              tr_adj_src(:) = tr_adj_src(:) * num_Rayleigh_V
              tr_chi = tr_chi * num_Rayleigh_V
            case( Rayleigh_R )
              tr_adj_src(:) = tr_adj_src(:) * num_Rayleigh_R
              tr_chi = tr_chi * num_Rayleigh_R
            case( Love_T )
              tr_adj_src(:) = tr_adj_src(:) * num_Love_T
              tr_chi = tr_chi * num_Love_T
            case default
              stop 'error ipick_type unknown'
          end select
       endif

        ! combine adjoint sources from different measurement windows
       if (mod(imeas,2)==1) then
          adj_syn_all(:) = adj_syn_all(:) + tr_adj_src(:)   ! imeas = 1,3,5,7
          all_chi = all_chi + tr_chi
       else
          adj_syn_all(:) = adj_syn_all(:) + am_adj_src(:)   ! imeas = 2,4,6,8
          all_chi = all_chi + am_chi
       endif

        ! combine CC-reconstructed records
       if (imeas >= 7) then
          recon_cc_all(istart:istart+nlen-1) = recon_cc_all(istart:istart+nlen-1) + syn_dtw_mt(1:nlen)
       else
          recon_cc_all(istart:istart+nlen-1) = recon_cc_all(istart:istart+nlen-1) + syn_dtw_cc(1:nlen)
       endif

     endif ! COMPUTE_ADJOINT_SOURCE

      ! CHT: (re-)set to multitaper parameters, if originally specified
      if (is_mtm0 == 1) then
         imeas = imeas0
         is_mtm = is_mtm0
      endif

   enddo ! nmeas

    !----------------------------
    ! write out the adjoint source for the trace (STA.NI.CMP) by combining contribution from all wins

    if (COMPUTE_ADJOINT_SOURCE) then

    ! write out the CC-reconstructed data from synthetics
       if (OUTPUT_MEASUREMENT_FILES) &
            call dwsac1(trim(file_prefix2)//'.recon.sac',recon_cc_all,npts,t0,dt)

      ! OPTIONAL: A conservative choice is to filter the adjoint source,
      !   since higher frequencies could enter from the tapering operations.
      ! Note: time_window in mt_adj.f90 tapers the windows.

      ! note also:
      ! measurements are done on filtered synthetics F(s) and filtered data F(d), such that DeltaT
      ! is given for filtered data & synthetics.
      ! then kernels,
      ! i.e. for a traveltime measurement: DeltaT = 1/N * int  F(d/dt s) F(ds)
      ! should contain this filter as well.
      !
      ! when we construct the adjoint source here,it is initially a filtered version
      ! as well F(s_adj) since we use/depend on filtered synthetics F(s).
      ! however, for kernel simulations, we do run with a reconstructed forward wavefield,
      ! which is unfiltered (only filter there is by source half-time), but we want to convolve
      !  K = int F*(s_adj) F(s)
      ! using the same (bandpass) filter F() as used for filtereing data & synthetics in the meausurements
      ! We can write the kernel expression as K = int F*{F* (s_adj)}  s
      ! thus we should apply the filter F() twice on the adjoint source
      !
      ! why is this important? the filter, like bandpassing, is usually acausal, that is, it can
      ! introduce a slight phase-shift to the data. but, phase-shifts is what we are interested in
      ! and invert for. so, filtering might affect our inversions...

      ! we do use a bandpass filter here again on the adjoint source. this is slightly different
      ! to the transfer function filter in SAC used initially to filter data & synthetics.
      ! but this seems to be the best and fairly easy what we can do here...
      call bandpass(adj_syn_all,npts,dt,fstart0,fend0) ! sac butterworth filter

      ! cut and interpolate to match time-stepping for SEM
      ! NOTE: This can leave a non-zero value to start the record,
      !       which is NOT GOOD for the SEM simulation.
      call interpolate_syn(adj_syn_all,t0,dt,npts,tt,dtt,nn)

      ! Taper the start of the adjoint source, since cutting the record
      ! may have left a non-zero value to start the record,
      ! which is not good for the SEM simulation.
      itmax = int(TSHORT/dtt)
      call taper_start(adj_syn_all,nn,itmax)

      ! output the adjoint source (or ray density) as ASCII or SAC format
      print *, 'writing adjoint source to file for the full seismogram'
      if( DO_RAY_DENSITY_SOURCE ) then
        call dwascii(trim(adj_file_prefix)//'.density.adj',adj_syn_all,nn,tt,dtt)
      else
        call dwascii(trim(adj_file_prefix)//'.adj',adj_syn_all,nn,tt,dtt)
        ! LQY add the sum of chi values (total misfit), weights included if DO_WEIGHTING on
        open(14,file='window_chi_sum',status='unknown')
        write(14,*) all_chi
        close(14)
      endif

    endif

  enddo ! npairs

  close(11)  ! read: MEASUREMENT.WINDOWS
  close(12)  ! write: window_index
  close(13)  ! write: window_chi

 end subroutine measure_adj

  ! subroutine measure_adj_rf_data(data,tstart,tend,t0,tp,dt,npts,f0,tshift,net,sta,chan_dat,&
  !                               bandname,adj_r_tw, adj_z_tw)
  !   use decon_mod
  !   use interpolation_mod, only : PI

  !   implicit none
  !   integer, intent(in)                                :: npts
  !   character(len=MAX_STRING_LEN), intent(in)          :: bandname
  !   double precision, dimension(npts), intent(in)      :: data
  !   character(len=10), intent(in)                      :: net, sta,chan_dat
  !   integer                                            :: nstart, nend, nb, i, n
  !   double precision, intent(in)                       :: tstart, tend,t0, dt, tp
  !   real, intent(in)                                   :: f0, tshift
  !   double precision, dimension(npts)                  :: adj_r, adj_z, syn
  !   double precision, dimension(npts), intent(inout)   :: adj_r_tw, adj_z_tw
  !   character(len=MAX_STRING_LEN)                      :: adj_file_prefix
  !   double precision                                   :: fac

  !   adj_r_tw(:) = 0.
  !   adj_z_tw(:) = 0.
  !   syn(:) = 0.

  !   nstart = floor((-tstart + tshift)/dt+1)
  !   nend = floor((tend + tshift)/dt+1)
  !   nb = floor((tp - tshift - t0)/dt+1)
  !   adj_r = -data(1:npts)
  !   syn(int(tshift/dt+1)) = 1.
  !   adj_z = 0.
  !   call deconit(syn, syn, npts, real(dt), tshift, f0, 10, 0.001, 0, adj_z)
  !   adj_z = -adj_z

  !   n = 1
  !   ! Cosine taper
  !   do i = nstart, nend
  !     fac = 1. - cos(PI*(n-1)/(nend-nstart))**10
  !     adj_r_tw(nb+i) = adj_r(i) * fac
  !     adj_z_tw(nb+i) = adj_z(i) * fac
  !     n = n+1
  !   enddo

  !   ! write windowed adjoint source
  !   adj_file_prefix = trim(net)//'.'//trim(sta)//'.'//trim(chan_dat)//'R'
  !   call dwsac1(trim(OUTPUT_FILES)//'/../SEM/'//trim(adj_file_prefix)//'.adj.sac'//'.'//trim(bandname),adj_r_tw,npts,t0,dt)
  !   adj_file_prefix = trim(net)//'.'//trim(sta)//'.'//trim(chan_dat)//'Z'
  !   call dwsac1(trim(OUTPUT_FILES)//'/../SEM/'//trim(adj_file_prefix)//'.adj.sac'//'.'//trim(bandname),adj_z_tw,npts,t0,dt)

  ! end subroutine measure_adj_rf_data


  subroutine meas_adj_conv_diff(datr, datz, synr, synz, tstart, tend, tp, npts,&
                                tshift, f0, maxit, minderr, freq_min, freq_max, &
                                window_chi, adj_r_tw, adj_z_tw, sta)
    !============= MJ: measure adjoint sourece for ||Dr*Sz - Dz*Sr|| =====================
    use utils, only: zeros_dp
    use signal, only: myconvolution_dp, bandpass_dp
    use decon_mod, only: deconit

    implicit none
    integer, intent(in)                                :: npts, maxit
    double precision, intent(in)                       :: tshift, f0, minderr, tp, freq_min, freq_max
    double precision, dimension(NCHI), intent(inout)   :: window_chi
    double precision, dimension(npts)                  :: r_rev, z_rev, synr_shift, synz_shift
    double precision, dimension(npts), intent(in)      :: datr, datz, synr, synz
    double precision, intent(in)                       :: tstart, tend
    character(len=MAX_STRING_LEN)                      :: adj_file_prefix, sta
    double precision, dimension(:), allocatable        :: conv_den, conv_num, conv_full, conv_diff,&
                                                          conv1, conv2, tmp, adj_r, adj_z, &
                                                          data_tw, synt_tw, obj_cc, obj_cc_tw, zrf
    double precision, dimension(:), allocatable, intent(out) :: adj_r_tw, adj_z_tw
    double precision                                   :: fac, e
    integer                                            :: n, nc, nb, nstart, nend,i, nstop

    ! n = int((tend - tstart)/SPECFEM_DT)
    nc = floor((2*(tp + SPECFEM_T0)-tshift)/ SPECFEM_DT + 1)
    nb = floor((tp + tstart + SPECFEM_T0)/SPECFEM_DT+1)

    call deconit(synz, synz, real(SPECFEM_DT), 10., real(f0), 10, 0.001, 0, zrf)
    ! call myconvolution_dp(synr, datz, conv1, 1)
    ! call myconvolution_dp(synz, datr, conv2, 1)
    call deconit(synr, synz, real(SPECFEM_DT), real(tshift), real(f0), maxit, real(minderr), 0, conv1)
    call deconit(datr, datz, real(SPECFEM_DT), real(tshift), real(f0), maxit, real(minderr), 0, conv2)
    ! conv_diff = (conv1 - conv2)/maxval(zrf)
    obj_cc = (conv1 - conv2)/maxval(zrf)
    ! call myconvolution_dp(datz, synz, tmp, 1)
    ! call deconit(conv_diff, tmp, real(SPECFEM_DT), real(tshift), real(f0), maxit, real(minderr), 1, obj_cc)
    ! conv_full = zeros_dp(npts)
    ! conv_full = tmp(nc:npts+nc-1)
    e = npts * SPECFEM_DT - tshift

    ! call reverse(synr, npts, r_rev)
    ! call reverse(synz, npts, z_rev)
    synr_shift = 0.
    synz_shift = 0.
    synr_shift(1:npts-nb+1) = synr(nb:npts)
    synz_shift(1:npts-nb+1) = synz(nb:npts)
    call reverse(synr_shift, npts, r_rev)
    call reverse(synz_shift, npts, z_rev)
    ! radial adjoint source
    ! call myconvolution_dp(conv_full, z_rev, conv_den, 1)
    call deconit(obj_cc(1:npts), z_rev, real(SPECFEM_DT), real(e), real(f0), maxit, real(minderr), 1, adj_r)

    ! vertical adjoint source
    call myconvolution_dp(-obj_cc(1:npts), r_rev, conv_num, 1)
    call myconvolution_dp(z_rev, z_rev, conv_den, 1)
    ! call myconvolution_dp(tmp, conv_full, conv_den, 1)
    call deconit(conv_num, conv_den, real(SPECFEM_DT), real(e), real(f0), maxit, real(minderr), 1, adj_z)

    ! adj_file_prefix = 'MX.'//trim(sta)//'.BXR'
    ! call dwsac1(trim(OUTPUT_FILES)//'/../SEM/'//trim(adj_file_prefix)//'.objcc.sac',obj_cc(1:npts),npts,-tshift,SPECFEM_DT)
    ! adj_file_prefix = 'MX.'//trim(sta)//'.BXZ'
    ! call dwsac1(trim(OUTPUT_FILES)//'/../SEM/'//trim(adj_file_prefix)//'.adj.sac',adj_z(1:npts),npts,-tshift,SPECFEM_DT)
    ! adj_file_prefix = 'MX.'//trim(sta)//'.BXR'
    ! call dwsac1(trim(OUTPUT_FILES)//'/../SEM/'//trim(adj_file_prefix)//'.adj.sac',adj_r(1:npts),npts,-tshift,SPECFEM_DT)

    nstart = floor((tstart + tshift)/SPECFEM_DT+1)
    nend = floor((tend + tshift)/SPECFEM_DT+1)
    ! nstop = min(nend, npts)

    adj_r_tw = zeros_dp(npts)
    adj_z_tw = zeros_dp(npts)
    data_tw = zeros_dp(npts)
    synt_tw = zeros_dp(npts)
    obj_cc_tw = zeros_dp(npts)
    n = 1
    ! Cosine taper
    do i = nstart, nend
      fac = 1. - cos(PI*(n-1)/(nend-nstart))**10
      data_tw(i) = conv1(i) * fac
      synt_tw(i) = conv2(i) * fac
      obj_cc_tw(i) = obj_cc(i) * fac
      adj_r_tw(nb+i) = adj_r(i) * fac
      adj_z_tw(nb+i) = adj_z(i) * fac
      n = n+1
    enddo
    ! call bandpass_dp(obj_cc_tw,npts,SPECFEM_DT,real(freq_min),real(freq_max),IORD)
    ! call bandpass_dp(adj_r_tw,npts,SPECFEM_DT,real(freq_min),real(freq_max),IORD)
    ! call bandpass_dp(adj_z_tw,npts,SPECFEM_DT,real(freq_min),real(freq_max),IORD)
             
    ! write windowed adjoint source
    window_chi = 0.
    window_chi(13) = 0.5 * sum( data_tw**2 )
    window_chi(14) = 0.5 * sum( synt_tw**2 )
    window_chi(15) = 0.5 * sum( obj_cc_tw**2 )
    window_chi(16) = (nend - nstart)*SPECFEM_DT
    window_chi(17) = 0.5 * sum( conv1**2 )
    window_chi(18) = 0.5 * sum( conv2**2 )
    window_chi(19) = 0.5 * sum( (obj_cc)**2 )
    window_chi(20) = npts*SPECFEM_DT

  end subroutine meas_adj_conv_diff

  subroutine measure_adj_rf(data,syn,synr,synz,tstart,tend,t0,tp,dt,npts,f0,tshift,maxit,minderr,&
                            window_chi, adj_r_tw,adj_z_tw, net, sta)
    use decon_mod, only : deconit
    use signal, only : myconvolution
    use fwat_constants, only : PI
    use utils, only: zeros_dp, zeros

    implicit none
    integer, intent(in)                                :: npts, maxit
    double precision, dimension(npts)                  :: r_rev, z_rev, diff_data, &
                                                          synr_shift, synz_shift, data_tw, &
                                                          synt_tw, data_norm, syn_norm
    double precision, dimension(npts), intent(in)      :: data, syn, synr, synz
    double precision, intent(in)                       :: tstart, tend,t0, dt, tp
    real, intent(in)                                   :: minderr, f0, tshift
    double precision, dimension(:), allocatable, intent(out)      :: adj_r_tw, adj_z_tw
    double precision, dimension(NCHI), intent(inout)   :: window_chi
    integer                                            :: nn,nb,nstart,nend,n,i
    real, dimension(:), allocatable                    :: tmp_n, tmp_d
    double precision, dimension(:), allocatable        :: adj_z, zrf, adj_r
    character(len=MAX_STRING_LEN)                      :: adj_file_prefix
    character(len=*), intent(in) :: net, sta
    double precision                                   :: fac
    real                                               :: e

    call deconit(synz, synz, real(dt), 10., f0, 10, 0.001, 0, zrf)
    data_norm = data/maxval(zrf)
    syn_norm = syn/maxval(zrf)
    synr_shift = 0.
    synz_shift = 0.
    r_rev = 0.
    z_rev = 0.
    diff_data = 0.
    data_tw = 0.
    synt_tw = 0.
    adj_r_tw = zeros_dp(npts)
    adj_z_tw = zeros_dp(npts)
    ! Align syn waveform with P
    nb = floor((tp - tshift - t0)/dt+1)
    synr_shift(1:npts-nb+1) = synr(nb:npts)
    synz_shift(1:npts-nb+1) = synz(nb:npts)

    ! Calculate RF adjoint source in time domain
    e = real(npts * dt - tshift)
    diff_data = syn_norm-data_norm
    nn = npts*2-1
    ! tmp_n = zeros(nn)
    ! tmp_d = zeros(nn)
    ! adj_z = zeros_dp(nn)
    call reverse(synr_shift, npts, r_rev)
    call reverse(synz_shift, npts, z_rev)
    call deconit(diff_data, z_rev, real(dt), e, f0, maxit, minderr, 1, adj_r)
    call myconvolution(-real(diff_data),real(r_rev),tmp_n,1)
    call myconvolution(real(z_rev),real(z_rev),tmp_d,1)
    call deconit(dble(tmp_n), dble(tmp_d), real(dt), e, f0, maxit, minderr, 1, adj_z)
    nstart = floor((-tstart + tshift)/dt+1)
    nend = floor((tend + tshift)/dt+1)

    n = 1
    ! Cosine taper
    do i = nstart, nend
      fac = 1. - cos(PI*(n-1)/(nend-nstart))**10
      data_tw(i) = data_norm(i) * fac
      synt_tw(i) = syn_norm(i) * fac
      adj_r_tw(nb+i) = adj_r(i) * fac
      adj_z_tw(nb+i) = adj_z(i) * fac
      n = n+1
    enddo

    ! write windowed adjoint source
    ! adj_file_prefix = trim(net)//'.'//trim(sta)//'.BXR'
    ! call dwsac1(trim(OUTPUT_FILES)//'/../SEM/'//trim(adj_file_prefix)//'.adj.sac.F1.5',adj_r,npts,t0,dt)
    ! adj_file_prefix = trim(net)//'.'//trim(sta)//'.BXZ'
    ! call dwsac1(trim(OUTPUT_FILES)//'/../SEM/'//trim(adj_file_prefix)//'.adj.sac.F1.5',adj_r,npts,t0,dt)
    ! window_chi(:) = 0.

    ! compute integrated waveform difference, normalized by duration of the record
    ! NOTE: (1) this is for the FULL record, not the windowed record
    !       (2) for comparison with waveform_chi, we include the 0.5 factor
    !       (3) we might want to include dt as an integration factor (also for waveform_chi),
    !           but the ratio (d-s)^2 / d^2 avoids the need for dt, nstep, or length of record
    window_chi(13) = 0.5 * sum( data_tw**2 )
    window_chi(14) = 0.5 * sum( synt_tw**2 )
    window_chi(15) = 0.5 * sum( (synt_tw-data_tw)**2 )
    window_chi(16) = (nend - nstart +1) *dt
    window_chi(17) = 0.5 * sum( data_norm**2 )
    window_chi(18) = 0.5 * sum( syn_norm**2 )
    window_chi(19) = 0.5 * sum( (syn_norm-data_norm)**2 )
    window_chi(20) = npts*dt
  end subroutine measure_adj_rf

  subroutine rotate_adj_src

    implicit none

    character(len=150) :: bazch, zfile,tfile,rfile,efile,nfile

    double precision :: baz
    integer :: nptst, nptsr, npts
    logical :: r_exist, t_exist, z_exist
    double precision :: t0t,dtt,t0r,dtr, t0,dt, costh, sinth
    !integer, parameter :: NDIM = 40000  ! check ma_constants.f90
    double precision, dimension(NDIM_MA) :: zdata, rdata, tdata, edata, ndata

    call getarg(1,bazch)
    call getarg(2,zfile)
    call getarg(3,tfile)
    call getarg(4,rfile)
    call getarg(5,efile)
    call getarg(6,nfile)

    if (trim(bazch)=='' .or. trim(zfile) =='' .or. trim(tfile)=='' .or. &
              trim(rfile) == '' .or. trim(efile) =='' .or. trim(nfile) =='') then
      stop 'rotate_adj_src baz(radian!) zfile tfile rfile efile nfile'
    endif

    read(bazch,*) baz

    inquire(file=trim(tfile),exist=t_exist)
    inquire(file=trim(rfile),exist=r_exist)
    inquire(file=trim(zfile),exist=z_exist)

    ! initialization
    rdata = 0; tdata = 0

    ! at least one file (T,R,Z) should be present
    if (.not. t_exist .and. .not. r_exist) then
      if (.not. z_exist) stop 'At least one file should exist: zfile, tfile, rfile'
    ! need to read Z comp adjoint source for [to,dt,npts]
      write(*,*) "read Z file"
      call drascii(zfile,zdata,npts,t0,dt)
    endif

    ! read in T file
    if (t_exist) then
      call drascii(tfile,tdata,nptst,t0t,dtt)
      write(*,*) "read T file"
    endif
    ! read in R file
    if (r_exist) then
      call drascii(rfile,rdata,nptsr,t0r,dtr)
      write(*,*) "read R file"
    endif

    ! check consistency of t0,dt,npts
    if (t_exist .and. r_exist) then
      if (abs(t0t-t0r)>1.0e-2 .or. abs(dtt-dtr)>1.0e-2 .or. nptst /= nptsr) &
                stop 'check t0 and npts'
    endif
    if (t_exist) then
      t0 =t0t; dt = dtt; npts = nptst
    else if (r_exist) then
      t0 =t0r; dt = dtr; npts = nptsr
    endif
    if (NDIM_MA < npts) stop 'Increase NDIM_MA'

    ! rotate T,R to E,N based on baz (note in radian!)
    costh = cos(baz)
    sinth = sin(baz)
    edata(1:npts) = -costh * tdata(1:npts) - sinth * rdata(1:npts)
    ndata(1:npts) =  sinth * tdata(1:npts) - costh * rdata(1:npts)

    ! write E,N files
    call dwascii(efile,edata,npts,t0,dt)
    call dwascii(nfile,ndata,npts,t0,dt)

    ! write Z file if did not exist
    if (.not. z_exist) then
      tdata = 0.
      call dwascii(zfile,zdata,npts,t0,dt)
    endif

  end subroutine rotate_adj_src

  subroutine reverse(src, nt, rev)
      integer :: nt, i
      double precision, dimension(nt) :: src, rev
      rev(:) = 0.
      do i = 1, nt
          rev(i) = src(nt-i+1)
      enddo
  end subroutine
end module measure_adj_mod
