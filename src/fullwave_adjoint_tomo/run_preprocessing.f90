!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
!
subroutine run_preprocessing(model,evtset,ievt,simu_type,icmt)

  use specfem_par, only: CUSTOM_REAL, MAX_STRING_LEN, OUTPUT_FILES, IIN, network_name, station_name,nrec, myrank, t0
  use shared_input_parameters, only: NSTEP, DT, NPROC

  use specfem_par_elastic, only: ELASTIC_SIMULATION
  use specfem_par_acoustic, only: ispec_is_acoustic, ACOUSTIC_SIMULATION
  use fullwave_adjoint_tomo_par
  use fwat_input
  use fwat_utils
  use FKTimes_mod
  use my_mpi
  use measure_adj_mod, MAX_NDIM => NDIM ! module from preproc_measure_adj
  use telestf_mod      ! module form teleseis_stf
  use collect_data, only: collect_seismograms_d, collect_stf_deconvolution, collect_chi
  use preproc_subs, only: cal_fktimes, get_rf_times, pre_proc_tele_elastic,read_local_stf,&
                          average_amp_scale, pre_proc_rf_elastic, pre_proc_tele_cd_elastic,&
                          get_zdat_pca
  use preproc_measure_adj_subs, only: meas_adj_tele, meas_adj_rf, meas_adj_tele_cd,&
                                      meas_adj_noise, meas_adj_leq, sum_adj_source

  
  implicit none

  character(len=MAX_STRING_LEN)                            :: model 
  character(len=MAX_STRING_LEN)                            :: evtset 
  character(len=MAX_STRING_LEN)                            :: simu_type
  integer                                                  :: irec, irec_local, ispec
  integer                                                  :: ievt,icomp,iwin,icmt,ier
  double precision                                         :: fstart0,fend0
  character(len=MAX_STRING_LEN)                            :: bandname !,adjfile
  integer                                                  :: iband, nflt
  double precision, dimension(nrec)                        :: baz_all
  double precision, dimension(NCHI)                        :: window_chi
  double precision                                         :: tr_chi, am_chi
  ! distribute stations evenly to local proc for measure_adj
  integer                                                  :: my_nrec_local
  integer                                                  :: nrec_local_max
  ! for collecting SEM synthetics to one large array
  real(kind=CUSTOM_REAL), dimension(nrec,NSTEP,3)          :: glob_sem_disp
  ! Kai added for amplitude scaling of data gather 
  real                                                     :: avgamp
  ! for calculating STF for TeleFWI
  real(kind=4), dimension(nrec,NSTEP,NRCOMP)               :: glob_dat_tw,glob_syn_tw,glob_syn
  real(kind=4), dimension(nrec,NSTEP)                      :: glob_ff
  character(len=256), dimension(nrec)                      :: glob_stnm
  real(kind=4), dimension(nrec)                            :: ttp,tb,te
  real(kind=4) , dimension(NSTEP,NRCOMP)                   :: stf_array
  ! for writing window_chi of noise/tele by one proc
  character(len=256),dimension(:), allocatable             :: glob_file_prefix0
  character(len=256),dimension(:), allocatable             :: glob_net 
  character(len=256),dimension(:), allocatable             :: glob_sta  
  character(len=256),dimension(:,:), allocatable           :: glob_chan_dat   
  integer           ,dimension(:,:), allocatable           :: glob_num_win   
  double precision, dimension(:,:,:), allocatable          :: glob_tstart
  double precision, dimension(:,:,:), allocatable          :: glob_tend
  double precision, dimension(:,:,:,:), allocatable        :: glob_window_chi
  double precision, dimension(:,:,:), allocatable          :: glob_tr_chi
  double precision, dimension(:,:,:), allocatable          :: glob_am_chi
  double precision, dimension(:,:,:), allocatable          :: glob_T_pmax_dat
  double precision, dimension(:,:,:), allocatable          :: glob_T_pmax_syn
  integer         , dimension(:,:,:), allocatable          :: glob_imeas
  integer                                                  :: chi_fileid, win_fileid
  ! for writing STATIONS_***_ADJOINT
  character(len=MAX_STRING_LEN)                            :: dummystring, sta_name,net_name
  double precision                                         :: stlat,stlon,stele,stbur
  ! for writing adjoint sources
  integer                                                  :: num_adj
  double precision, dimension(3,MAX_NDIM)                  :: adj_syn_all_sum
  real(kind=CUSTOM_REAL)                                   :: total_misfit,total_misfit_reduced

  !
  if (myrank==0) then
    ! write(OUT_FWAT_LOG,*) 'This is run_preprocessing ...' 
    call write_timestamp_log(OUT_FWAT_LOG, 'This is run_preprocessing ...')
    flush(OUT_FWAT_LOG)
  endif
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !-------- collect synthetic from different PROCs to one large array   ------------------- 
  call collect_seismograms_d(glob_sem_disp)
  call bcast_all_cr(glob_sem_disp,nrec*NSTEP*3)
  if ( leq_par%CMT3D_INV ) then
    if (myrank==0) write(*,*) 'running forward simulatoin for ', icmt,'-th CMT_PAR:',trim(leq_par%CMT_PAR(icmt))
  endif
  if (myrank==0 .and. icmt==0) write(*,*) 'Now distridute ',nrec,' stations to ',NPROC,'procs'
  nrec_local_max=ceiling(real(nrec)/real(NPROC))
  my_nrec_local=0
  do irec_local = 1, nrec_local_max
    irec = myrank*nrec_local_max+irec_local
    if(irec<=nrec) my_nrec_local=my_nrec_local+1
  enddo
  if (icmt==0) write(*,*) 'Rank ',myrank,' has ',my_nrec_local,' stations for measure_adj '
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (simu_type =='leq') then 
    leq_par%FLEX_NWIN=10
  else
    leq_par%FLEX_NWIN=1
  endif
  allocate(glob_file_prefix0(nrec))
  allocate(glob_net(nrec))
  allocate(glob_sta(nrec)) 
  allocate(glob_chan_dat(nrec,NRCOMP))   
  allocate(glob_num_win(nrec,NRCOMP))  
  allocate(glob_tstart(nrec,NRCOMP,leq_par%FLEX_NWIN))
  allocate(glob_tend(nrec,NRCOMP,leq_par%FLEX_NWIN))  
  allocate(glob_window_chi(nrec,NRCOMP,NCHI,leq_par%FLEX_NWIN))
  allocate(glob_tr_chi(nrec,NRCOMP,leq_par%FLEX_NWIN))
  allocate(glob_am_chi(nrec,NRCOMP,leq_par%FLEX_NWIN))
  allocate(glob_T_pmax_dat(nrec,NRCOMP,leq_par%FLEX_NWIN))
  allocate(glob_T_pmax_syn(nrec,NRCOMP,leq_par%FLEX_NWIN))
  allocate(glob_imeas(nrec,NRCOMP,leq_par%FLEX_NWIN)) 

  if (simu_type == 'rf') then
    nflt = rf_par%NGAUSS
  else
    nflt = NUM_FILTER
  endif
  do iband=1,nflt
    if (simu_type /= 'rf') then
      fstart0=1./LONG_P(iband)
      fend0=1./SHORT_P(iband)
      call get_band_name(SHORT_P(iband),LONG_P(iband),bandname)
      ! write(bandname,'(a1,i3.3,a2,i3.3)') 'T',int(SHORT_P(iband)),'_T',int(LONG_P(iband))
      !********* 
      if (myrank==0) chi_fileid=30+iband; win_fileid=100+iband
      if (myrank==0.and.ievt==1) then
        open(chi_fileid,file='misfits/'//trim(model)//'.'//trim(evtset)//'_'//trim(bandname)&
                              //'_window_chi',status='unknown',iostat=ier)
        if (leq_par%CMT3D_INV) then
          open(win_fileid,file='misfits/'//trim(model)//'.'//trim(evtset)//'_'//trim(bandname)&
                                  //'.win',status='unknown',iostat=ier)
        endif
      endif
    else
      fstart0=1./LONG_P(1)
      fend0=1./SHORT_P(1)
      write(bandname,'(a1,f3.1)') 'F',rf_par%f0(iband)
      if (myrank==0) chi_fileid=30+iband
      if (myrank==0.and.ievt==1) then
        open(chi_fileid,file='misfits/'//trim(model)//'.'//trim(evtset)//'_'//trim(bandname)&
              //'_window_chi',status='unknown',iostat=ier)
      endif
    endif ! simu_type /= 'rf'
    glob_num_win=0
    glob_tstart=0.
    glob_tend=0.
    glob_window_chi=0.
    glob_tr_chi=0.
    glob_am_chi=0.
    glob_T_pmax_dat=0.
    glob_T_pmax_syn=0.
    glob_imeas=0
    adj_syn_all_sum=0.
    baz_all=0.
    glob_dat_tw=0.
    glob_syn_tw=0.
    glob_syn=0.
    glob_ff=0.
    total_misfit=0.
    !*********

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !========================= Calculate RFs for RFFWI =================================
    ! Mijian added at Mar 2022
    !===================================================================================
    ! Only 1 filter is accessiable for RF. Set up different gaussian factor for multi-scale inversion
    if (trim(simu_type) == 'rf') then
      !==== Get P arrival time ====
      call get_rf_times(glob_sem_disp, ttp, tb, te)
      !==== Calculate syn RF and pre-processing ====
      if(my_nrec_local > 0) then
        do irec_local = 1, my_nrec_local
          irec = myrank*nrec_local_max+irec_local
          if (ELASTIC_SIMULATION) then
            call pre_proc_rf_elastic(ievt, irec, iband, glob_sem_disp, ttp,&
                                        baz_all, glob_dat_tw, glob_syn_tw, glob_syn)
          endif
        enddo
      endif
      !==== collect and broadcast data =======
      call synchronize_all()
      call collect_stf_deconvolution(glob_stnm, glob_dat_tw, glob_syn_tw, glob_ff, glob_syn, nrec_local_max, my_nrec_local)
      call synchronize_all()
      call bcast_all_cr(glob_dat_tw,nrec*NSTEP*NRCOMP)
      call bcast_all_cr(glob_syn_tw,nrec*NSTEP*NRCOMP)
      call bcast_all_cr(glob_syn,nrec*NSTEP*NRCOMP)
      !==== Kai added to obtain a robust average amplitude scaling factor, R for RFs ====
      call average_amp_scale(glob_dat_tw, 2, avgamp)
      if(myrank==0) write(*,*) 'avgamp of data gather:', avgamp
      !==== Measure adjoint source for RFs ====
      if(my_nrec_local > 0) then
        do irec_local = 1, my_nrec_local
          irec = myrank*nrec_local_max+irec_local
          if (ELASTIC_SIMULATION) then
            if ((NSTEP * dt - t0 - ttp(irec))<te(irec)) stop 'The TW_AFTER is too large'
            call meas_adj_rf(ievt, irec, iband, bandname, glob_dat_tw, glob_syn_tw, glob_syn,&
                              tb(irec),te(irec), avgamp, ttp, window_chi, total_misfit, &
                              glob_net, glob_sta, glob_chan_dat, glob_tstart, glob_tend, &
                              glob_window_chi, glob_tr_chi, glob_am_chi, glob_num_win)
          endif
        enddo ! end loop irec
      endif ! end nrec >0
    endif ! end trim(simu_type) == 'rf'
    !=================== end of meas_adj for RFFWI =======================================
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !================= Calculate STF for TeleFWI =======================================
    if(trim(simu_type)=="tele") then
        !write(*,*) "Rank ",myrank,' has ',nrec_local,' receivers:'
      call cal_fktimes(ievt, ttp, tb, te)
      if (my_nrec_local > 0.) then
        do irec_local = 1, my_nrec_local
          irec = myrank*nrec_local_max+irec_local
          if (ELASTIC_SIMULATION) then
            call pre_proc_tele_elastic(ievt, irec, glob_sem_disp, tb(irec),te(irec), fstart0, fend0, bandname, &
                                       baz_all, glob_stnm, glob_dat_tw, glob_syn_tw, glob_ff, ttp)
          endif
          if (ACOUSTIC_SIMULATION) then
              if (ispec_is_acoustic(ispec)) then
              !! compute adjoint source according to cost L2 function
                write(*,*) 'To be added '
                stop
              endif
          endif
        enddo ! end loop irec
      endif ! my_nrec_local > 0.
      call synchronize_all()
      !---------------- collect info from different procs -------------------
      call collect_stf_deconvolution(glob_stnm,glob_dat_tw,glob_syn_tw,glob_ff,glob_syn,nrec_local_max,my_nrec_local)
      if (myrank == 0) then
        stf_array=0.
        ! calculate STF by deconvolution
        if (tele_par%USE_LOCAL_STF) then
          call read_local_stf(ievt, stf_array)
        else
          call seis_pca(glob_stnm,glob_dat_tw,glob_syn_tw,glob_ff,nrec,NSTEP,real(-t0),real(DT),stf_array)  
        endif
        write(*,*)'Finished calculating teleseismic STF.'
      endif
      call synchronize_all()
      call bcast_all_cr(stf_array,NSTEP*NRCOMP)
      call bcast_all_cr(glob_dat_tw,nrec*NSTEP*NRCOMP)
      call bcast_all_cr(glob_syn_tw,nrec*NSTEP*NRCOMP)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      !!! Kai added to obtain a robust average amplitude scaling factor
      call average_amp_scale(glob_dat_tw, 1, avgamp)
      if(myrank==0) write(*,*) 'avgamp of data gather:', avgamp
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! measure 
      if(my_nrec_local > 0) then
        do irec_local = 1, my_nrec_local
          irec = myrank*nrec_local_max+irec_local
          if (ELASTIC_SIMULATION) then
            call meas_adj_tele(ievt, irec, bandname, glob_dat_tw, glob_syn_tw,tb(irec),te(irec), &
                              stf_array, avgamp, ttp, window_chi,tr_chi,am_chi,total_misfit,&
                              glob_file_prefix0, glob_net, glob_sta,  glob_chan_dat, glob_tstart, &
                              glob_tend, glob_window_chi, glob_tr_chi, glob_am_chi,  glob_T_pmax_dat, &
                              glob_T_pmax_syn, glob_imeas, glob_num_win)
          endif
        enddo ! end loop irec
      endif ! end my_nrec_local >0
    endif ! end simu_type == 'tele'
    !================= end of meas_adj for TeleFWI =======================================
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !======================= For Convolutional Difference FWI ============================
    if (trim(simu_type) == 'telecd') then
      write(*,*) "Rank ",myrank,' has ',my_nrec_local,' receivers:'
      !==== Get P arrival time ====
      call get_rf_times(glob_sem_disp, ttp, tb, te)
      if(my_nrec_local > 0) then
        do irec_local = 1, my_nrec_local
          irec = myrank*nrec_local_max+irec_local
          if (ELASTIC_SIMULATION) then
            call pre_proc_tele_cd_elastic(ievt, irec, glob_sem_disp, fstart0, fend0, bandname, &
                                          baz_all, glob_dat_tw, glob_syn_tw, ttp)
          endif
        enddo
      endif
      call synchronize_all()
      call get_zdat_pca(glob_dat_tw, glob_ff)
      call collect_stf_deconvolution(glob_stnm, glob_dat_tw, glob_syn_tw, glob_ff, glob_syn, nrec_local_max, my_nrec_local)
      call synchronize_all()
      call bcast_all_cr(glob_dat_tw,nrec*NSTEP*NRCOMP)
      call bcast_all_cr(glob_syn_tw,nrec*NSTEP*NRCOMP)
      !==== Kai added to obtain a robust average amplitude scaling factor of Z comp ====
      call average_amp_scale(glob_dat_tw, 1, avgamp)
      if(myrank==0) write(*,*) 'avgamp of data gather:', avgamp
      if(my_nrec_local > 0) then
        do irec_local = 1, my_nrec_local
          irec = myrank*nrec_local_max+irec_local
          if (ELASTIC_SIMULATION) then
            !=============================================================================================
            call meas_adj_tele_cd(ievt, irec, bandname,baz_all,ttp,tb(irec),te(irec),glob_dat_tw,&
                                glob_syn_tw, glob_ff, avgamp, window_chi,total_misfit, &
                                glob_net, glob_sta,glob_chan_dat, glob_tstart, &
                                glob_tend, glob_window_chi, glob_tr_chi, glob_am_chi, glob_num_win)
            !==============================================================================================
          endif
          if (ACOUSTIC_SIMULATION) then
              if (ispec_is_acoustic(ispec)) then
              !! compute adjoint source according to cost L2 function
                write(*,*) 'To be added '
                stop
              endif
          endif
        enddo ! end loop irec
      endif ! end nrec >0 for noise
    endif ! end trim(simu_type) == 'telecd'
    !================= end of meas_adj for TeleCDFWI =======================================
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !======================= For noise FWI ===============================================
    if(my_nrec_local > 0 .and. simu_type=="noise") then
      !write(*,*) "Rank ",myrank,' has ',nrec_local,' receivers:'
      do irec_local = 1, my_nrec_local
        irec= myrank*nrec_local_max+irec_local
        
        if (ELASTIC_SIMULATION) then
          !====================================================================================================
          call meas_adj_noise(ievt, irec, glob_sem_disp, iband, fstart0, fend0, bandname, &
                              baz_all, glob_dat_tw, glob_syn_tw, &
                              window_chi,tr_chi,am_chi,total_misfit,&
                              glob_file_prefix0, glob_net, glob_sta,  glob_chan_dat, glob_tstart, &
                              glob_tend, glob_window_chi, glob_tr_chi, glob_am_chi,  glob_T_pmax_dat, &
                              glob_T_pmax_syn, glob_imeas, glob_num_win)
          !====================================================================================================
        endif
        if (ACOUSTIC_SIMULATION) then
          if (ispec_is_acoustic(ispec)) then

          !! compute adjoint source according to cost L2 function
            write(*,*) 'To be added '
            stop
          endif
        endif
      enddo ! end loop irec
    endif ! end nrec >0 for noise
    !====================== end of meas_adj noise FWI   ===========================================
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !======================= For local earthquake (leq) FWI ===============================================
    if(my_nrec_local > 0 .and. simu_type=="leq") then
      !write(*,*) "Rank ",myrank,' has ',my_nrec_local,' receivers for simu_type: ',trim(simu_type)
      do irec_local = 1, my_nrec_local
        irec= myrank*nrec_local_max+irec_local
        if (ELASTIC_SIMULATION) then
          !====================================================================================================
          call meas_adj_leq(icmt, ievt, irec, glob_sem_disp, iband, fstart0, fend0, bandname, &
                            baz_all, glob_dat_tw, glob_syn_tw, &
                            window_chi,tr_chi,am_chi,total_misfit,&
                            glob_file_prefix0, glob_net, glob_sta,  glob_chan_dat, glob_tstart, &
                            glob_tend, glob_window_chi, glob_tr_chi, glob_am_chi,  glob_T_pmax_dat, &
                            glob_T_pmax_syn, glob_imeas, glob_num_win)      
          !====================================================================================================
        endif
        if (ACOUSTIC_SIMULATION) then
          if (ispec_is_acoustic(ispec)) then
            !! compute adjoint source according to cost L2 function
            write(*,*) 'To be added '
            stop
          endif
        endif
      enddo ! end loop irec
    endif ! end nrec >0 for noise
    !====================== end of meas_adj leq FWI   ===========================================
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !======================= collect and write ============================
    if (icmt==0) then
      total_misfit_reduced=0.
      call sum_all_all_cr(total_misfit, total_misfit_reduced)
      if(myrank==0) then
        write(*,*) 'Total misfit is ',total_misfit_reduced
        write(OUT_FWAT_LOG,*) 'Total misfit of ', trim(bandname), ' is ',total_misfit_reduced
        flush(OUT_FWAT_LOG)
      endif
      
      call synchronize_all()
      ! glob_file_prefix0 seems not to be used
      call collect_chi(glob_file_prefix0, glob_net, glob_sta, glob_chan_dat, glob_tstart, &
                        glob_tend, glob_window_chi, glob_tr_chi, glob_am_chi, glob_T_pmax_dat, &
                        glob_T_pmax_syn, glob_imeas, glob_num_win, nrec_local_max, my_nrec_local)
      call synchronize_all()
      ! write out to one file
      if (myrank==0) then
        write(*,*) 'write window_chi to file for ',trim(bandname),' ...'
        open(14,file='src_rec/STATIONS_'//trim(acqui_par%evtid_names(ievt)),iostat=ier)
        open(15,file='src_rec/STATIONS_'//trim(acqui_par%evtid_names(ievt))//'_ADJOINT',iostat=ier)
        do irec=1,nrec
          read(14,"(a)",iostat=ier) dummystring
          dummystring = trim(dummystring)
          read(dummystring, *) sta_name, net_name, stlat, stlon, stele, stbur
          do icomp=1,NRCOMP
            if (leq_par%CMT3D_INV) then
              dummystring=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
                //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
                //'.obs.sac'//'.'//trim(bandname)
              if ( glob_num_win(irec,icomp)>0 ) then
                write(win_fileid,'(a)') trim(dummystring)
              endif
              dummystring=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
                //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
                //'.syn.sac'//'.'//trim(bandname)
              if ( glob_num_win(irec,icomp)>0 ) then
                write(win_fileid,'(a)') trim(dummystring)
                write(win_fileid,'(I10)') glob_num_win(irec,icomp)
              endif
            endif ! CMT3D_INV
            if ( glob_num_win(irec,icomp)>0 ) then
              do iwin=1,glob_num_win(irec,icomp)
                if (glob_tend(irec,icomp,iwin)>0.00001) then
                  write(chi_fileid,'(a20,a8,a4,a5,i4,i4,2e14.6,20e14.6,2e14.6,2f14.6)') &
                    acqui_par%evtid_names(ievt),glob_sta(irec),glob_net(irec),&
                    glob_chan_dat(irec,icomp),irec,glob_imeas(irec,icomp,iwin),&
                    glob_tstart(irec,icomp,iwin),glob_tend(irec,icomp,iwin),glob_window_chi(irec,icomp,:,iwin), &
                    glob_tr_chi(irec,icomp,iwin),glob_am_chi(irec,icomp,iwin),glob_T_pmax_dat(irec,icomp,iwin), &
                    glob_T_pmax_syn(irec,icomp,iwin)
                  if (leq_par%CMT3D_INV) write(win_fileid,'(2f11.4)') glob_tstart(irec,icomp,iwin),glob_tend(irec,icomp,iwin) 
                endif
              enddo ! iwin
              num_adj=1
              flush(chi_fileid)
            endif ! nwin > 0
          enddo ! icomp
          if (num_adj>0) then ! if only there is an adjoint source in any component, write it to STATIONS_ADJOINT 
              write(15,'(a10,1x,a10,4e18.6)') &
                      trim(sta_name),trim(net_name), &
                      sngl(stlat),sngl(stlon),sngl(stele),sngl(stbur)
          endif
        enddo ! irec  
        close(14)
        close(15)
      endif
    endif ! end if icmt
  enddo ! end loop iband
  !======================================================================
  deallocate(glob_file_prefix0)
  deallocate(glob_net)
  deallocate(glob_sta) 
  deallocate(glob_chan_dat)   
  deallocate(glob_num_win)  
  deallocate(glob_tstart)
  deallocate(glob_tend)  
  deallocate(glob_window_chi)
  deallocate(glob_tr_chi)
  deallocate(glob_am_chi)
  deallocate(glob_T_pmax_dat)
  deallocate(glob_T_pmax_syn)
  deallocate(glob_imeas) 


  !***************** sum adjoint sources over bands **********************
  if(icmt==0) then
    if (myrank==0) write(*,*) 'sum and rotate adjoint sources to ZNE ...'
    if(my_nrec_local > 0) then
      do irec_local = 1, my_nrec_local
        irec=myrank*nrec_local_max+irec_local
        call sum_adj_source(ievt, irec, simu_type, adj_syn_all_sum, baz_all, stf_array, avgamp, num_adj)       
      enddo
    endif ! end if nrec_local
  endif
  call synchronize_all()
  if(myrank==0) write(*,*) 'Finished preprocessing here.'
  
end subroutine run_preprocessing

subroutine dif1(signal, dt, nt)
  use ma_constants, only : NDIM
  implicit none
  double precision, intent(in)                      :: dt
  integer, intent(in)                               :: nt
  double precision, dimension(NDIM), intent(inout)  :: signal
  double precision, dimension(NDIM)                 :: temp_signal
  integer                                           :: i

  temp_signal(1:nt)=signal(1:nt)

  do i=2,nt-1
    signal(i) = 0.5 * ( temp_signal(i+1) - temp_signal(i-1) ) / dt
  enddo

  signal(1)=0.
  signal(nt)=0.


end subroutine dif1

!================================================================================

!===============================================
!  subroutine detrend(array,npt1,t01,dt1)
!  use ma_constants, only : NDIM
!  implicit none
!  double precision :: array(NDIM) 
!  double precision t01,dt1,t
!  integer :: npt1,i,ier
!  real                :: a,b,siga,sigb,sig,cc
!  
!  call lifite(t01,dt1,array,npt1,a,b,siga,sigb,sig,cc) 
!  !write(*,*)'slop and its STD: ', a,siga
!  !write(*,*)'intercept and its STD: ', b,sigb
!  !write(*,*)'data STD: ', sig
!  !write(*,*)'corr coef: ', cc
!  !write(*,*)'npt1=',npt1
!
!
!  do i=1,npt1
!      t=t01+(i-1)*dt1
!      array(i)=array(i)-b-a*t
!  enddo
!  end subroutine detrend


