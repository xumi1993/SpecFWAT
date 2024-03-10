module preproc_subs 
  use specfem_par, only: cr => CUSTOM_REAL, MAX_LEN => MAX_STRING_LEN,&
                              network_name, station_name,nrec,myrank,&
                              t0,NSTEP,DT,IIN,OUTPUT_FILES
  use fullwave_adjoint_tomo_par
  use fwat_input
  use fwat_utils, only: rotate_ZNE_to_ZRT
  use measure_adj_mod
  use sacio

  implicit none
  integer                                                :: ier,irec
contains

subroutine read_fktimes(ievt, ttp, tb, te)
  character(len=MAX_LEN)                             :: datafile, dummystring
  character(len=256),dimension(:), allocatable       :: fk_netwk,fk_stnm
  real(kind=4) , dimension(:), allocatable           :: fk_ttp,fk_tb,fk_te 
  integer                                            :: fk_irec,fk_nrec,ievt
  real                                               :: dummy
  real(kind=4)                                       :: ttp(nrec),tb(nrec),te(nrec)

  if (myrank == 0) then
    datafile='src_rec/FKtimes_'//trim(acqui_par%evtid_names(ievt))
    open(unit=IIN,file=trim(datafile),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(datafile))
        ! reads all stations
    fk_nrec=0
    do
      read(IIN,*,iostat=ier) dummystring,dummystring,dummy
      if (ier /= 0) exit
      fk_nrec=fk_nrec+1
    enddo
    ! close receiver file
    close(IIN)
    allocate(fk_netwk(fk_nrec)) 
    allocate(fk_stnm(fk_nrec)) 
    allocate(fk_ttp(fk_nrec)) 
    allocate(fk_tb(fk_nrec)) 
    allocate(fk_te(fk_nrec)) 
    open(unit=IIN,file=trim(datafile),status='old',action='read',iostat=ier)
    do irec=1,fk_nrec
      if (tele_par%TW_AFTER.eq.0.) then
        if (irec==1) write(*,*) 'Read TW_BEFORE and TW_AFTER from FKtimes...'
        read(IIN,*,iostat=ier) fk_netwk(irec),fk_stnm(irec),fk_ttp(irec),fk_tb(irec),fk_te(irec)
      else
        if (irec==1) write(*,*) 'Read TW_BEFORE and TW_AFTER from FWAT.PAR'
        read(IIN,*,iostat=ier) fk_netwk(irec),fk_stnm(irec),fk_ttp(irec)
      endif
    enddo
    close(IIN)
    ! find the right ttp among all stations (Note STATIONS_FILTERED might be
    ! less than STATIONS/FKtimes, thus we need to search for the right ttp )
    do irec=1,nrec
      do fk_irec=1,fk_nrec
        if (trim(network_name(irec))==trim(fk_netwk(fk_irec)) .and. &
          trim(station_name(irec))==trim(fk_stnm(fk_irec)) ) then
          ttp(irec)=fk_ttp(fk_irec)
          if (tele_par%TW_AFTER.eq.0.) then
            tb(irec)=fk_tb(fk_irec)
            te(irec)=fk_te(fk_irec)
          else
            tb(irec)=tele_par%TW_BEFORE
            te(irec)=tele_par%TW_AFTER
          endif
          write(*,*)'netwk,stnm,ttp,tb,te=',trim(network_name(irec)),'.',trim(station_name(irec)), &
                                              ttp(irec),tb(irec),te(irec)
          if( (ttp(irec)-t0+te(irec)) > (-t0+(NSTEP-1)*DT) ) then
            write(*,*)'ttp exceed data range'
            stop
          endif
        endif
      enddo
    enddo
    deallocate(fk_netwk)
    deallocate(fk_stnm)
    deallocate(fk_ttp)
    deallocate(fk_tb)
    deallocate(fk_te)
  endif
  call bcast_all_cr(ttp,nrec)
  call bcast_all_cr(tb,nrec)
  call bcast_all_cr(te,nrec)
end subroutine

subroutine get_rf_times(glob_sem_disp, ttp, tb, te)
  real(kind=4)                                     :: ttp(nrec),tb(nrec),te(nrec)
  integer                                          :: irec, mloc(1)
  real(kind=cr), dimension(nrec,NSTEP,3)           :: glob_sem_disp
  real(kind=cr), dimension(NSTEP)                  :: seismo_z
  ! real(kind=4), dimension(:), allocatable          :: corr1, corr2


  if (myrank == 0) then
    do irec = 1, nrec
      !---- Need lowpass filter here? -----
      seismo_z(:)=glob_sem_disp(irec,:,3)
      ! allocate(corr1(NSTEP))
      ! allocate(corr2(NSTEP))
      ! call mycorrelation(seismo_z,seismo_z,NSTEP,NSTEP,corr1,0)
      ! call mycorrelation(seismo_z,corr1,NSTEP,NSTEP,corr2,0)
      mloc = maxloc(seismo_z)
      ttp(irec) = mloc(1) * dt - t0
      if (ttp(irec) + t0 < rf_par%rf_tshift) stop 'Time length before P arrival must be greater than RF_TSHIFT.'
      write(*,*) trim(network_name(irec))//'.'//trim(station_name(irec)),', RF times=',ttp(irec)
      tb(irec)=tele_par%TW_BEFORE
      te(irec)=tele_par%TW_AFTER
      ! deallocate(corr1)
      ! deallocate(corr2)
    enddo
  endif
  call bcast_all_cr(ttp,nrec)
  call bcast_all_cr(tb,nrec)
  call bcast_all_cr(te,nrec)
end subroutine

subroutine read_local_stf(ievt, stf_array)
  implicit none
  double precision, dimension(NDIM)                :: datarray
  double precision                                 :: t01,dt1
  real(kind=CUSTOM_REAL), dimension(NSTEP, NRCOMP) :: stf_array
  character(len=MAX_STRING_LEN)                    :: fname
  integer                                          :: ievt, icomp,npt1
  
  do icomp=1,NRCOMP
    datarray = 0.
    fname = 'src_rec/STF_'//trim(acqui_par%evtid_names(ievt))//&
            '_'//trim(RCOMPS(icomp))//'.sac'
    print *, trim(fname)
    call drsac1(trim(fname),datarray,npt1,t01,dt1)
    stf_array(1:NSTEP,icomp) = real(datarray(1:NSTEP))
  enddo

end subroutine read_local_stf

subroutine pre_proc_tele_cd_elastic(ievt, irec, glob_sem_disp, fstart0, fend0, bandname, &
                                     baz_all, glob_dat_tw, glob_syn_tw)
  implicit None

  integer                                          :: ievt, irec, NDIM_CUT
  logical                                          :: findfile
  character(len=MAX_STRING_LEN)                    :: datafile,adjfile,bandname
  character(len=10)                                :: net,sta,chan_dat
  integer                                          :: yr,jda,ho,mi, npt1
  double precision                                 :: t01,dt1,sec,dist,az,baz,slat,slon,&
                                                      fstart0,fend0
  real(kind=cr), dimension(3,NSTEP)                :: seismo_syn
  double precision, dimension(nrec)                :: baz_all
  real(kind=4), dimension(nrec,NSTEP,NRCOMP)       :: glob_dat_tw,glob_syn_tw
  real(kind=cr), dimension(nrec,NSTEP,3)           :: glob_sem_disp
  double precision, dimension(NDIM)                :: datr_inp,synr_inp,datz_inp,synz_inp,&
                                                      datr_inp_bp,synr_inp_bp,&
                                                      datz_inp_bp,synz_inp_bp,datarray
  
  ! RCOMPS(1) = Z; RCOMPS(2) = R
  seismo_syn(:,:) = 0.
  ! Read Z component
  datafile='./'//trim(acqui_par%in_dat_path(ievt))//'/'//trim(network_name(irec))//'.'&
          //trim(station_name(irec))//'.'//trim(CH_CODE)//'Z.sac' 
  inquire(file=trim(datafile),exist=findfile)
  if ( findfile ) then
    call drsac1(trim(datafile),datarray,npt1,t01,dt1)
    call get_sacfile_header(trim(datafile),yr,jda,ho,mi,sec,net,sta, &
                              chan_dat,dist,az,baz,slat,slon)
    baz_all(irec)=baz
  endif
  datz_inp(1:npt1)=datarray(1:npt1)

  ! Rotation
  if (trim(dat_coord)=='ZRT') then
    call rotate_ZNE_to_ZRT(glob_sem_disp(irec,:,3),glob_sem_disp(irec,:,2),&
    glob_sem_disp(irec,:,1),seismo_syn(1,:),seismo_syn(2,:),seismo_syn(3,:),NSTEP,real(baz)) 
  else
    seismo_syn(1,:)=glob_sem_disp(irec,:,3)
    seismo_syn(2,:)=glob_sem_disp(irec,:,2)
    seismo_syn(3,:)=glob_sem_disp(irec,:,1)
  endif
  synz_inp(1:NSTEP)=dble(seismo_syn(1,1:NSTEP))
  synr_inp(1:NSTEP)=dble(seismo_syn(2,1:NSTEP))

  ! Read component 2
  datafile='./'//trim(acqui_par%in_dat_path(ievt))//'/'//trim(network_name(irec))//'.'&
          //trim(station_name(irec))//'.'//trim(CH_CODE)//'R.sac' 
  inquire(file=trim(datafile),exist=findfile)
  if ( findfile ) then
    call drsac1(trim(datafile),datarray,npt1,t01,dt1)
  endif
  datr_inp(1:npt1)=datarray(1:npt1)

  ! rtrend
  NDIM_CUT=NSTEP
  call detrend(datr_inp,NDIM_CUT)
  call detrend(datz_inp,NDIM_CUT)
  call detrend(synr_inp,NDIM_CUT)
  call detrend(synz_inp,NDIM_CUT)
  ! remean
  datz_inp(1:NDIM_CUT)=datz_inp(1:NDIM_CUT)-sum(datz_inp(1:NDIM_CUT))/NDIM_CUT
  datr_inp(1:NDIM_CUT)=datr_inp(1:NDIM_CUT)-sum(datr_inp(1:NDIM_CUT))/NDIM_CUT
  synz_inp(1:NDIM_CUT)=synz_inp(1:NDIM_CUT)-sum(synz_inp(1:NDIM_CUT))/NDIM_CUT
  synr_inp(1:NDIM_CUT)=synr_inp(1:NDIM_CUT)-sum(synr_inp(1:NDIM_CUT))/NDIM_CUT
  ! filter
  datz_inp_bp=datz_inp
  datr_inp_bp=datr_inp
  synz_inp_bp=synz_inp
  synr_inp_bp=synr_inp
  call bandpass(datz_inp_bp,NDIM_CUT,DT,fstart0,fend0)
  call bandpass(datr_inp_bp,NDIM_CUT,DT,fstart0,fend0)
  call bandpass(synz_inp_bp,NDIM_CUT,DT,fstart0,fend0)
  call bandpass(synr_inp_bp,NDIM_CUT,DT,fstart0,fend0)
  if (VERBOSE_MODE) then
    adjfile=trim(OUTPUT_FILES)//'/syn.'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)&
            //'Z.sac'//'.'//trim(bandname) 
    call dwsac1(trim(adjfile),synz_inp_bp,NDIM_CUT,dble(-T0),dble(DT))
    adjfile=trim(OUTPUT_FILES)//'/dat.'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)&
            //'Z.sac'//'.'//trim(bandname) 
    call dwsac1(trim(adjfile),datz_inp_bp,NDIM_CUT,dble(-T0),dble(DT))
    adjfile=trim(OUTPUT_FILES)//'/syn.'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)&
            //'R.sac'//'.'//trim(bandname) 
    call dwsac1(trim(adjfile),synr_inp_bp,NDIM_CUT,dble(-T0),dble(DT))
    adjfile=trim(OUTPUT_FILES)//'/dat.'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)&
            //'R.sac'//'.'//trim(bandname) 
    call dwsac1(trim(adjfile),datr_inp_bp,NDIM_CUT,dble(-T0),dble(DT))
  endif
  glob_dat_tw(irec,1:NDIM_CUT,1)=datz_inp_bp(1:NDIM_CUT)
  glob_dat_tw(irec,1:NDIM_CUT,2)=datr_inp_bp(1:NDIM_CUT)
  glob_syn_tw(irec,1:NDIM_CUT,1)=synz_inp_bp(1:NDIM_CUT)
  glob_syn_tw(irec,1:NDIM_CUT,2)=synr_inp_bp(1:NDIM_CUT)

end subroutine

subroutine pre_proc_tele_elastic(ievt, irec, glob_sem_disp, win_tb,win_te, fstart0, fend0, bandname, &
                                 baz_all, glob_stnm, glob_dat_tw, glob_syn_tw, glob_ff, ttp)

  use telestf_mod

  integer                                          :: icomp, ievt, irec, NDIM_CUT
  logical                                          :: findfile
  double precision                                 :: t01,dt1, fstart0, fend0
  integer                                          :: yr,jda,ho,mi, npt1
  double precision                                 :: t0_inp,t1_inp,dt_inp
  double precision, dimension(NDIM)                :: datarray
  double precision, dimension(NDIM)                :: dat_inp, syn_inp
  double precision, dimension(NDIM)                :: dat_inp_bp, syn_inp_bp
  real(kind=cr), dimension(nrec,NSTEP,3)           :: glob_sem_disp
  real(kind=cr)                                    :: seismo_syn(3,NSTEP)                           
  real(kind=4), dimension(NSTEP)                   :: one_seismo_dat,one_seismo_syn
  double precision                                 :: sec,dist,az,baz,slat,slon
  character(len=10)                                :: net,sta,chan_dat
  character(len=MAX_LEN)                           :: bandname,adjfile,datafile
  double precision                                 :: baz_all(nrec)
  real                                             :: win_tb,win_te
  character(len=256), dimension(nrec)              :: glob_stnm
  real(kind=4), dimension(nrec,NSTEP,NRCOMP), intent(inout)       :: glob_dat_tw,glob_syn_tw
  real(kind=4), dimension(nrec,NSTEP)              :: glob_ff
  real(kind=4)                                     :: ttp(nrec)

  seismo_syn(:,:)=0.d0 !glob_sem_disp(irec,:,:)

  do icomp=1,NRCOMP
    datafile='./'//trim(acqui_par%in_dat_path(ievt))//'/'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))//'.sac' 
    !write(*,*)'myrank,datafile=',myrank,trim(datafile)
    inquire(file=trim(datafile),exist=findfile)
    if ( .not. findfile ) then
      stop 'No such sac file in the fwat_data'
    else
      call drsac1(trim(datafile),datarray,npt1,t01,dt1)
      call get_sacfile_header(trim(datafile),yr,jda,ho,mi,sec,net,sta, &
                              chan_dat,dist,az,baz,slat,slon)
      
      baz_all(irec)=baz
      !write(*,*)'myrank,datafile,npt1,t01,dt1,NSTEP,t0,DT,dist,ttp= ',myrank,trim(datafile),npt1,t01,dt1,&
      !                                   NSTEP,t0,DT,dist,ttp(irec)
      ! if(abs(dt1-dT)>0.0001) stop 'delta of data does NOT match DT of SEM synthetics'
      if(abs(dt1-dT)>0.0001 .or. npt1 /= NSTEP) then
        call interpolate_syn(datarray, t01,dt1,npt1,-dble(t0),dt,NSTEP)
      endif
      ! only do the rotation for one time
      if (icomp.eq.1 ) then
        if (trim(dat_coord)=='ZRT') then
          call rotate_ZNE_to_ZRT(glob_sem_disp(irec,:,3),glob_sem_disp(irec,:,2),&
          glob_sem_disp(irec,:,1),seismo_syn(1,:),seismo_syn(2,:),seismo_syn(3,:),NSTEP,real(baz)) 
        else
          seismo_syn(1,:)=glob_sem_disp(irec,:,3)
          seismo_syn(2,:)=glob_sem_disp(irec,:,2)
          seismo_syn(3,:)=glob_sem_disp(irec,:,1)
        endif
      endif
      !****************  pre-process dat and syn ****************************
      dt_inp=DT !
      t0_inp=-t0 
      t1_inp=-t0+(NSTEP-1)*DT
      NDIM_CUT=NSTEP !(t1_inp-t0_inp)/dt_inp+1
      dat_inp(:)=0.
      syn_inp(:)=0.
      dat_inp(1:npt1)=datarray(1:npt1)
      syn_inp(1:NSTEP)=dble(seismo_syn(icomp,1:NSTEP))
      !!! Filter
      ! rtrend
      !call detrend(dat_inp,NDIM_CUT,t0_inp,dt_inp)
      !call detrend(dat_inp,NDIM_CUT,t0_inp,dt_inp)
      call detrend(dat_inp,NDIM_CUT)
      call detrend(syn_inp,NDIM_CUT)
      ! remean
      dat_inp(1:NDIM_CUT)=dat_inp(1:NDIM_CUT)-sum(dat_inp(1:NDIM_CUT))/NDIM_CUT
      syn_inp(1:NDIM_CUT)=syn_inp(1:NDIM_CUT)-sum(syn_inp(1:NDIM_CUT))/NDIM_CUT
      dat_inp_bp=dat_inp
      syn_inp_bp=syn_inp
      call bandpass(dat_inp_bp,NDIM_CUT,dt_inp,fstart0,fend0)
      call bandpass(syn_inp_bp,NDIM_CUT,dt_inp,fstart0,fend0)
      ! rtrend
      !call detrend(dat_inp_bp,NDIM_CUT,t0_inp,dt_inp)
      !call detrend(syn_inp_bp,NDIM_CUT,t0_inp,dt_inp)
      call detrend(dat_inp_bp,NDIM_CUT)
      call detrend(syn_inp_bp,NDIM_CUT)
      ! remean
      dat_inp_bp(1:NDIM_CUT)=dat_inp_bp(1:NDIM_CUT)-sum(dat_inp_bp(1:NDIM_CUT))/NDIM_CUT
      syn_inp_bp(1:NDIM_CUT)=syn_inp_bp(1:NDIM_CUT)-sum(syn_inp_bp(1:NDIM_CUT))/NDIM_CUT

      if (VERBOSE_MODE) then
        adjfile=trim(OUTPUT_FILES)//'/syn.'//trim(network_name(irec))//'.'&
                //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
                //'.sac'//'.'//trim(bandname) 
        call dwsac1(trim(adjfile),syn_inp_bp,NDIM_CUT,t0_inp,dt_inp)
        adjfile=trim(OUTPUT_FILES)//'/dat.'//trim(network_name(irec))//'.'&
                //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
                //'.sac'//'.'//trim(bandname) 
        call dwsac1(trim(adjfile),dat_inp_bp,NDIM_CUT,t0_inp,dt_inp)
      endif

      one_seismo_dat(1:NDIM_CUT)=dat_inp_bp(1:NDIM_CUT)
      glob_dat_tw(irec,1:NDIM_CUT,icomp)=dat_inp_bp(1:NDIM_CUT)
      one_seismo_syn(1:NDIM_CUT)=syn_inp_bp(1:NDIM_CUT)
      glob_syn_tw(irec,1:NDIM_CUT,icomp)=syn_inp_bp(1:NDIM_CUT)
      glob_stnm(irec)=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
                      //trim(station_name(irec))//'.'//trim(CH_CODE)& 
                      //'.sac'//'.'//trim(bandname)
      if (icomp==1)  then !!! get only vertical deconved traces
        write(*,*)'run time_deconv for ',trim(glob_stnm(irec))//' on window: ',win_tb,win_te
        call time_iterdeconv(glob_stnm,one_seismo_dat,one_seismo_syn,glob_dat_tw,glob_syn_tw,glob_ff,&
                              ttp(irec),win_tb,win_te,irec,nrec,NSTEP,real(-t0),real(DT),fstart0,fend0)
      endif
    endif ! end findfile
  enddo ! end icomp  
end subroutine pre_proc_tele_elastic


subroutine average_amp_scale(glob_dat_tw, icomp, avgamp)
  real, intent(inout)                              :: avgamp
  real                                             :: avgamp0
  integer                                          :: igood, icomp
  real(kind=4), dimension(nrec,NSTEP,NRCOMP)       :: glob_dat_tw

  ! use only Z component for amplitude scale
  avgamp0=0.
  do irec =1 ,nrec
    avgamp0=avgamp0+maxval(abs(glob_dat_tw(irec,:,icomp))) 
  enddo
  avgamp0=avgamp0/nrec
  avgamp=0
  igood=0
  do irec =1, nrec
    if ((maxval(abs(glob_dat_tw(irec,:,icomp)))-avgamp0)<0.2*avgamp0) then
      avgamp=avgamp+maxval(abs(glob_dat_tw(irec,:,icomp)))
      igood=igood+1
    endif
  enddo
  avgamp=avgamp/igood
end subroutine average_amp_scale


subroutine pre_proc_rf_elastic(ievt, irec, igaus, glob_sem_disp,&
                               ttp, baz_all, glob_dat_tw, glob_syn_tw, glob_syn)
  use decon_mod

  type(sachead)                                    :: head
  real(kind=cr), dimension(nrec,NSTEP,3)           :: glob_sem_disp
  double precision, dimension(nrec)                :: baz_all
  real(kind=cr), dimension(3,NSTEP)                :: seismo_syn
  integer                                          :: ievt, irec, igaus, icomp
  character(len=MAX_LEN)                           :: bandname,adjfile,datafile
  real(kind=4), dimension(nrec,NSTEP,NRCOMP)       :: glob_dat_tw,glob_syn_tw,glob_syn
  double precision                                 :: t01,dt1,sec,dist,az,baz,slat,slon
  integer                                          :: yr,jda,ho,mi,npt1,nshift, flag
  character(len=10)                                :: net,sta,chan_dat
  double precision, dimension(NSTEP)               :: rfi!, zrf
  double precision, dimension(NDIM)                :: datarray, inp_r, inp_z, inp_r_bp, inp_z_bp
  logical                                          :: findfile
  real(kind=4)                                     :: ttp(nrec)

  ! seismo_syn(:,:) = glob_sem_disp(irec, :, :)
  seismo_syn(:,:)=0.d0 
  icomp = 2
  write(bandname,'(a1,f3.1)') 'F',rf_par%f0(igaus)
  datafile='./'//trim(acqui_par%in_dat_path(ievt))//'/'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//'R.'&
              //trim(bandname)//'.rf.sac'
  inquire(file=trim(datafile),exist=findfile)
  if ( findfile ) then
    write(*,*) 'Read rf ',trim(datafile)
    call drsac1(trim(datafile),datarray,npt1,t01,dt1)
    call get_sacfile_header(trim(datafile),yr,jda,ho,mi,sec,net,sta, &
                                chan_dat,dist,az,baz,slat,slon)
    ! Resample RF data to the same dimension as the syn RF.
    ! NOTE: the begin time (b) in RF data should be the -time_shift.
    if(abs(dt1-dT)>0.0001 .or. npt1 /= NSTEP) then
      call interpolate_syn(datarray, t01,dt1,npt1,-dble(rf_par%RF_TSHIFT),dt,NSTEP)
    endif
    nshift = floor((ttp(irec)-rf_par%RF_TSHIFT)/dt)
    glob_dat_tw(irec, 1:NSTEP, icomp) = datarray(1:NSTEP)
    baz_all(irec) = baz
    call rotate_ZNE_to_ZRT(glob_sem_disp(irec,:,3),glob_sem_disp(irec,:,2),&
                glob_sem_disp(irec,:,1),seismo_syn(1,:),seismo_syn(2,:),seismo_syn(3,:),NSTEP,real(baz))
    write(*,*) 'rotate synthetics from ZNE to ZRT'
    inp_r(:) = 0.
    inp_z(:) = 0.
    inp_r(1:NSTEP) = dble(seismo_syn(2,:))
    inp_z(1:NSTEP) = dble(seismo_syn(1,:))
    inp_r_bp = inp_r
    inp_z_bp = inp_z
    call bandpass(inp_r_bp,NSTEP,dble(dt),dble(1/tele_par%LONG_P(1)),dble(1/tele_par%SHORT_P(1)))
    call bandpass(inp_z_bp,NSTEP,dble(dt),dble(1/tele_par%LONG_P(1)),dble(1/tele_par%SHORT_P(1)))
    if (VERBOSE_MODE) then
      call sacio_newhead(head, sngl(dt), NSTEP, sngl(-t0))
      head%knetwk = trim(network_name(irec))
      head%kstnm = trim(station_name(irec))
      head%kcmpnm = trim(CH_CODE)//'Z'
      adjfile=trim(OUTPUT_FILES)//'/syn.'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//'Z.sac'
      call sacio_writesac(adjfile, head, dble(seismo_syn(1,:)), flag)
      adjfile=trim(OUTPUT_FILES)//'/wsyn.'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//'Z.sac'
      call sacio_writesac(adjfile, head, inp_z_bp(1:NSTEP), flag)
      head%kcmpnm = trim(CH_CODE)//'R'
      adjfile=trim(OUTPUT_FILES)//'/syn.'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//'R.sac'
      call sacio_writesac(adjfile, head, dble(seismo_syn(2,:)), flag)
      adjfile=trim(OUTPUT_FILES)//'/wsyn.'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//'R.sac'
      call sacio_writesac(adjfile, head, inp_r_bp(1:NSTEP), flag)
    endif
    call deconit(inp_r_bp(1:NSTEP), inp_z_bp(1:NSTEP), NSTEP,&
                  real(DT), rf_par%rf_tshift, rf_par%f0(igaus),&
                  rf_par%maxit, rf_par%minderr, 0, rfi)
    write(*,*) 'Deconvolution finished'
    glob_syn_tw(irec, 1:NSTEP, icomp) = rfi(1:NSTEP)
    ! call deconit(inp_z_bp(1:NSTEP), inp_z_bp(1:NSTEP), NSTEP, real(dt), RF_TSHIFT, &
                ! f0(igaus), 10, 0.001, 0, zrf)
    ! glob_syn_tw(irec, 1:NSTEP, icomp) = glob_syn_tw(irec, 1:NSTEP, icomp)/maxval(zrf)
    ! glob_dat_tw(irec, 1:NSTEP, icomp) = glob_dat_tw(irec, 1:NSTEP, icomp)/maxval(zrf)
    if (VERBOSE_MODE) then
      call sacio_newhead(head, sngl(dt), NSTEP, -rf_par%rf_tshift)
      head%knetwk = trim(network_name(irec))
      head%kstnm = trim(station_name(irec))
      head%kcmpnm = trim(CH_CODE)//'R'
      adjfile=trim(OUTPUT_FILES)//'/dat.'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//'R.rf.sac'//'.'//trim(bandname)
      call sacio_writesac(adjfile, head, dble(glob_dat_tw(irec, :, icomp)), flag)
      ! call dwsac1(trim(adjfile),dble(glob_dat_tw(irec, :, icomp)),NSTEP,-dble(rf_tshift),dt)
      adjfile=trim(OUTPUT_FILES)//'/syn.'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//'R.rf.sac'//'.'//trim(bandname)
      call sacio_writesac(adjfile, head, dble(glob_syn_tw(irec, :, icomp)), flag)
    endif
    glob_syn(irec, 1:NSTEP, 1) = inp_z_bp(1:NSTEP)
    glob_syn(irec, 1:NSTEP, 2) = inp_r_bp(1:NSTEP)
  else
    print *, 'Cannot open '//trim(datafile)
    stop
  endif
end subroutine

end module
