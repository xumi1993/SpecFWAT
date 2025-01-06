module preproc_measure_adj_subs
  use specfem_par, only: cr => CUSTOM_REAL, MAX_LEN =>MAX_STRING_LEN, OUTPUT_FILES, IIN, &
                          network_name, station_name,nrec, nrec_local, &
                          ispec_selected_rec, number_receiver_global, myrank, t0,islice_selected_rec
  use shared_input_parameters, only: NSTEP, DT, NPROC, SUPPRESS_UTM_PROJECTION
  use fullwave_adjoint_tomo_par
  use fwat_input
  use my_mpi
  use fwat_utils
  use measure_adj_mod, MAX_NDIM => NDIM
  ! module from preproc_flexwin
  use seismo_variables 
  use interpolation_mod 
  use decon_mod
  use sacio

  implicit none

  real(kind=cr), parameter :: target_dt_noise=0.02, max_duration=20

contains

  subroutine meas_adj_tele(ievt, irec, bandname, glob_dat_tw, glob_syn_tw, win_tb, win_te, &
                         stf_array, avgamp, ttp, window_chi,tr_chi,am_chi,total_misfit,&
                         glob_file_prefix0, glob_net, glob_sta,  glob_chan_dat, glob_tstart, &
                         glob_tend, glob_window_chi, glob_tr_chi, glob_am_chi,  glob_T_pmax_dat, &
                         glob_T_pmax_syn, glob_imeas, glob_num_win)
    
    type(sachead)                                    :: head
    logical                                          :: findfile
    integer                                          :: yr,jda,ho,mi, npt1, irec, icomp, ievt, &
                                                        out_imeas, flag
    double precision                                 :: t01,dt1,sec,dist,az,baz,slat,slon,tstart,tend
    double precision, dimension(MAX_NDIM)                :: datarray, dat_inp_bp, syn_inp_bp, adj_syn_all
    character(len=256),dimension(nrec)               :: glob_net, glob_sta
    character(len=256),dimension(nrec,NRCOMP)        :: glob_chan_dat
    character(len=MAX_STRING_LEN)                    :: bandname,adjfile,datafile,file_prefix0
    character(len=10)                                :: net,sta,chan_dat
    real(kind=4), dimension(:), allocatable          :: tmpl
    real(kind=4), dimension(nrec,NSTEP,NRCOMP)       :: glob_dat_tw,glob_syn_tw
    real(kind=4), dimension(nrec)                    :: ttp
    real                                             :: avgamp, win_tb, win_te
    double precision, dimension(nrec,NRCOMP,1)       :: glob_tstart, glob_tend, glob_tr_chi, &
                                                        glob_am_chi, glob_T_pmax_dat, glob_T_pmax_syn
    integer, dimension(nrec,NRCOMP)                  :: glob_num_win
    double precision, dimension(nrec,NRCOMP,NCHI,1)  :: glob_window_chi
    character(len=256),dimension(nrec)               :: glob_file_prefix0
    integer         , dimension(nrec,NRCOMP,1)       :: glob_imeas
    double precision, dimension(NCHI)                :: window_chi
    double precision                                 :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn
    real(kind=CUSTOM_REAL)                           :: total_misfit,misfit
    real(kind=4) , dimension(NSTEP,NRCOMP)           :: stf_array


    do icomp=1,NRCOMP
      datafile='./'//trim(acqui_par%in_dat_path(ievt))//'/'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))//'.sac' 
      inquire(file=trim(datafile),exist=findfile)
      if ( findfile ) then
        call drsac1(trim(datafile),datarray,npt1,t01,dt1)
        call get_sacfile_header(trim(datafile),yr,jda,ho,mi,sec,net,sta, &
                                  chan_dat,dist,az,baz,slat,slon)
        dat_inp_bp=0.
        syn_inp_bp=0.
        dat_inp_bp(1:NSTEP)=dble(glob_dat_tw(irec,1:NSTEP,icomp))
        syn_inp_bp(1:NSTEP)=dble(glob_syn_tw(irec,1:NSTEP,icomp))
        allocate(tmpl(NSTEP))
        call myconvolution(glob_syn_tw(irec,:,icomp),stf_array(:,icomp),NSTEP,NSTEP,tmpl,0)
        syn_inp_bp(1:NSTEP)=dble(tmpl)*dt1
        deallocate(tmpl)
        if (VERBOSE_MODE) then
          call sacio_newhead(head, sngl(dt1), npt1, sngl(t01))
          head%knetwk = trim(network_name(irec))
          head%kstnm = trim(station_name(irec))
          head%kcmpnm = trim(CH_CODE)//trim(RCOMPS(icomp))
          adjfile=trim(OUTPUT_FILES)//'/wdat.'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
              //'.sac'//'.'//trim(bandname)
          call sacio_writesac(adjfile, head, dat_inp_bp(1:NSTEP), flag)
          adjfile=trim(OUTPUT_FILES)//'/wsyn.'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
              //'.sac'//'.'//trim(bandname) 
          call sacio_writesac(adjfile, head, syn_inp_bp(1:NSTEP), flag)
        endif
        !!**************** measure_adj ******************************************
        glob_num_win(irec,icomp)=1
        tstart=ttp(irec)-win_tb
        tend=ttp(irec)+win_te
        tstart = max(tstart,-t0)
        tend = min(tend, -t0+(NSTEP-1)*DT)
        !if (trim(station_name(irec))=='CC53') then
        chan_dat = trim(CH_CODE)//trim(RCOMPS(icomp))
        net = trim(network_name(irec))
        sta = trim(station_name(irec))
        call measure_adj_fwat(dat_inp_bp,syn_inp_bp,tstart,tend,dble(-t0),dble(DT),NSTEP,net,sta,&
                            chan_dat,window_chi,tr_chi,am_chi,&
                            T_pmax_dat,T_pmax_syn,adj_syn_all,file_prefix0,out_imeas,bandname)
        glob_file_prefix0(irec)=trim(file_prefix0)
        glob_net(irec)=trim(net)
        glob_sta(irec)=trim(sta)
        glob_chan_dat(irec,icomp)=trim(CH_CODE)//trim(RCOMPS(icomp)) !trim(chan_dat)
        glob_tstart(irec,icomp,1)=tstart
        glob_tend(irec,icomp,1)=tend
        glob_window_chi(irec,icomp,1:NCHI,1)=window_chi(1:NCHI)
        if (out_imeas==2) then
          glob_tr_chi(irec,icomp,1)=tr_chi/avgamp/avgamp *DT
          glob_am_chi(irec,icomp,1)=am_chi/avgamp/avgamp *DT
        else
          glob_tr_chi(irec,icomp,1)=tr_chi
          glob_am_chi(irec,icomp,1)=am_chi
        endif
          glob_T_pmax_dat(irec,icomp,1)=T_pmax_dat
          glob_T_pmax_syn(irec,icomp,1)=T_pmax_syn
          glob_imeas(irec,icomp,1)=out_imeas
        if (out_imeas==2) then
          misfit=glob_window_chi(irec,icomp,15,1)/avgamp/avgamp *DT
        else
          misfit=glob_window_chi(irec,icomp,15,1)
        endif
        total_misfit=total_misfit + misfit
        write(*,*)'stnm,chan,misfit:',trim(glob_sta(irec)),trim(glob_chan_dat(irec,icomp)),misfit
        !endif
      endif
    enddo 
  end subroutine meas_adj_tele

  subroutine meas_adj_tele_cd(ievt, irec, bandname,baz_all,ttp,win_tb,win_te,glob_dat_tw,&
                            glob_syn_tw, glob_ff, avgamp, window_chi,total_misfit, &
                            glob_net, glob_sta,glob_chan_dat, glob_tstart, &
                            glob_tend, glob_window_chi, glob_tr_chi, glob_am_chi, glob_num_win)
    implicit none
    ! type(sachead)                                            :: head
    ! logical                                                  :: findfile
    character(len=MAX_STRING_LEN)                            :: adjfile,bandname
    character(len=10)                                        :: net,sta,chan_dat
    integer                                                  :: irec, ievt, nc, nmax, ier
    real(kind=CUSTOM_REAL)                                   :: win_tb, win_te,total_misfit,&
                                                                misfit,avgamp
    character(len=256),dimension(nrec)                       :: glob_net,glob_sta
    character(len=256),dimension(nrec,NRCOMP)                :: glob_chan_dat   
    integer           ,dimension(nrec,NRCOMP)                :: glob_num_win   
    double precision, dimension(nrec,NRCOMP,NCHI,1)          :: glob_window_chi
    double precision, dimension(NCHI)                        :: window_chi
    double precision, dimension(nrec,NRCOMP,1)               :: glob_tr_chi,glob_am_chi,&
                                                                glob_tstart,glob_tend
    real(kind=4), dimension(nrec,NSTEP,NRCOMP)               :: glob_dat_tw,glob_syn_tw
    real(kind=4), dimension(nrec,NSTEP)                      :: glob_ff
    real(kind=4), dimension(nrec)                            :: ttp
    double precision                                         :: tstart,tend, app_half_dura, t00
    double precision, dimension(nrec)                        :: baz_all
    ! real(kind=CUSTOM_REAL), dimension(3,NSTEP)               :: seismo_syn
    real(kind=CUSTOM_REAL), dimension(:), allocatable        :: conv1,conv2,conv_full
    double precision, dimension(MAX_NDIM)                    :: datr_inp_bp,synr_inp_bp,&
                                                                datz_inp_bp,synz_inp_bp
    type(sachead)                                            :: head

    ! setup data and syn
    datr_inp_bp=0.
    synr_inp_bp=0.
    datz_inp_bp=0.
    synz_inp_bp=0.
    datz_inp_bp(1:NSTEP)=dble(glob_dat_tw(irec,1:NSTEP,1))/avgamp
    synz_inp_bp(1:NSTEP)=dble(glob_syn_tw(irec,1:NSTEP,1))
    datr_inp_bp(1:NSTEP)=dble(glob_dat_tw(irec,1:NSTEP,2))/avgamp
    synr_inp_bp(1:NSTEP)=dble(glob_syn_tw(irec,1:NSTEP,2))
    ! ----- get time shift between datz and synz -------
    ! app_half_dura = (abs(maxloc(abs(conv_full), dim=1) - NSTEP) + 1)* dt
    ! print *, 'net,sta,app_half_dura:', trim(network_name(irec)), trim(station_name(irec)), app_half_dura

    block
      integer, dimension(:), allocatable :: max_idx
      integer :: nttp, ndura, i, nb
      real(kind=cr), dimension(:), allocatable :: cut_conv, tmp_dat
      real(kind=cr) :: max_amp

      tmp_dat = real(datz_inp_bp(1:NSTEP))
      nttp = (ttp(irec)+T0)/DT+1
      ndura = (ttp(irec)+T0+max_duration)/DT+1
      cut_conv = abs(tmp_dat(nttp:ndura))
      max_idx = find_maxima(cut_conv)
      max_amp = maxval(cut_conv)
      do i = 1, size(max_idx)
        if (cut_conv(max_idx(i))>0.4*max_amp) then
          app_half_dura = max_idx(i)*DT
          exit
        end if
      end do
      print *, 'net,sta,app_half_dura:', trim(network_name(irec)),'_',trim(station_name(irec)), app_half_dura
    end block

    if (VERBOSE_MODE) then
      nmax = maxloc(synz_inp_bp, dim=1)
      t00 = nmax * DT + app_half_dura
      call sacio_newhead(head, sngl(DT), NSTEP, sngl(-t00))
      head%knetwk = trim(network_name(irec))
      head%kstnm = trim(station_name(irec))
      head%kcmpnm = trim(CH_CODE)//'Z'
      call myconvolution(real(datr_inp_bp(1:NSTEP)),real(synz_inp_bp(1:NSTEP)),&
                          NSTEP,NSTEP,conv1,1)
      adjfile=trim(OUTPUT_FILES)//'/wconv.'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)&
              //'Z.sac'//'.'//trim(bandname)
      call sacio_writesac(adjfile, head, dble(conv1(nmax:nmax+NSTEP-1)*DT), ier)      
      ! call dwsac1(trim(adjfile),dble(conv1(nmax:nmax+NSTEP-1)*DT),NSTEP,-t00,dble(DT))
      head%kcmpnm = trim(CH_CODE)//'R'
      call myconvolution(real(datz_inp_bp(1:NSTEP)),real(synr_inp_bp(1:NSTEP)),&
                          NSTEP,NSTEP,conv2,1)
      adjfile=trim(OUTPUT_FILES)//'/wconv.'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)&
              //'R.sac'//'.'//trim(bandname)
      call sacio_writesac(adjfile, head, dble(conv2(nmax:nmax+NSTEP-1)*DT), ier)
      ! call dwsac1(trim(adjfile),dble(conv2(nmax:nmax+NSTEP-1)*DT),NSTEP,-t00,dble(DT))
      deallocate(conv1,conv2)
    endif

    glob_num_win(irec,1)=1
    glob_num_win(irec,2)=0
    tstart=ttp(irec)-win_tb
    tend=ttp(irec)+win_te
    tstart = max(tstart,-T0)
    tend = min(tend, -T0+(NSTEP-1)*DT)
    chan_dat = trim(CH_CODE)
    net = trim(network_name(irec))
    sta = trim(station_name(irec))
    call meas_adj_conv_diff(datr_inp_bp, datz_inp_bp, synr_inp_bp, synz_inp_bp,&
                            tstart, tend, app_half_dura, NSTEP, net,sta,&
                            chan_dat, bandname, window_chi)
    glob_net(irec)=trim(net) !trim(net)
    glob_sta(irec)=trim(sta) !trim(sta)
    glob_chan_dat(irec,1)=trim(chan_dat)
    glob_tstart(irec,1,1)=tstart
    glob_tend(irec,1,1)=tend
    glob_window_chi(irec,1,1:NCHI,1)=window_chi(1:NCHI)
    misfit = real(window_chi(15))*src_weight(ievt) * DT
    glob_tr_chi(irec,1,1)=misfit
    glob_am_chi(irec,1,1)=misfit
    total_misfit = total_misfit+misfit
    write(*,*)'stnm,chan,misfit:',trim(net)//'.'//trim(sta), ',cross_conv_diff', real(misfit)
  end subroutine meas_adj_tele_cd

  subroutine meas_adj_rf(ievt, irec, igaus, bandname, glob_dat_tw, glob_syn_tw, glob_syn,&
                       win_tb, win_te, avgamp, ttp, window_chi, total_misfit, glob_net, glob_sta,&
                       glob_chan_dat, glob_tstart, glob_tend, glob_window_chi, glob_tr_chi, &
                       glob_am_chi, glob_num_win)
    
    implicit none
    integer                                          :: irec, ievt, igaus
    double precision, dimension(NSTEP)               :: dat_inp_bp, syn_inp_bp, adj_r_tw, adj_z_tw
    real(kind=4), dimension(nrec,NSTEP,NRCOMP)       :: glob_dat_tw,glob_syn_tw,glob_syn
    character(len=256),dimension(nrec)               :: glob_net, glob_sta
    character(len=256),dimension(nrec,NRCOMP)        :: glob_chan_dat
    double precision, dimension(NSTEP)               :: synr_bp, synz_bp
    double precision, dimension(nrec,NRCOMP,NCHI,1)  :: glob_window_chi
    double precision, dimension(nrec,NRCOMP,1)       :: glob_tstart, glob_tend
    double precision, dimension(nrec,NRCOMP,1)       :: glob_tr_chi
    double precision, dimension(nrec,NRCOMP,1)       :: glob_am_chi
    double precision, dimension(NCHI)                :: window_chi
    integer, dimension(nrec,NRCOMP)                  :: glob_num_win
    real(kind=4), dimension(nrec)                    :: ttp
    real                                             :: avgamp, win_tb, win_te
    character(len=MAX_STRING_LEN)                    :: bandname
    character(len=10)                                :: net,sta,chan_dat
    real(kind=CUSTOM_REAL)                           :: total_misfit, misfit

    synz_bp = 0.
    synr_bp = 0.
    synz_bp(:) = dble(glob_syn(irec,1:NSTEP,1))
    synr_bp(:) = dble(glob_syn(irec,1:NSTEP,2))
    dat_inp_bp = 0.
    syn_inp_bp = 0.
    glob_num_win(irec, 2) = 1
    dat_inp_bp = dble(glob_dat_tw(irec,1:NSTEP,2))
    syn_inp_bp = dble(glob_syn_tw(irec,1:NSTEP,2))
    net = trim(network_name(irec))
    sta = trim(station_name(irec))
    chan_dat = trim(CH_CODE)
    ! if (RF_ADJ_SRC_TYPE == 1) then
    !     call measure_adj_rf_data(dat_inp_bp,dble(win_tb),dble(win_te),dble(-t0),dble(ttp(irec)),&
    !                              dble(DT),NSTEP,f0(igaus),rf_tshift,net,sta,chan_dat,&
    !                              bandname,adj_r_tw, adj_z_tw)
    ! else
    call measure_adj_rf(dat_inp_bp,syn_inp_bp,synr_bp,synz_bp,dble(win_tb),dble(win_te),&
                dble(-t0),dble(ttp(irec)),dble(DT),NSTEP,rf_par%f0(igaus),&
                rf_par%rf_tshift,rf_par%maxit,rf_par%minderr,&
                net,sta,chan_dat, bandname, window_chi,adj_r_tw,adj_z_tw)
    glob_window_chi(irec,2,1:NCHI,1) = window_chi(1:NCHI)
    ! misfit = real(window_chi(15)/avgamp/avgamp *DT)
    misfit = real(window_chi(15))
    glob_window_chi(irec,2,3,1) = misfit* src_weight(ievt)*DT
    glob_tr_chi(irec,2,1)=misfit 
    glob_am_chi(irec,2,1)=misfit 
    total_misfit = total_misfit + misfit
    ! endif
    glob_net(irec)=trim(net)
    glob_sta(irec)=trim(sta)
    glob_chan_dat(irec,2)=trim(CH_CODE)//'R' !trim(chan_dat)
    glob_tstart(irec,2,1)=dble(ttp(irec)-win_tb)
    glob_tend(irec,2,1)=dble(ttp(irec)+win_te)
    write(*,*)'stnm,chan,misfit:',trim(net)//'.'//trim(sta),chan_dat, misfit
  end subroutine meas_adj_rf

  subroutine meas_adj_noise(ievt, irec, glob_sem_disp,iband, fstart0, fend0, bandname, &
                         baz_all, glob_dat_tw, glob_syn_tw, &
                         window_chi,tr_chi,am_chi,total_misfit,&
                         glob_file_prefix0, glob_net, glob_sta,  glob_chan_dat, glob_tstart, &
                         glob_tend, glob_window_chi, glob_tr_chi, glob_am_chi,  glob_T_pmax_dat, &
                         glob_T_pmax_syn, glob_imeas, glob_num_win)
                   
    !**********************  from run_preprocessing.f90 * *****************************
    integer                                                  :: irec, ierr
    integer                                                  :: ievt,icomp
  
    ! for reading dat in sac format
    type(sachead)                                            :: head
    logical                                                  :: findfile
    character(len=MAX_STRING_LEN)                            :: datafile,file_prefix0
    double precision, dimension(MAX_NDIM)                        :: datarray
    integer                                                  :: npt1,npt1_inp
    double precision                                         :: t01,t01_p,dt1
    ! for rotation
    real(kind=CUSTOM_REAL), dimension(3,NSTEP)               :: seismo_syn                       
    ! for preprocessing 
    double precision                                         :: t0_inp,t1_inp,dt_inp
    integer                                                  :: win_b, win_e ! for norm in signal window
    double precision, dimension(MAX_NDIM)                        :: dat_inp, syn_inp 
    double precision, dimension(MAX_NDIM)                        :: dat_inp_bp, syn_inp_bp
    double precision                                         :: fstart0,fend0
    ! for measurement 
    double precision                                         :: tstart,tend
    integer                                                  :: out_imeas
    integer                                                  :: NDIM_CUT
    double precision, dimension(MAX_NDIM)                        :: adj_syn_all
    character(len=MAX_STRING_LEN)                            :: bandname,adjfile
    integer                                                  :: yr,jda,ho,mi        
    integer                                                  :: iband
    double precision                                         :: sec,dist,az,baz,slat,slon
    double precision, dimension(nrec)                        :: baz_all
    character(len=10)                                        :: net,sta,chan_dat
    double precision, dimension(NCHI)                        :: window_chi
    double precision                                         :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn
    ! for collecting SEM synthetics to one large array
    real(kind=CUSTOM_REAL), dimension(nrec,NSTEP,3)          :: glob_sem_disp
    ! for calculating STF for TeleFWI
    real(kind=4), dimension(nrec,NSTEP,NRCOMP)               :: glob_dat_tw,glob_syn_tw
    ! for writing window_chi of noise/tele/leq by one proc
    character(len=256),dimension(nrec)                       :: glob_file_prefix0
    character(len=256),dimension(nrec)                       :: glob_net 
    character(len=256),dimension(nrec)                       :: glob_sta  
    character(len=256),dimension(nrec,NRCOMP)                :: glob_chan_dat   
    integer           ,dimension(nrec,NRCOMP)                :: glob_num_win   
    double precision, dimension(nrec,NRCOMP,1)               :: glob_tstart
    double precision, dimension(nrec,NRCOMP,1)               :: glob_tend
    double precision, dimension(nrec,NRCOMP,NCHI,1)          :: glob_window_chi
    double precision, dimension(nrec,NRCOMP,1)               :: glob_tr_chi
    double precision, dimension(nrec,NRCOMP,1)               :: glob_am_chi
    double precision, dimension(nrec,NRCOMP,1)               :: glob_T_pmax_dat
    double precision, dimension(nrec,NRCOMP,1)               :: glob_T_pmax_syn
    integer         , dimension(nrec,NRCOMP,1)               :: glob_imeas
    ! for writing STATIONS_***_ADJOINT
    real(kind=CUSTOM_REAL)                                   :: total_misfit,misfit

    seismo_syn(:,:)=0.d0 !glob_sem_disp(irec,:,:)
    do icomp=1,NRCOMP
      datafile='./'//trim(acqui_par%in_dat_path(ievt))//'/'//trim(network_name(irec))//'.'&
                 //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))//'.sac' 
      inquire(file=trim(datafile),exist=findfile)
      if ( findfile ) then
        call drsac1(trim(datafile),datarray,npt1,t01,dt1)
        call get_sacfile_header(trim(datafile),yr,jda,ho,mi,sec,net,sta, &
                                chan_dat,dist,az,baz,slat,slon)
        baz_all(irec)=baz
        !write(*,*)'myrank,datafile,npt1,t01,dt1,NSTEP,t0,DT,dist= ',myrank,trim(datafile),npt1,t01,dt1,&
        !                                   NSTEP,t0,DT,dist
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
        !**************** 1. pre-process dat and syn ****************************
        ! data are processed with procedures in the order of those shown
        ! in seis_process/process_data.pl
        dt_inp=dble(target_dt_noise) !!!! Here, I set dt=0.01 to interpolate data and syn
        t0_inp=-10.
        t1_inp=-t0+(NSTEP-1)*DT
        t01_p=t01+dt_inp ! interp from the second point to avoid large value of the first point ofor dif1
        npt1=floor(NSTEP*DT/dt1) ! In case of long data (large npt1), we choose part of the data
        npt1_inp=(npt1-1)*dt1/dt_inp
        NDIM_CUT=(t1_inp-t0_inp)/dt_inp+1
        dat_inp(:)=0.
        syn_inp(:)=0.
        dat_inp(1:npt1)=datarray(1:npt1)
        syn_inp(1:NSTEP)=dble(seismo_syn(icomp,1:NSTEP))
        ! interp
        if (npt1_inp>MAX_NDIM) call exit_mpi(myrank,'NDIM too small!!!') 
        call interpolate_syn(dat_inp,t01,dt1,npt1,t01_p,dt_inp,npt1_inp)
        call interpolate_syn(syn_inp,-t0,DT,NSTEP,t0_inp,dt_inp,NDIM_CUT)

        ! time derivative: CCFs==>EGFs
        if (.not. noise_par%SUPPRESS_EGF) then
          call dif1(dat_inp,dt_inp,NDIM_CUT)
          dat_inp=-1.0*dat_inp
          write(*,*) 'CCFs==>EGFs...'
        endif
        ! Cut
        call interpolate_syn(dat_inp,t01_p,dt_inp,npt1_inp,t0_inp,dt_inp,NDIM_CUT)

        !!! Filter
        ! rtrend
        !call detrend(dat_inp,NDIM_CUT,t0_inp,dt_inp)
        !call detrend(syn_inp,NDIM_CUT,t0_inp,dt_inp)
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
        if(.not. noise_par%USE_NEAR_OFFSET)then!BinHe added
          tstart=(noise_par%GROUPVEL_MAX(iband)+noise_par%GROUPVEL_MIN(iband))*LONG_P(iband)/2.0
          if(dist<tstart)then
            dat_inp_bp(:)=1.0e-10
            syn_inp_bp(:)=1.0e-10
          endif
        endif
        ! norm by maximum amplitude of signal window
        !!! ATTENTION: WK it is very important to add window for norm
        win_b=floor((dist/noise_par%GROUPVEL_MAX(iband)-LONG_P(iband)/2.+t0_inp)/dt_inp)!floor((dist/5+t0_inp)/dt_inp)
        win_e=floor((dist/noise_par%GROUPVEL_MIN(iband)+LONG_P(iband)/2.-t0_inp)/dt_inp)!floor((dist/2-t0_inp)/dt_inp)
        win_b=max(win_b,1)
        win_e=min(win_e,NDIM_CUT)
        dat_inp_bp=dat_inp_bp/maxval(abs(dat_inp_bp(win_b:win_e)))*maxval(abs(syn_inp_bp(win_b:win_e)))
        !!**************** 2. window ******************************************
        glob_num_win(irec,icomp)=1
        tstart=dist/noise_par%GROUPVEL_MAX(iband)-LONG_P(iband)/2.!BinHe modified advised by Kai
        tend=dist/noise_par%GROUPVEL_MIN(iband)+LONG_P(iband)/2.  !BinHe modified advised by Kai
        tstart = max(tstart,t0_inp)
        tend = min(tend, -t0+(NSTEP-2)*DT)

        if (VERBOSE_MODE) then
          call sacio_newhead(head, sngl(dt_inp), NDIM_CUT, sngl(t0_inp))
          head%knetwk = trim(net)
          head%kstnm = trim(sta)
          head%kcmpnm = trim(chan_dat)
          head%az = az
          head%baz = baz
          head%dist = dist
          head%stla = slat
          head%stlo = slon
          head%t1 = tstart
          head%t2 = tend
          adjfile=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
              //'.syn.sac'//'.'//trim(bandname)
          call sacio_writesac(trim(adjfile), head, syn_inp_bp(1:NDIM_CUT), ierr)
          adjfile=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
              //'.obs.sac'//'.'//trim(bandname) 
          call sacio_writesac(trim(adjfile), head, dat_inp_bp(1:NDIM_CUT), ierr)
        endif
        !!**************** 3. measure_adj ******************************************
        !if (trim(station_name(irec))=='CC53') then
        chan_dat = trim(CH_CODE)//trim(RCOMPS(icomp))
        net = trim(network_name(irec))
        sta = trim(station_name(irec))
        call measure_adj_fwat(dat_inp_bp,syn_inp_bp,tstart,tend,t0_inp,dt_inp,NDIM_CUT,net,sta,&
                              chan_dat,window_chi,tr_chi,am_chi,&
                              T_pmax_dat,T_pmax_syn,adj_syn_all,file_prefix0,out_imeas,bandname)
        glob_file_prefix0(irec)=trim(file_prefix0)
        glob_net(irec)=trim(net) !trim(net)
        glob_sta(irec)=trim(sta) !trim(sta)
        glob_chan_dat(irec,icomp)=trim(CH_CODE)//trim(RCOMPS(icomp)) !trim(chan_dat)
        glob_tstart(irec,icomp,1)=tstart
        glob_tend(irec,icomp,1)=tend
        glob_window_chi(irec,icomp,1:NCHI,1)=window_chi(1:NCHI)
        glob_tr_chi(irec,icomp,1)=tr_chi * src_weight(ievt)
        glob_am_chi(irec,icomp,1)=am_chi * src_weight(ievt)
        glob_T_pmax_dat(irec,icomp,1)=T_pmax_dat
        glob_T_pmax_syn(irec,icomp,1)=T_pmax_syn
        glob_imeas(irec,icomp,1)=out_imeas

        misfit=tr_chi
        total_misfit=total_misfit + misfit
        !endif
      endif ! end findfile
    enddo ! end icomp
  end subroutine

  subroutine meas_adj_leq(icmt, ievt, irec, glob_sem_disp,iband, fstart0, fend0, bandname, &
                        baz_all, glob_dat_tw, glob_syn_tw, &
                        window_chi,tr_chi,am_chi,total_misfit,&
                        glob_file_prefix0, glob_net, glob_sta,  glob_chan_dat, glob_tstart, &
                        glob_tend, glob_window_chi, glob_tr_chi, glob_am_chi,  glob_T_pmax_dat, &
                        glob_T_pmax_syn, glob_imeas, glob_num_win)

    !**********************  from run_preprocessing.f90 * *****************************
    integer                                                  :: irec
    integer                                                  :: ievt,icomp,iwin,icmt,ier
  
    ! for reading dat in sac format
    logical                                                  :: findfile
    character(len=MAX_STRING_LEN)                            :: datafile,file_prefix0
    double precision, dimension(MAX_NDIM)                        :: datarray
    integer                                                  :: npt1
    double precision                                         :: t01,dt1
    ! sacio from sacio.f90
    type(sachead) :: head
    ! for rotation
    real(kind=CUSTOM_REAL), dimension(3,NSTEP)               :: seismo_syn                       
    ! for preprocessing 
    double precision                                         :: t0_inp,t1_inp,dt_inp
    double precision, dimension(MAX_NDIM)                        :: dat_inp, syn_inp 
    double precision, dimension(MAX_NDIM)                        :: dat_inp_bp, syn_inp_bp
    double precision                                         :: fstart0,fend0
    ! for measurement 
    double precision                                         :: tstart,tend
    integer                                                  :: out_imeas
    integer                                                  :: NDIM_CUT
    double precision, dimension(MAX_NDIM)                        :: adj_syn_all
    character(len=MAX_STRING_LEN)                            :: bandname
    integer                                                  :: yr,jda,ho,mi        
    integer                                                  :: iband
    double precision                                         :: sec,dist,az,baz,slat,slon
    double precision, dimension(nrec)                        :: baz_all
    character(len=10)                                        :: net,sta,chan_dat
    double precision, dimension(NCHI)                        :: window_chi
    double precision                                         :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn
    ! for collecting SEM synthetics to one large array
    real(kind=CUSTOM_REAL), dimension(nrec,NSTEP,3)          :: glob_sem_disp
    ! for calculating STF for TeleFWI
    real(kind=4), dimension(nrec,NSTEP,NRCOMP)               :: glob_dat_tw,glob_syn_tw
    ! for writing window_chi of noise/tele/leq by one proc
    character(len=256),dimension(nrec)                       :: glob_file_prefix0
    character(len=256),dimension(nrec)                       :: glob_net 
    character(len=256),dimension(nrec)                       :: glob_sta  
    character(len=256),dimension(nrec,NRCOMP)                :: glob_chan_dat   
    integer           ,dimension(nrec,NRCOMP)                :: glob_num_win   
    double precision, dimension(nrec,NRCOMP,1)               :: glob_tstart
    double precision, dimension(nrec,NRCOMP,1)               :: glob_tend
    double precision, dimension(nrec,NRCOMP,NCHI,1)          :: glob_window_chi
    double precision, dimension(nrec,NRCOMP,1)               :: glob_tr_chi
    double precision, dimension(nrec,NRCOMP,1)               :: glob_am_chi
    double precision, dimension(nrec,NRCOMP,1)               :: glob_T_pmax_dat
    double precision, dimension(nrec,NRCOMP,1)               :: glob_T_pmax_syn
    integer         , dimension(nrec,NRCOMP,1)               :: glob_imeas
    ! for writing STATIONS_***_ADJOINT
    character(len=MAX_STRING_LEN)                            :: dummystring
    real(kind=CUSTOM_REAL)                                   :: total_misfit,misfit

    seismo_syn(:,:)=0.d0 !glob_sem_disp(irec,:,:)
    call exit_mpi(myrank,'****** FWAT of LEQ is still under development by Kai Wang. Please do NOT use this !!! ******') 
    do icomp=1,NRCOMP
      datafile='./'//trim(acqui_par%in_dat_path(ievt))//'/'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))//'.sac' 
      inquire(file=trim(datafile),exist=findfile)
      if ( findfile ) then
        call read_parameter_file_flexwin()
        call drsac1(trim(datafile),datarray,npt1,t01,dt1)
        call get_sacfile_header(trim(datafile),yr,jda,ho,mi,sec,net,sta, &
                                chan_dat,dist,az,baz,slat,slon)
        !------------------------------
        ! The following sac headers are required by ttimes in FLEXWIN            
        b=t01
        npts=npt1

        call sacio_readhead(trim(datafile), head, ier)
        if(ier /= 0) then
          write(*,*)'Error reading SAC header: ', trim(datafile)
          call exit(-1)
        endif
        evla = head%evla
        evlo = head%evlo
        stla = head%stla
        stlo = head%stlo
        evdp = head%evdp
        kstnm = head%kstnm
        kcmpnm = head%kcmpnm
        knetwk = head%knetwk

        if (BODY_WAVE_ONLY) then
          P_pick = head%t1
          S_pick = head%t2
        endif

        ! calculate distances and azimuths
        call distaz(stla,stlo,evla,evlo,azimuth,backazimuth,dist_deg,dist_km)

        ! Frequency limits may be conditional on station or event information
        ! so call user function to modify them if required
        WIN_MIN_PERIOD=SHORT_P(iband)
        WIN_MAX_PERIOD=LONG_P(iband)
        call modify_T0_T1_on_condition
        !------------------------

        baz_all(irec)=baz
        !write(*,*)'myrank,datafile,npt1,t01,dt1,NSTEP,t0,DT,dist= ',myrank,trim(datafile),npt1,t01,dt1,&
        !                                   NSTEP,t0,DT,dist
        if(abs(dt1-dT)>0.0001) stop 'delta of data does NOT match DT of SEM synthetics'
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
        !**************** 1. pre-process dat and syn ****************************
        ! data are processed with procedures in the order of those shown
        ! in seis_process/process_data.pl
        dt_inp=DT 
        t0_inp=-t0
        t1_inp=-t0+(NSTEP-1)*DT
        NDIM_CUT=(t1_inp-t0_inp)/dt_inp+1
        dat_inp(:)=0.
        syn_inp(:)=0.
        dat_inp(1:npt1)=datarray(1:npt1)
        syn_inp(1:NSTEP)=dble(seismo_syn(icomp,1:NSTEP))
        !!! Filter
        ! rtrend
        !call detrend(dat_inp,NDIM_CUT,t0_inp,dt_inp)
        !call detrend(syn_inp,NDIM_CUT,t0_inp,dt_inp)
        call detrend(dat_inp,NDIM_CUT)
        call detrend(syn_inp,NDIM_CUT)
        ! remean
        dat_inp(1:NDIM_CUT)=dat_inp(1:NDIM_CUT)-sum(dat_inp(1:NDIM_CUT))/NDIM_CUT
        syn_inp(1:NDIM_CUT)=syn_inp(1:NDIM_CUT)-sum(syn_inp(1:NDIM_CUT))/NDIM_CUT
        dat_inp_bp=dat_inp
        syn_inp_bp=syn_inp
        ! band pass filter
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
          dummystring=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
          //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
          //'.syn.sac'//'.'//trim(bandname) 
          call dwsac0(trim(dummystring),syn_inp_bp,NDIM_CUT,t0_inp,dt_inp,knetwk,kstnm, &
                      kcmpnm,dist,az,baz,dble(stla),dble(stlo),dble(evla),dble(evlo),dble(evdp))

          dummystring=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
          //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
          //'.obs.sac'//'.'//trim(bandname) 
          call dwsac0(trim(dummystring),dat_inp_bp,NDIM_CUT,t0_inp,dt_inp,knetwk,kstnm, &
                      kcmpnm,dist,az,baz,dble(stla),dble(stlo),dble(evla),dble(evlo),dble(evdp))
        endif
        if (leq_par%CMT3D_INV) then  
          if (icmt>0) then
            dummystring=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
            //'.syn.sac'//'.'//trim(bandname)//'.'//trim(leq_par%CMT_PAR(icmt)) 
            call dwsac0(trim(dummystring),syn_inp_bp,NDIM_CUT,t0_inp,dt_inp,knetwk,kstnm, &
                      kcmpnm,dist,az,baz,dble(stla),dble(stlo),dble(evla),dble(evlo),dble(evdp))
          else
            dummystring=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
            //'.obs.sac'//'.'//trim(bandname) 
            call dwsac0(trim(dummystring),dat_inp_bp,NDIM_CUT,t0_inp,dt_inp,knetwk,kstnm, &
                      kcmpnm,dist,az,baz,dble(stla),dble(stlo),dble(evla),dble(evlo),dble(evdp))
            dummystring=trim(OUTPUT_FILES)//'/'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
            //'.syn.sac'//'.'//trim(bandname) 
            call dwsac0(trim(dummystring),syn_inp_bp,NDIM_CUT,t0_inp,dt_inp,knetwk,kstnm, &
                      kcmpnm,dist,az,baz,dble(stla),dble(stlo),dble(evla),dble(evlo),dble(evdp))
          endif 
        endif
        if (icmt==0) then
          !!**************** 2. flexwin ******************************************
          obs_lp(1:MAX_NDIM)=dat_inp_bp(1:MAX_NDIM)
          synt_lp(1:MAX_NDIM)=syn_inp_bp(1:MAX_NDIM)
          call select_windows_stalta2()
          glob_num_win(irec,icomp)=num_win
          if ( num_win>0 ) then
            do iwin=1,num_win
              tstart = win_start(iwin)
              tend = win_end(iwin)
              if ( .not. leq_par%CMT3D_INV ) then
              !!**************** 3. measure_adj ******************************************
              out_imeas=0
              call measure_adj_fwat(dat_inp_bp,syn_inp_bp,tstart,tend,t0_inp,dt_inp,NDIM_CUT,net,sta,&
                                chan_dat,window_chi,tr_chi,am_chi,&
                                T_pmax_dat,T_pmax_syn,adj_syn_all,file_prefix0,out_imeas,bandname)
              endif
              glob_file_prefix0(irec)=trim(file_prefix0)
              glob_net(irec)=trim(net)
              glob_sta(irec)=trim(sta)
              glob_chan_dat(irec,icomp)=trim(CH_CODE)//trim(RCOMPS(icomp)) !trim(chan_dat)
              glob_tstart(irec,icomp,iwin)=tstart
              glob_tend(irec,icomp,iwin)=tend
              glob_window_chi(irec,icomp,1:NCHI,iwin)=window_chi(1:NCHI)
              glob_tr_chi(irec,icomp,iwin)=tr_chi
              glob_am_chi(irec,icomp,iwin)=am_chi
              glob_T_pmax_dat(irec,icomp,iwin)=T_pmax_dat
              glob_T_pmax_syn(irec,icomp,iwin)=T_pmax_syn
              glob_imeas(irec,icomp,iwin)=out_imeas

              misfit=tr_chi
              total_misfit=total_misfit + misfit
            enddo
          endif ! nwin > 0
        endif ! icmt==0
      endif ! end findfile
    enddo ! end icomp  
  end subroutine

  subroutine sum_adj_source(ievt, irec, simu_type, adj_syn_all_sum, baz_all, stf_array, avgamp, num_adj)

    character(len=MAX_STRING_LEN)                    :: adjfile,bandname,datafile,simu_type
    integer                                          :: icomp, iband, npt1, num_adj, irec, ier, ievt, nflt
    double precision, dimension(3,MAX_NDIM)              :: adj_syn_all_sum
    real(kind=CUSTOM_REAL), dimension(3,MAX_NDIM)        :: adj_arrays_zne
    double precision, dimension(MAX_NDIM)                :: datarray
    double precision                                 :: baz_all(nrec)
    real                                             :: avgamp
    logical                                          :: findfile
    real(kind=4), dimension(:), allocatable          :: tmpl
    double precision                                 :: t01,dt1,adj_amp_max
    real(kind=4) , dimension(NSTEP,NRCOMP)           :: stf_array
    character(len=1), dimension(3)                   :: comp_name


    num_adj = 0
    adj_arrays_zne = 0.
    if (simu_type == 'rf') then
      nflt = rf_par%NGAUSS
    else
      nflt = NUM_FILTER
    endif
    do icomp=1,NRCOMP
      adjfile=trim(OUTPUT_FILES)//'/../SEM/'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
              //'.adj' 
      adj_syn_all_sum(icomp,:)=0.
      datarray=0.
      ! find maximum amplitude of adjoint source at diferent period bands
      adj_amp_max=0.0
      do iband=1,nflt
        if (simu_type == 'rf') then
          write(bandname,'(a1,f3.1)') 'F',rf_par%f0(iband)
        else
          call get_band_name(SHORT_P(iband),LONG_P(iband),bandname)
          ! write(bandname,'(a1,i3.3,a2,i3.3)') 'T',int(SHORT_P(iband)),'_T',int(LONG_P(iband))
        endif
        datafile=trim(OUTPUT_FILES)//'/../SEM/'//trim(network_name(irec))//'.'&
                //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
                //'.adj.sac'//'.'//trim(bandname) 
        inquire(file=trim(datafile),exist=findfile)
        if ( findfile ) then
          !write(*,*)trim(datafile)
          call drsac1(trim(datafile),datarray,npt1,t01,dt1)
          if (simu_type=="tele") then
            datarray=datarray/avgamp ! normalize by maximum abs of the data gather
            allocate(tmpl(npt1))
            call myconvolution(real(datarray),stf_array(:,icomp),npt1,npt1,tmpl,0)
            datarray(1:npt1)=dble(tmpl)*dt1
            deallocate(tmpl)
          endif
          if(maxval(abs(datarray(:))) > adj_amp_max ) adj_amp_max=maxval(abs(datarray(:)))
        endif ! findfile
      enddo ! iband=1,nflt
      !  
      do iband=1,nflt
        if (simu_type == 'rf') then
          write(bandname,'(a1,f3.1)') 'F',rf_par%f0(iband)
        else
          call get_band_name(SHORT_P(iband),LONG_P(iband),bandname)
          ! write(bandname,'(a1,i3.3,a2,i3.3)') 'T',int(SHORT_P(iband)),'_T',int(LONG_P(iband))
        endif
        datafile=trim(OUTPUT_FILES)//'/../SEM/'//trim(network_name(irec))//'.'&
                //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(RCOMPS(icomp))&
                //'.adj.sac'//'.'//trim(bandname) 
        inquire(file=trim(datafile),exist=findfile)
        if ( findfile ) then
          !write(*,*)trim(datafile)
          call drsac1(trim(datafile),datarray,npt1,t01,dt1)
          if (.not. VERBOSE_MODE) then
            open(unit=1234, iostat=ier, file=trim(datafile), status='old')
            if (ier == 0) close(1234, status='delete') 
          endif
          if (simu_type=="tele") then
            datarray=datarray/avgamp ! normalize by maximum abs of the data gather
            allocate(tmpl(npt1))
            call mycorrelation(real(datarray),stf_array(:,icomp),npt1,npt1,tmpl,0)
            datarray(1:npt1)=dble(tmpl)*dt1
            if (maxval(abs(datarray)).ne.0) then
              datarray=datarray/maxval(abs(stf_array(:,icomp)))
            endif
            deallocate(tmpl)
          endif
          if(noise_par%ADJ_SRC_NORM.and.simu_type=="noise".and.maxval(abs(datarray)).ne.0.) then
            adj_syn_all_sum(icomp,:)=adj_syn_all_sum(icomp,:)+datarray(:)/maxval(abs(datarray)) &
                                      *adj_amp_max
          else
            adj_syn_all_sum(icomp,:)=adj_syn_all_sum(icomp,:)+datarray(:)
          endif
          num_adj=num_adj+1
        endif ! findfile
      enddo ! iband=1,nflt
    enddo ! icomp=1,NRCOMP
    adj_arrays_zne=0.
    adj_syn_all_sum = adj_syn_all_sum * src_weight(ievt)
    if (trim(dat_coord)=='ZRT') then
      call rotate_ZRT_to_ZNE(real(adj_syn_all_sum(1,:)),real(adj_syn_all_sum(2,:)),real(adj_syn_all_sum(3,:)), &
              adj_arrays_zne(1,:),adj_arrays_zne(2,:),adj_arrays_zne(3,:),MAX_NDIM,real(baz_all(irec)))      
    else
      adj_arrays_zne=adj_syn_all_sum
    endif
    if (num_adj>0) then
      if (SUPPRESS_UTM_PROJECTION) then
        comp_name(1) = 'Z'
        comp_name(2) = 'Y'
        comp_name(3) = 'X'
      else
        comp_name(1) = 'Z'
        comp_name(2) = 'N'
        comp_name(3) = 'E'
      endif
      adjfile=trim(OUTPUT_FILES)//'/../SEM/'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(comp_name(1))//'.adj' 
      call dwascii(trim(adjfile),dble(adj_arrays_zne(1,:)),npt1,t01,dt1)
      adjfile=trim(OUTPUT_FILES)//'/../SEM/'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(comp_name(2))//'.adj' 
      call dwascii(trim(adjfile),dble(adj_arrays_zne(2,:)),npt1,t01,dt1)
      adjfile=trim(OUTPUT_FILES)//'/../SEM/'//trim(network_name(irec))//'.'&
              //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(comp_name(3))//'.adj' 
      call dwascii(trim(adjfile),dble(adj_arrays_zne(3,:)),npt1,t01,dt1)
      if (VERBOSE_MODE) then
        adjfile=trim(OUTPUT_FILES)//'/../SEM/'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(comp_name(1))//'.adj.sac' 
        call dwsac1(trim(adjfile),dble(adj_arrays_zne(1,:)),npt1,t01,dt1) 
        adjfile=trim(OUTPUT_FILES)//'/../SEM/'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(comp_name(2))//'.adj.sac' 
        call dwsac1(trim(adjfile),dble(adj_arrays_zne(2,:)),npt1,t01,dt1) 
            adjfile=trim(OUTPUT_FILES)//'/../SEM/'//trim(network_name(irec))//'.'&
            //trim(station_name(irec))//'.'//trim(CH_CODE)//trim(comp_name(3))//'.adj.sac' 
        call dwsac1(trim(adjfile),dble(adj_arrays_zne(3,:)),npt1,t01,dt1) 
      endif           
    endif
  end subroutine sum_adj_source

end module preproc_measure_adj_subs

