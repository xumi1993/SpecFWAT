module collect_data

use specfem_par, only: seismograms_d,number_receiver_global, myrank, nrec, nrec_local, &
                       CUSTOM_REAL, islice_selected_rec, MAX_STRING_LEN
use shared_input_parameters, only: NSTEP, DT, NPROC, SUPPRESS_UTM_PROJECTION
use fwat_input 
use fullwave_adjoint_tomo_par 
use my_mpi

implicit none

integer                                        :: irank, irec_local, tag, nsta_irank, &
                                                  icomp, irec, ier
integer, dimension(MPI_STATUS_SIZE)            :: status

contains

subroutine collect_seismograms_d(glob_sem_disp)
  integer, dimension(:), allocatable       :: irec_global 
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: sem_disp_loc
  real(kind=CUSTOM_REAL), dimension(nrec,NSTEP,3)  :: glob_sem_disp

  do irank=1,NPROC-1
    if (myrank == 0) then
      nsta_irank=0
      do irec = 1,  nrec
        if (islice_selected_rec(irec) == irank) nsta_irank = nsta_irank + 1
      enddo
      if (nsta_irank > 0) then
        allocate(irec_global(nsta_irank),stat=ier)
        allocate(sem_disp_loc(nsta_irank,NSTEP,3))
        !! data to receive
        irec_local=0
        tag   = 2001 !MPI_ANY_TAG

        call MPI_RECV(sem_disp_loc, nsta_irank*NSTEP*3, MPI_REAL, irank, tag, my_local_mpi_comm_world, status,  ier)
        call MPI_RECV(irec_global, nsta_irank, MPI_INTEGER, irank, tag, my_local_mpi_comm_world, status,  ier)
        do icomp=1,3
          do irec_local = 1, nsta_irank
            glob_sem_disp(irec_global(irec_local),:,icomp) = sem_disp_loc(irec_local,:, icomp)
          enddo
        enddo
        deallocate(sem_disp_loc)
        deallocate(irec_global)
      endif
    else
      if (myrank == irank .and. nrec_local > 0) then
        !write(*,*)'myrank,irank,nrec_local=',myrank,irank,nrec_local,number_receiver_global
        allocate(irec_global(nrec_local),stat=ier)
        allocate(sem_disp_loc(nrec_local,NSTEP,3))

        do irec_local = 1, nrec_local
          irec_global(irec_local) = number_receiver_global(irec_local)
          !! choose the rigth seismograms_*
          do icomp=1,3
              sem_disp_loc(irec_local,:,icomp)=seismograms_d(icomp,irec_local,:)
          enddo
        enddo

        tag    = 2001
        call MPI_SEND(sem_disp_loc,nrec_local*NSTEP*3,MPI_REAL,0,tag,my_local_mpi_comm_world,ier)
        call MPI_SEND(irec_global,nrec_local,MPI_INTEGER,0,tag,my_local_mpi_comm_world,ier)
        deallocate(sem_disp_loc)
        deallocate(irec_global)

      endif
    endif
    ! not sure if need this sync
    call synchronize_all()
  enddo ! end loop irank
    !!! WK added: now collect data on rank 0
  if (myrank==0) then
    allocate(irec_global(nrec_local),stat=ier)
    do irec_local = 1, nrec_local
      irec_global(irec_local) = number_receiver_global(irec_local)
      !! choose the rigth seismograms_*
      do icomp=1,3
        glob_sem_disp(irec_global(irec_local),:,icomp)=seismograms_d(icomp,irec_local,:)
      enddo
    enddo
    deallocate(irec_global)
  endif
  call synchronize_all()
end subroutine collect_seismograms_d

subroutine collect_stf_deconvolution(glob_stnm, glob_dat_tw, glob_syn_tw, glob_ff, glob_syn, nrec_local_max, my_nrec_local)
  character(len=MAX_STRING_LEN), dimension(nrec)  :: glob_stnm
  character(len=MAX_STRING_LEN), dimension(:), allocatable :: stnm_loc
  real(kind=4), dimension(nrec,NSTEP,NRCOMP)  :: glob_dat_tw,glob_syn_tw, glob_syn
  real(kind=4), dimension(:,:,:), allocatable :: dat_tw_loc,syn_tw_loc, syn_loc
  real(kind=4), dimension(nrec,NSTEP)  :: glob_ff
  real(kind=4), dimension(:,:), allocatable :: ff_loc
  integer :: nrec_local_max, my_nrec_local

  do irank=1,NPROC-1
    if (myrank == 0) then
      nsta_irank=0
      !do irec = 1,  nrec
      !   if (islice_selected_rec(irec) == irank) nsta_irank = nsta_irank + 1
      !enddo
      do irec_local = 1, nrec_local_max
          irec = irank*nrec_local_max+irec_local
          if(irec<=nrec) nsta_irank=nsta_irank+1
      enddo

      !write(*,*)'myrank,irank,nsta_irank=',myrank,irank,nsta_irank
      if (nsta_irank > 0) then
        allocate(stnm_loc(nsta_irank))
        allocate(dat_tw_loc(nsta_irank,NSTEP,NRCOMP))
        allocate(syn_tw_loc(nsta_irank,NSTEP,NRCOMP))
        allocate(syn_loc(nsta_irank,NSTEP,NRCOMP))
        allocate(ff_loc(nsta_irank,NSTEP))
        !! data to receive
        irec_local=0
        tag   = 2001 !MPI_ANY_TAG

        call MPI_RECV(stnm_loc,  nsta_irank*256, MPI_CHARACTER, irank, tag, my_local_mpi_comm_world, status,  ier)
        call MPI_RECV(dat_tw_loc, nsta_irank*NSTEP*NRCOMP, MPI_REAL, irank, tag, my_local_mpi_comm_world, status,  ier)
        call MPI_RECV(syn_tw_loc, nsta_irank*NSTEP*NRCOMP, MPI_REAL, irank, tag, my_local_mpi_comm_world, status,  ier)
        call MPI_RECV(syn_loc, nsta_irank*NSTEP*NRCOMP, MPI_REAL, irank, tag, my_local_mpi_comm_world, status,  ier)
        call MPI_RECV(ff_loc, nsta_irank*NSTEP, MPI_REAL, irank, tag, my_local_mpi_comm_world, status,  ier)
        do icomp=1,NRCOMP
          do irec_local = 1, nsta_irank
            irec = irank*nrec_local_max+irec_local
            glob_stnm(irec) = trim(stnm_loc(irec_local))
            glob_dat_tw(irec,:,icomp) = dat_tw_loc(irec_local,:, icomp)
            glob_syn_tw(irec,:,icomp) = syn_tw_loc(irec_local,:, icomp)
            glob_syn(irec,:,icomp) = syn_loc(irec_local,:, icomp)
            glob_ff(irec,:) = ff_loc(irec_local,:)
          enddo
        enddo
        deallocate(stnm_loc)
        deallocate(dat_tw_loc)
        deallocate(syn_tw_loc)
        deallocate(syn_loc)
        deallocate(ff_loc)
      endif
    else
      if (myrank == irank .and. my_nrec_local > 0) then
        !write(*,*)'myrank,irank,nrec_local=',myrank,irank,my_nrec_local
        allocate(stnm_loc(my_nrec_local))
        allocate(dat_tw_loc(my_nrec_local,NSTEP,NRCOMP))
        allocate(syn_tw_loc(my_nrec_local,NSTEP,NRCOMP))
        allocate(syn_loc(my_nrec_local,NSTEP,NRCOMP))
        allocate(ff_loc(my_nrec_local,NSTEP))

        do irec_local = 1, my_nrec_local
          irec = irank*nrec_local_max+irec_local
          !! choose the rigth seismograms_*
          do icomp=1,NRCOMP
          stnm_loc(irec_local)=trim(glob_stnm(irec))
          dat_tw_loc(irec_local,:,icomp)=glob_dat_tw(irec,:, icomp)
          syn_tw_loc(irec_local,:,icomp)=glob_syn_tw(irec,:, icomp)
          syn_loc(irec_local,:,icomp)=glob_syn(irec,:, icomp)
          ff_loc(irec_local,:)=glob_ff(irec,:)
          enddo
        enddo

        tag = 2001
        call MPI_SEND(stnm_loc,  my_nrec_local*256, MPI_CHARACTER, 0, tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(dat_tw_loc,  my_nrec_local*NSTEP*NRCOMP, MPI_REAL, 0, tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(syn_tw_loc,  my_nrec_local*NSTEP*NRCOMP, MPI_REAL, 0, tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(syn_loc,  my_nrec_local*NSTEP*NRCOMP, MPI_REAL, 0, tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(ff_loc,  my_nrec_local*NSTEP, MPI_REAL, 0, tag, my_local_mpi_comm_world, ier)
        deallocate(stnm_loc)
        deallocate(dat_tw_loc)
        deallocate(syn_tw_loc)
        deallocate(syn_loc)
        deallocate(ff_loc)
      endif
    endif
    ! not sure if need this sync
    call synchronize_all()
  enddo ! end loop irank
end subroutine collect_stf_deconvolution

subroutine collect_chi(glob_file_prefix0, glob_net, glob_sta,  glob_chan_dat, glob_tstart, &
                       glob_tend, glob_window_chi, glob_tr_chi, glob_am_chi,  glob_T_pmax_dat, &
                       glob_T_pmax_syn, glob_imeas, glob_num_win, nrec_local_max, my_nrec_local)
  use measure_adj_mod

  character(len=256),dimension(nrec) :: glob_file_prefix0
  character(len=256),dimension(nrec) :: glob_net 
  character(len=256),dimension(nrec) :: glob_sta  
  character(len=256),dimension(nrec,NRCOMP) :: glob_chan_dat   
  integer         , dimension(nrec,NRCOMP)  :: glob_num_win   
  double precision, dimension(nrec,NRCOMP,leq_par%FLEX_NWIN)  :: glob_tstart
  double precision, dimension(nrec,NRCOMP,leq_par%FLEX_NWIN)  :: glob_tend
  double precision, dimension(nrec,NRCOMP,NCHI,leq_par%FLEX_NWIN)  :: glob_window_chi
  double precision, dimension(nrec,NRCOMP,leq_par%FLEX_NWIN)  :: glob_tr_chi
  double precision, dimension(nrec,NRCOMP,leq_par%FLEX_NWIN)  :: glob_am_chi
  double precision, dimension(nrec,NRCOMP,leq_par%FLEX_NWIN)  :: glob_T_pmax_dat
  double precision, dimension(nrec,NRCOMP,leq_par%FLEX_NWIN)  :: glob_T_pmax_syn
  integer         , dimension(nrec,NRCOMP,leq_par%FLEX_NWIN)  :: glob_imeas
  character(len=256),dimension(:), allocatable :: file_prefix0_loc
  character(len=256),dimension(:), allocatable :: net_loc 
  character(len=256),dimension(:), allocatable :: sta_loc  
  character(len=256),dimension(:,:), allocatable :: chan_dat_loc   
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: tstart_loc
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: tend_loc
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable  :: window_chi_loc
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: tr_chi_loc
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: am_chi_loc
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: T_pmax_dat_loc
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: T_pmax_syn_loc
  integer               , dimension(:,:,:), allocatable  :: imeas_loc 
  integer               , dimension(:,:), allocatable  :: num_win_loc
  integer :: nrec_local_max, my_nrec_local,iwin

  do irank=1,NPROC-1
    if (myrank == 0) then
      ! count the receiver in slice irank
      !write(*,*)'NPROC,myrank,irank=',NPROC,myrank,irank
      nsta_irank=0
      !do irec = 1,  nrec
      !   if (islice_selected_rec(irec) == irank) nsta_irank = nsta_irank + 1
      !enddo
      do irec_local=1,nrec_local_max
        irec=irank*nrec_local_max+irec_local
        if(irec<=nrec) nsta_irank=nsta_irank+1
      enddo
      if (nsta_irank > 0) then
        allocate(file_prefix0_loc(nsta_irank))
        allocate(net_loc(nsta_irank))
        allocate(sta_loc(nsta_irank))
        allocate(chan_dat_loc(nsta_irank,NRCOMP))
        allocate(num_win_loc(nsta_irank,NRCOMP))
        allocate(tstart_loc(nsta_irank,NRCOMP,leq_par%FLEX_NWIN))
        allocate(tend_loc(nsta_irank,NRCOMP,leq_par%FLEX_NWIN))
        allocate(window_chi_loc(nsta_irank,NRCOMP,NCHI,leq_par%FLEX_NWIN))
        allocate(tr_chi_loc(nsta_irank,NRCOMP,leq_par%FLEX_NWIN))
        allocate(am_chi_loc(nsta_irank,NRCOMP,leq_par%FLEX_NWIN))
        allocate(T_pmax_dat_loc(nsta_irank,NRCOMP,leq_par%FLEX_NWIN))
        allocate(T_pmax_syn_loc(nsta_irank,NRCOMP,leq_par%FLEX_NWIN))
        allocate(imeas_loc(nsta_irank,NRCOMP,leq_par%FLEX_NWIN))
        !! data to receive
        irec_local=0
        tag   = 2001 !MPI_ANY_TAG

        call MPI_RECV(file_prefix0_loc,  nsta_irank*256, MPI_CHARACTER, irank, tag, my_local_mpi_comm_world, status,  ier)
        call MPI_RECV(net_loc,  nsta_irank*256, MPI_CHARACTER, irank, tag, my_local_mpi_comm_world, status,ier)
        call MPI_RECV(sta_loc,  nsta_irank*256, MPI_CHARACTER, irank, tag, my_local_mpi_comm_world, status,ier)
        call MPI_RECV(chan_dat_loc,nsta_irank*256*NRCOMP,MPI_CHARACTER,irank,tag,my_local_mpi_comm_world,status,ier)
        call MPI_RECV(num_win_loc,nsta_irank*NRCOMP,MPI_INTEGER,irank,tag,my_local_mpi_comm_world,status,ier)
        call MPI_RECV(tstart_loc,nsta_irank*NRCOMP*leq_par%FLEX_NWIN,MPI_REAL,irank, tag,my_local_mpi_comm_world, status,ier)
        call MPI_RECV(tend_loc,nsta_irank*NRCOMP*leq_par%FLEX_NWIN,MPI_REAL,irank,tag,my_local_mpi_comm_world,status,ier)
        call MPI_RECV(window_chi_loc,nsta_irank*NRCOMP*NCHI*leq_par%FLEX_NWIN,&
                      MPI_REAL,irank,tag,my_local_mpi_comm_world,status,ier)
        call MPI_RECV(tr_chi_loc,nsta_irank*NRCOMP*leq_par%FLEX_NWIN,&
                      MPI_REAL,irank,tag,my_local_mpi_comm_world,status,ier)
        call MPI_RECV(am_chi_loc,nsta_irank*NRCOMP*leq_par%FLEX_NWIN,&
                      MPI_REAL,irank,tag,my_local_mpi_comm_world,status,ier)
        call MPI_RECV(T_pmax_dat_loc, nsta_irank*NRCOMP*leq_par%FLEX_NWIN,&
                      MPI_REAL, irank, tag, my_local_mpi_comm_world,status,ier)
        call MPI_RECV(T_pmax_syn_loc, nsta_irank*NRCOMP*leq_par%FLEX_NWIN,&
                      MPI_REAL, irank, tag, my_local_mpi_comm_world,status,ier)
        call MPI_RECV(imeas_loc, nsta_irank*NRCOMP*leq_par%FLEX_NWIN,&
                      MPI_INTEGER,irank,tag,my_local_mpi_comm_world,status,ier)
        do icomp=1,NRCOMP
          do irec_local = 1, nsta_irank
            do iwin=1,leq_par%FLEX_NWIN
              irec=irank*nrec_local_max+irec_local
              glob_file_prefix0(irec) = trim(file_prefix0_loc(irec_local))
              glob_net(irec) = trim(net_loc(irec_local))
              glob_sta(irec) = trim(sta_loc(irec_local))
              glob_chan_dat(irec,icomp) = trim(chan_dat_loc(irec_local,icomp))
              glob_num_win(irec,icomp) = num_win_loc(irec_local,icomp)
              glob_tstart(irec, icomp,iwin) = dble(tstart_loc(irec_local, icomp,iwin))
              glob_tend(irec, icomp,iwin) = dble(tend_loc(irec_local, icomp,iwin))
              glob_window_chi(irec, icomp,:,iwin) = dble(window_chi_loc(irec_local, icomp,:,iwin))
              glob_tr_chi(irec, icomp,iwin) = dble(tr_chi_loc(irec_local, icomp,iwin))
              glob_am_chi(irec, icomp,iwin) = dble(am_chi_loc(irec_local, icomp,iwin))
              glob_T_pmax_dat(irec, icomp,iwin) = dble(T_pmax_dat_loc(irec_local, icomp,iwin))
              glob_T_pmax_syn(irec, icomp,iwin) = dble(T_pmax_syn_loc(irec_local, icomp,iwin))
              glob_imeas(irec, icomp,iwin) = imeas_loc(irec_local, icomp,iwin)
            enddo
          enddo
        enddo
        deallocate(file_prefix0_loc)
        deallocate(net_loc)
        deallocate(sta_loc)
        deallocate(chan_dat_loc)
        deallocate(num_win_loc)
        deallocate(tstart_loc)
        deallocate(tend_loc)
        deallocate(window_chi_loc)
        deallocate(tr_chi_loc)
        deallocate(am_chi_loc)
        deallocate(T_pmax_dat_loc)
        deallocate(T_pmax_syn_loc)
        deallocate(imeas_loc)
      endif
    else
      if (myrank == irank .and. my_nrec_local > 0) then
        !write(*,*)'myrank,irank,nrec_local=',myrank,irank,nrec_local,number_receiver_global
        allocate(file_prefix0_loc(my_nrec_local))
        allocate(net_loc(my_nrec_local))
        allocate(sta_loc(my_nrec_local))
        allocate(chan_dat_loc(my_nrec_local,NRCOMP))
        allocate(num_win_loc(my_nrec_local,NRCOMP))
        allocate(tstart_loc(my_nrec_local,NRCOMP,leq_par%FLEX_NWIN))
        allocate(tend_loc(my_nrec_local,NRCOMP,leq_par%FLEX_NWIN))
        allocate(window_chi_loc(my_nrec_local,NRCOMP,NCHI,leq_par%FLEX_NWIN))
        allocate(tr_chi_loc(my_nrec_local,NRCOMP,leq_par%FLEX_NWIN))
        allocate(am_chi_loc(my_nrec_local,NRCOMP,leq_par%FLEX_NWIN))
        allocate(T_pmax_dat_loc(my_nrec_local,NRCOMP,leq_par%FLEX_NWIN))
        allocate(T_pmax_syn_loc(my_nrec_local,NRCOMP,leq_par%FLEX_NWIN))
        allocate(imeas_loc(my_nrec_local,NRCOMP,leq_par%FLEX_NWIN))

        do irec_local = 1, my_nrec_local
          irec=irank*nrec_local_max+irec_local
          !! choose the rigth seismograms_*
          do icomp=1,NRCOMP
            do iwin=1,leq_par%FLEX_NWIN
              file_prefix0_loc(irec_local)=trim(glob_file_prefix0(irec))
              net_loc(irec_local)=trim(glob_net(irec))
              sta_loc(irec_local)=trim(glob_sta(irec))
              chan_dat_loc(irec_local,icomp)=trim(glob_chan_dat(irec,icomp))
              num_win_loc(irec_local,icomp)=glob_num_win(irec,icomp)
              tstart_loc(irec_local,icomp,iwin)=real(glob_tstart(irec, icomp,iwin))
              tend_loc(irec_local,icomp,iwin)=real(glob_tend(irec,icomp,iwin))
              window_chi_loc(irec_local,icomp,1:NCHI,iwin)=real(glob_window_chi(irec, icomp, 1:NCHI,iwin))
              tr_chi_loc(irec_local,icomp,iwin)=real(glob_tr_chi(irec,icomp,iwin))
              am_chi_loc(irec_local,icomp,iwin)=real(glob_am_chi(irec,icomp,iwin))
              T_pmax_dat_loc(irec_local,icomp,iwin)=real(glob_T_pmax_dat(irec, icomp,iwin))
              T_pmax_syn_loc(irec_local,icomp,iwin)=real(glob_T_pmax_syn(irec, icomp,iwin))
              imeas_loc(irec_local,icomp,iwin)=glob_imeas(irec, icomp,iwin)
            enddo
          enddo
        enddo

        tag    = 2001
        call MPI_SEND(file_prefix0_loc,  my_nrec_local*256, MPI_CHARACTER, 0, tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(net_loc,my_nrec_local*256, MPI_CHARACTER, 0, tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(sta_loc,my_nrec_local*256, MPI_CHARACTER, 0, tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(chan_dat_loc,my_nrec_local*256*NRCOMP, MPI_CHARACTER, 0, tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(num_win_loc,my_nrec_local*NRCOMP, MPI_INTEGER, 0, tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(tstart_loc,my_nrec_local*NRCOMP*leq_par%FLEX_NWIN, MPI_REAL,0,tag,my_local_mpi_comm_world,ier)
        call MPI_SEND(tend_loc,my_nrec_local*NRCOMP*leq_par%FLEX_NWIN, MPI_REAL,0,tag, my_local_mpi_comm_world,ier)
        call MPI_SEND(window_chi_loc,my_nrec_local*NRCOMP*NCHI*leq_par%FLEX_NWIN,MPI_REAL,0,tag,my_local_mpi_comm_world,ier)
        call MPI_SEND(tr_chi_loc,my_nrec_local*NRCOMP*leq_par%FLEX_NWIN,MPI_REAL,0,tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(am_chi_loc,my_nrec_local*NRCOMP*leq_par%FLEX_NWIN,MPI_REAL,0,tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(T_pmax_dat_loc,my_nrec_local*NRCOMP*leq_par%FLEX_NWIN, MPI_REAL, 0, tag, my_local_mpi_comm_world, ier)
        call MPI_SEND(T_pmax_syn_loc,my_nrec_local*NRCOMP*leq_par%FLEX_NWIN, MPI_REAL,0,tag,my_local_mpi_comm_world, ier)
        call MPI_SEND(imeas_loc,my_nrec_local*NRCOMP*leq_par%FLEX_NWIN, MPI_INTEGER,0,tag,my_local_mpi_comm_world, ier)
        deallocate(file_prefix0_loc)
        deallocate(net_loc)
        deallocate(sta_loc)
        deallocate(chan_dat_loc)
        deallocate(num_win_loc)
        deallocate(tstart_loc)
        deallocate(tend_loc)
        deallocate(window_chi_loc)
        deallocate(tr_chi_loc)
        deallocate(am_chi_loc)
        deallocate(T_pmax_dat_loc)
        deallocate(T_pmax_syn_loc)
        deallocate(imeas_loc)
      endif
    endif
    ! not sure if need this sync
    call synchronize_all()
  enddo ! end loop irank

end subroutine collect_chi
end module