module tele_data
  use config
  use ma_constants
  use common_lib, only: get_band_name, rotate_R_to_NE_dp, dwascii, get_gauss_fac, mkdir
  use signal, only: bandpass_dp, interpolate_syn_dp, detrend, demean, &
                    myconvolution_dp, time_deconv, mycorrelation_dp, interpolate_func_dp
  use syn_data, only: SynData, average_amp_scale
  use obs_data, only: ObsData
  use input_params, fpar => fwat_par_global
  use fk_coupling
  use fwat_mpi
  use utils, only: zeros_dp, zeros, interp1
  use sacio
  use logger, only: log
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0,OUTPUT_FILES
  use decon_mod, only: deconit
  implicit none

  integer, private :: ier

  type, extends(SynData) :: TeleData
    real(kind=cr), dimension(:), pointer :: ttp ! travel time of P
    real(kind=dp), dimension(:), allocatable :: tstart, tend
    real(kind=dp), dimension(:,:), allocatable :: stf_array
    real(kind=dp), dimension(:,:), allocatable :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn
    real(kind=dp), dimension(:,:,:), allocatable :: seismo_dat, seismo_syn, window_chi
    character(len=MAX_STR_CHI), dimension(:), allocatable :: sta, net
    real(kind=cr) :: baz, az, avgamp
    integer :: ttp_win
    contains
    procedure :: preprocess, calc_fktimes, semd2sac, interp_data_to_syn
    procedure, private :: deconv_for_stf, seis_pca, measure_adj
    procedure :: finalize
  end type TeleData

contains

  subroutine semd2sac(this, ievt)
    class(TeleData), intent(inout) :: this
    integer, intent(in) :: ievt
    integer :: irec, irec_local, icomp
    logical :: is_conv_stf_local
    real(kind=dp), dimension(:), allocatable :: stf_array, seismo_syn
    real(kind=cr), dimension(:), allocatable :: bazi
    character(len=MAX_STRING_LEN) :: datafile
    type(sachead) :: header
    real(kind=cr) :: az, baz, dist, gcarc

    call this%init(ievt)

    call this%od%read_stations(ievt, .true.)
    call this%calc_fktimes()

    bazi = zeros(this%nrec)+this%baz
    call this%read(bazi)

    call mkdir(fpar%acqui%in_dat_path(this%ievt))
    call synchronize_all()

    if (.not. supp_stf) then
      if(worldrank == 0) call read_stf(stf_array)
      if (worldrank > 0) allocate(stf_array(fpar%sim%nstep))
      call bcast_all_dp(stf_array, fpar%sim%nstep)
    endif

    if (this%nrec_loc > 0) then
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
        ! read stf
        do icomp = 1, fpar%sim%NRCOMP
          if (.not. supp_stf) then
            ! convolve stf with seismogram
            call myconvolution_dp(this%data(:, icomp, irec), stf_array, seismo_syn, 0)
            seismo_syn = seismo_syn * fpar%sim%dt
          else
            seismo_syn = this%data(:, icomp, irec)
          endif
          
          ! write SAC files
          datafile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%od%netwk(irec))//'.'&
                     //trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)&
                     //trim(fpar%sim%RCOMPS(icomp))//'.sac'
          call sacio_newhead(header, real(fpar%sim%dt), fpar%sim%nstep, -real(T0))
          ! header%dist = dist
          header%az = this%az
          header%baz = this%baz
          header%stla = this%od%stla(irec)
          header%stlo = this%od%stlo(irec)
          header%stel = this%od%stel(irec)
          header%knetwk = this%od%netwk(irec)
          header%kstnm = this%od%stnm(irec)
          header%kcmpnm = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))
          header%t0 = this%ttp(irec)
          call sacio_writesac(datafile, header, seismo_syn, ier)
          if (ier /= 0) call exit_MPI(0, 'Error writing SAC file '//trim(datafile))
        enddo
      enddo
    endif
    call synchronize_all()

    contains
      subroutine read_stf(datarray)
        character(len=MAX_STRING_LEN) :: datafile
        type(sachead) :: header
        real(kind=dp), dimension(:), allocatable, intent(out) :: datarray
        real(kind=dp) :: thalf

        datafile = trim(SRC_REC_DIR)//'/STF_'//trim(fpar%acqui%evtid_names(this%ievt))//'.sac'
        call sacio_readsac(datafile, header, datarray, ier)
        if (ier /= 0) call exit_MPI(0, 'Error reading STF file '//trim(datafile))
        ! interpolate to the same time step
        if (abs(header%delta - fpar%sim%dt) > 1.0e-4 .or. header%npts /= fpar%sim%nstep) then
          thalf = (header%npts - 1) * header%delta / 2.0_dp
          datarray = interpolate_func_dp(datarray, dble(header%b), dble(header%delta), header%npts, &
                                          -thalf, dble(fpar%sim%dt), fpar%sim%nstep)
        endif
      end subroutine read_stf
  end subroutine semd2sac

  subroutine preprocess(this, ievt)
    class(TeleData), intent(inout) :: this
    integer, intent(in) :: ievt
    integer :: irec, irec_local, icomp
    real(kind=dp), dimension(:,:,:), allocatable :: seismo_dat,&
         seismo_syn, seismo_dat_glob, seismo_syn_glob
    real(kind=dp), dimension(:,:), allocatable :: seismo_stf, seismo_stf_glob
    real(kind=dp), dimension(:), allocatable :: seismo_inp, stf_local
    character(len=MAX_STRING_LEN) :: msg

    call this%init(ievt)

    ! read receiver information
    call this%od%read_stations(this%ievt)

    ! read observed data
    call this%od%read_obs_data()

    ! read syn data
    call this%read(this%od%baz)

    call get_band_name(fpar%sim%SHORT_P(1), fpar%sim%LONG_P(1), this%band_name)

    ! initialize misfits
    call this%wchi(1)%init(ievt, this%band_name)

    call this%get_comp_name_adj()

    call this%calc_fktimes()

    call synchronize_all()
    
    ! call log%write('Preprocessing', .true.)
    if (this%nrec_loc > 0) then
      this%seismo_dat = zeros_dp(fpar%sim%nstep, fpar%sim%NRCOMP, this%nrec_loc)
      this%seismo_syn = zeros_dp(fpar%sim%nstep, fpar%sim%NRCOMP, this%nrec_loc)
      seismo_stf = zeros_dp(fpar%sim%nstep, this%nrec_loc)
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
        do icomp = 1, fpar%sim%NRCOMP
          call this%interp_data_to_syn(this%od%data(:, icomp, irec), dble(this%od%tbeg(irec)),&
                                  dble(this%od%tarr(irec)), dble(this%ttp(irec)), dble(this%od%dt), seismo_inp)
          call detrend(seismo_inp)
          call demean(seismo_inp)
          call bandpass_dp(seismo_inp, fpar%sim%nstep, dble(fpar%sim%dt),&
                           1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          this%seismo_dat(:, icomp, irec_local) = seismo_inp(1:fpar%sim%nstep)

          this%seismo_syn(:, icomp, irec_local) = this%data(:, icomp, irec)
          ! call detrend(this%seismo_syn(:, icomp, irec_local))
          ! call demean(this%seismo_syn(:, icomp, irec_local))
          call bandpass_dp(this%seismo_syn(:, icomp, irec_local), fpar%sim%nstep, dble(fpar%sim%dt),&
                           1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          if (icomp == 1) then
            call this%deconv_for_stf(this%seismo_dat(:, 1, irec_local), this%seismo_syn(:, 1, irec_local),&
                                     this%ttp(irec), stf_local)
            seismo_stf(:, irec_local) = stf_local
          endif
        enddo
      enddo
    end if
    call synchronize_all()

    call this%assemble_3d(this%seismo_dat, seismo_dat_glob, fpar%sim%NRCOMP)
    call this%assemble_3d(this%seismo_syn, seismo_syn_glob, fpar%sim%NRCOMP)
    call this%assemble_2d(seismo_stf, seismo_stf_glob)

    if (worldrank == 0) then
      ! call log%write('Calculating STF', .true.)
      call this%seis_pca(seismo_dat_glob, seismo_syn_glob, seismo_stf_glob)
      this%avgamp = average_amp_scale(seismo_dat_glob, 1)
      if (IS_OUTPUT_PREPROC) then
        block
          type(sachead) :: header
          integer :: i
          call sacio_newhead(header, real(DT), NSTEP, -real(T0))
          do i = 1, fpar%sim%NRCOMP
            call sacio_writesac(trim(OUTPUT_FILES)//'/STF_'//trim(fpar%acqui%evtid_names(this%ievt))//&
                                '_'//trim(fpar%sim%RCOMPS(i))//'.sac', header, this%stf_array(:, i), ier)
            if (ier /= 0) call exit_MPI(0, 'Error writing STF file')
          enddo
        end block
      endif
    else
      allocate(this%stf_array(fpar%sim%nstep, fpar%sim%NRCOMP))
    endif
    call synchronize_all()
    call bcast_all_singlecr(this%avgamp)
    call bcast_all_dp(this%stf_array, fpar%sim%nstep*fpar%sim%NRCOMP)

    if (worldrank == 0) deallocate(seismo_dat_glob, seismo_syn_glob, seismo_stf_glob)
    
    call synchronize_all()

    call this%measure_adj()

    ! print *, shape(this%window_chi), shape(this%tr_chi), shape(this%am_chi)
    call this%wchi(1)%assemble_window_chi(this%window_chi, this%tr_chi, this%am_chi,&
                                          this%T_pmax_dat, this%T_pmax_syn, this%sta, this%net,&
                                          this%tstart, this%tend)

    if (worldrank == 0) then
      call this%wchi(1)%write()
      this%total_misfit(1) = this%wchi(1)%sum_chi(29)
      write(msg, '(a,f12.6)') 'Total misfit: of '//trim(this%band_name)//': ', this%total_misfit(1)
      call log%write(msg, .true.)
    endif
    call synchronize_all()
  end subroutine preprocess

  subroutine seis_pca(this, seismo_dat, seismo_syn, seismo_stf_glob)
    use spanlib, only: sl_pca
    class(TeleData), intent(inout) :: this
    real(kind=dp), dimension(:,:,:), intent(in) :: seismo_dat, seismo_syn
    real(kind=dp), dimension(:,:), intent(in) :: seismo_stf_glob
    real(kind=dp), dimension(:,:,:), allocatable :: recp_syn
    real(kind=dp), dimension(:), allocatable :: tmpl
    real(kind=cr), dimension(:,:), allocatable :: xeof, pc, ff
    real(kind=cr), dimension(:), allocatable :: eig, avgarr
    real(kind=cr) :: valid_avg, valid_sum, avgamp
    integer :: irec, icomp, valid_count, i

    this%stf_array = zeros_dp(fpar%sim%nstep, fpar%sim%NRCOMP)
    allocate(xeof(this%nrec, this%nrec))
    allocate(pc(fpar%sim%nstep, this%nrec))
    allocate(eig(this%nrec))
    ff = transpose(real(seismo_stf_glob(:, :)))
    call sl_pca(ff, this%nrec, xeof, pc, eig)

    ! build reconstructed seismograms
    recp_syn = zeros_dp(fpar%sim%nstep, fpar%sim%NRCOMP, this%nrec)
    do irec = 1, this%nrec
      do icomp = 1, fpar%sim%NRCOMP
        call myconvolution_dp(seismo_syn(:, icomp, irec), dble(pc(:, 1)), tmpl, 0)
        recp_syn(:, icomp, irec) = tmpl * fpar%sim%dt
      enddo
    enddo

    ! calculate amplitude correction
    do icomp = 1, fpar%sim%NRCOMP
      avgarr = zeros(this%nrec)
      avgarr = -1000.0
      valid_sum = 0.0
      valid_count = 0

      where (maxval(abs(recp_syn(:, icomp, :)), dim=1) /= 0.0)
        avgarr = sum(seismo_dat(:, icomp, :) * recp_syn(:, icomp, :), dim=1) / &
                 sum(recp_syn(:, icomp, :) * recp_syn(:, icomp, :), dim=1)
      end where
      do i = 1, size(avgarr)
        if (avgarr(i) /= -1000.0) then
          valid_sum = valid_sum + avgarr(i)
          valid_count = valid_count + 1
        endif
      enddo
      if (valid_count == 0) cycle

      valid_avg = valid_sum / valid_count
      valid_sum = 0.0
      valid_count = 0
      do i = 1, this%nrec
        if (abs(avgarr(i) - valid_avg) <= 0.2 .and. avgarr(i) /= -1000.0) then
          valid_sum = valid_sum + avgarr(i)
          valid_count = valid_count + 1
        endif
      enddo

      if (valid_count == 0 .or. valid_count <= 0.1 * real(this%nrec)) then
        call exit_MPI(0, 'Error: Too few valid amplitude values')
      endif
      avgamp = valid_sum / valid_count
      this%stf_array(:, icomp) = dble(pc(:, 1) * avgamp)
    enddo
    deallocate(xeof, pc, eig, recp_syn, ff, tmpl, avgarr)

  end subroutine seis_pca

  subroutine deconv_for_stf(this, data_num, data_den, tref, stf)
    class(TeleData), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: data_num, data_den
    real(kind=dp), dimension(:), allocatable, intent(out) :: stf
    real(kind=cr), dimension(:), allocatable :: stf_dp
    real(kind=cr), intent(in) :: tref
    real(kind=cr) :: tb, te, f0, time_shift
    integer :: nstep_cut
    real(kind=dp), dimension(:), allocatable :: data_num_win, data_den_win

    ! cut data to from fpar%sim%time_win(1) to fpar%sim%time_win(2)
    data_num_win = data_num(1:fpar%sim%nstep)
    data_den_win = data_den(1:fpar%sim%nstep)
    tb = tref + fpar%sim%time_win(1)
    te = tref + fpar%sim%time_win(2)
    time_shift = (NSTEP / 2) * DT
    nstep_cut = int((te - tb) / fpar%sim%dt) + 1
    call interpolate_syn_dp(data_num_win, -dble(T0), dble(fpar%sim%dt), fpar%sim%nstep, &
                            dble(tb), dble(fpar%sim%dt), nstep_cut)
    call interpolate_syn_dp(data_den_win, -dble(T0), dble(fpar%sim%dt), fpar%sim%nstep, &
                            dble(tb), dble(fpar%sim%dt), nstep_cut)
    call interpolate_syn_dp(data_num_win, dble(tb), dble(fpar%sim%dt), nstep_cut, &
                            -dble(T0), dble(fpar%sim%dt), fpar%sim%nstep)
    call interpolate_syn_dp(data_den_win, dble(tb), dble(fpar%sim%dt), nstep_cut, &
                            -dble(T0), dble(fpar%sim%dt), fpar%sim%nstep)
    if (maxval(abs(data_den_win)) < 1.0e-10) &
      call exit_MPI(0, 'Error: data_den_win is zero')
    call time_deconv(real(data_num_win),real(data_den_win),fpar%sim%dt,&
                     fpar%sim%nstep,NITER,stf_dp)
    stf = real(stf_dp)
    ! f0 = dble(get_gauss_fac(1/fpar%sim%SHORT_P(1))) * 4
    ! call deconit(data_num_win, data_den_win, fpar%sim%dt, time_shift, f0, NITER, 0.001, 1, stf)
    call bandpass_dp(stf, fpar%sim%nstep, dble(fpar%sim%dt),&
                     1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)

  end subroutine deconv_for_stf

  subroutine calc_fktimes(this)
    class(TeleData), intent(inout) :: this

    ! calculate fault rupture times
    real(kind=dp) :: stla, stlo
    real(kind=dp) :: xx, yy, zz
    integer :: irec

    call prepare_shm_array_cr_1d(this%ttp, this%nrec, this%ttp_win)

    call read_fk_model(fpar%acqui%evtid_names(this%ievt))

    this%baz = -phi_FK - 90.d0
    this%az = 90.d0 - phi_FK
    
    if (worldrank == 0) then
      do irec = 1, this%nrec
        stlo = dble(this%od%stlo(irec))
        stla = dble(this%od%stla(irec))
        if (.not. SUPPRESS_UTM_PROJECTION) then
          call utm_geo(stlo, stla, xx,yy, ILONGLAT2UTM)
        else
          xx = this%od%stlo(irec)
          yy = this%od%stla(irec)
        endif
        zz = this%od%stel(irec)
        call fktime(real(xx), real(yy), real(zz), this%ttp(irec))
      end do
      this%ttp = this%ttp - T0
    endif
    call sync_from_main_rank_cr_1d(this%ttp, this%nrec)
    call free_fk_arrays()
  end subroutine calc_fktimes

  subroutine finalize(this)
    class(TeleData), intent(inout) :: this

    call this%od%finalize()
    call this%wchi(1)%finalize()
    call free_shm_array(this%ttp_win)
    call free_shm_array(this%dat_win)
    deallocate(this%stf_array, stat=ier)
    deallocate(this%seismo_dat, stat=ier)
    deallocate(this%seismo_syn, stat=ier)
    deallocate(this%window_chi, stat=ier)
    deallocate(this%tr_chi, stat=ier)
    deallocate(this%am_chi, stat=ier)
    deallocate(this%T_pmax_dat, stat=ier)
    deallocate(this%T_pmax_syn, stat=ier)
    deallocate(this%sta, stat=ier)
    deallocate(this%net, stat=ier)
    deallocate(this%tstart, stat=ier)
    deallocate(this%tend, stat=ier)

  end subroutine finalize

  subroutine interp_data_to_syn(this, dat_in, tb, tarr, ttp, deltat, dat_out)
    class(TeleData), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: dat_in
    real(kind=dp), intent(in) :: tarr, deltat, tb, ttp
    real(kind=dp), dimension(:), allocatable, intent(out) :: dat_out
    real(kind=dp), dimension(:), allocatable :: time_in, time_out
    integer :: i
    
    time_in = [ (dble(i-1) * deltat, i=1, size(dat_in)) ] + tb - tarr + ttp
    time_out = [ (dble(i-1) * fpar%sim%dt, i=1, fpar%sim%nstep) ] - T0
    dat_out = interp1(time_in, dat_in, time_out, 0.0_dp)

  end subroutine interp_data_to_syn

  subroutine measure_adj(this)
    use measure_adj_mod, only: measure_adj_fwat
    use ma_variables, only: OUT_DIR
    class(TeleData), intent(inout) :: this
    type(sachead) :: header
    real(kind=dp) :: tstart, tend
    character(len=MAX_STRING_LEN) :: sacfile, file_prefix
    real(kind=dp), dimension(:), allocatable :: seismo_syn_local, tmpl
    real(kind=dp), dimension(:, :), allocatable :: adj_src
    real(kind=dp), dimension(NCHI) :: window_chi_local
    real(kind=dp), dimension(NDIM_MA) :: adj_syn_local, dat_meas, syn_meas
    ! real(kind=dp), dimension(:), allocatable :: adj_syn_local
    real(kind=dp) :: tr_chi_local, am_chi_local, T_pmax_dat_local, T_pmax_syn_local
    integer :: out_imeas, irec, icomp, irec_local
    
    ! allocate chi arrays
    OUT_DIR = trim(OUTPUT_FILES)
    allocate(this%window_chi(this%nrec_loc, NCHI, fpar%sim%NRCOMP))
    allocate(this%tr_chi(this%nrec_loc, fpar%sim%NRCOMP))
    allocate(this%am_chi(this%nrec_loc, fpar%sim%NRCOMP))
    allocate(this%T_pmax_dat(this%nrec_loc, fpar%sim%NRCOMP))
    allocate(this%T_pmax_syn(this%nrec_loc, fpar%sim%NRCOMP))
    allocate(this%sta(this%nrec_loc))
    allocate(this%net(this%nrec_loc))
    allocate(this%tstart(this%nrec_loc))
    allocate(this%tend(this%nrec_loc))
    if (this%nrec_loc > 0) then
      do irec_local = 1, this%nrec_loc
        irec = select_global_id_for_rec(irec_local)
        do icomp = 1, fpar%sim%NRCOMP
          ! convolution with STF
          call myconvolution_dp(this%seismo_syn(:, icomp, irec_local), this%stf_array(:, icomp), seismo_syn_local, 0)
          this%seismo_syn(:, icomp, irec_local) = seismo_syn_local * fpar%sim%dt
          if (IS_OUTPUT_PREPROC) then
            call sacio_newhead(header, real(DT), NSTEP, -real(T0))
            header%knetwk = trim(this%od%netwk(irec))
            header%kstnm = trim(this%od%stnm(irec))
            header%t0 = this%ttp(irec)
            ! header%kevnm = trim(fpar%acqui%evtid_names(this%ievt))
            sacfile = trim(OUTPUT_FILES)//'/wsyn.'//trim(this%od%netwk(irec))//&
                      '.'//trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
                      trim(fpar%sim%RCOMPS(icomp))//'.sac.'//trim(this%band_name)
            call sacio_writesac(sacfile, header, this%seismo_syn(:, icomp, irec_local), ier)
            sacfile = trim(OUTPUT_FILES)//'/wdat.'//trim(this%od%netwk(irec))//&
                      '.'//trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
                      trim(fpar%sim%RCOMPS(icomp))//'.sac.'//trim(this%band_name)
            call sacio_writesac(sacfile, header, this%seismo_dat(:, icomp, irec_local), ier)
            ! sacfile = trim(OUTPUT_FILES)//'/diff.'//trim(this%od%netwk(irec))//&
            !           '.'//trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
            !           trim(fpar%sim%RCOMPS(icomp))//'.sac.'//trim(this%band_name)
            ! call sacio_writesac(sacfile, header, this%seismo_syn(:, icomp, irec_local)- this%seismo_dat(:, icomp, irec_local), ier)
          endif

          ! measure adjoint
          tstart = this%ttp(irec) + fpar%sim%time_win(1)
          tend = this%ttp(irec) + fpar%sim%time_win(2)
          dat_meas(1:NSTEP) = this%seismo_dat(:, icomp, irec_local)
          syn_meas(1:NSTEP) = this%seismo_syn(:, icomp, irec_local)
          call measure_adj_fwat(dat_meas, syn_meas,&
                                tstart, tend, dble(fpar%sim%SHORT_P(1)), dble(fpar%sim%LONG_P(1)),&
                                this%od%netwk(irec), this%od%stnm(irec),&
                                trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp)), &
                                window_chi_local, tr_chi_local, am_chi_local,&
                                T_pmax_dat_local, T_pmax_syn_local, adj_syn_local, &
                                file_prefix, out_imeas)
          ! call measure_adj_tele(this%seismo_dat(:, icomp, irec_local), this%seismo_syn(:, icomp, irec_local),&
          !                        tstart, tend, window_chi_local, adj_syn_local)

          if (IS_OUTPUT_ADJ_SRC) then
            call sacio_newhead(header, real(DT), NSTEP, -real(T0))
            header%knetwk = trim(this%od%netwk(irec))
            header%kstnm = trim(this%od%stnm(irec))
            header%kcmpnm = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))
            sacfile = trim(fpar%acqui%out_fwd_path(this%ievt))//'/'//trim(ADJOINT_PATH)//&
                      '/'//trim(this%od%netwk(irec))//'.'//trim(this%od%stnm(irec))//&
                      '.'//trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))//&
                      '.adj.sac.'//trim(this%band_name)
            call sacio_writesac(sacfile, header, adj_syn_local(1:NSTEP), ier)
          endif
          adj_syn_local = adj_syn_local / this%avgamp
          call mycorrelation_dp(adj_syn_local, this%stf_array(:, icomp), tmpl, 0)
          adj_syn_local = tmpl * fpar%sim%dt
          adj_syn_local = adj_syn_local / maxval(abs(this%stf_array(:, icomp)))

          ! write adjoint source
          select case (icomp)
          case (1)
            call this%write_adj(adj_syn_local(1:NSTEP), trim(this%comp_name(1)), irec)
          case (2)
            adj_src = zeros_dp(fpar%sim%nstep, 2)
            call rotate_R_to_NE_dp(adj_syn_local(1:NSTEP), adj_src(:, 2), adj_src(:, 1), this%baz)
            call this%write_adj(adj_src(:, 1), trim(this%comp_name(2)), irec)
            call this%write_adj(adj_src(:, 2), trim(this%comp_name(3)), irec)
          end select
          ! save misfits
          this%window_chi(irec_local, :, icomp) = window_chi_local
          tr_chi_local = window_chi_local(15)
          am_chi_local = window_chi_local(15)
          this%tr_chi(irec_local, icomp) = fpar%acqui%src_weight(this%ievt)*tr_chi_local/this%avgamp/this%avgamp*dt
          this%am_chi(irec_local, icomp) = fpar%acqui%src_weight(this%ievt)*am_chi_local/this%avgamp/this%avgamp*dt
          ! this%T_pmax_dat(irec_local, icomp) = T_pmax_dat_local
          this%T_pmax_dat(irec_local, icomp) = (dt*NSTEP)/maxloc(abs(this%seismo_dat(:, icomp, irec_local)), dim=1)
          ! this%T_pmax_syn(irec_local, icomp) = T_pmax_syn_local
          this%T_pmax_syn(irec_local, icomp) = (dt*NSTEP)/maxloc(abs(this%seismo_syn(:, icomp, irec_local)), dim=1)
          if (icomp == 1) then
            this%sta(irec_local) = this%od%stnm(irec)
            this%net(irec_local) = this%od%netwk(irec)
            this%tstart(irec_local) = tstart
            this%tend(irec_local) = tend
          endif
        enddo
      enddo
    end if
    call synchronize_all()
  end subroutine measure_adj
end module tele_data