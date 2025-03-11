module tele_data
  use config
  use ma_constants
  use common_lib, only: get_band_name
  use signal, only: bandpass_dp, interpolate_syn_dp, detrend, demean, &
                    myconvolution_dp, time_deconv
  use syn_data, only: SynData, average_amp_scale
  use obs_data, only: ObsData
  use input_params, fpar => fwat_par_global
  use fk_coupling
  use fwat_mpi
  use utils, only: zeros_dp, zeros, interp1
  use sacio
  use logger, only: log
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, nrec, nrec_local,OUTPUT_FILES, &
                        number_receiver_global, ispec_selected_rec, islice_selected_rec

  implicit none

  integer, private :: ier

  type, extends(SynData) :: TeleData
    real(kind=cr), dimension(:), pointer :: ttp ! travel time of P
    real(kind=cr) :: baz, az
    integer :: ttp_win
    contains
    procedure :: preprocess, calc_fktimes, semd2sac!, calc_stf
    procedure, private :: deconv_for_stf, seis_pca, interp_data_to_syn
    procedure :: finalize
  end type TeleData

contains

  subroutine semd2sac(this, ievt, is_conv_stf)
    class(TeleData), intent(inout) :: this
    integer, intent(in) :: ievt
    logical, intent(in), optional :: is_conv_stf
    integer :: irec, irec_local, icomp
    logical :: is_conv_stf_local
    real(kind=dp), dimension(:), allocatable :: stf_array, seismo_syn
    character(len=MAX_STRING_LEN) :: datafile
    type(sachead) :: header
    real(kind=cr) :: az, baz, dist, gcarc

    if (present(is_conv_stf)) then
      is_conv_stf_local = is_conv_stf
    else
      is_conv_stf_local = .false.
    endif

    call this%read(ievt)

    call this%od%read_stations(this%ievt)
    
    call this%calc_fktimes()

    if (worldrank == 0) call system('mkdir -p '//trim(fpar%acqui%in_dat_path(this%ievt)))
    call synchronize_all()

    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        if (is_conv_stf_local) then
          ! read stf
          call read_stf(stf_array)
        endif
        do icomp = 1, fpar%sim%NRCOMP
          if (is_conv_stf_local) then
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
          ! header%gcarc = gcarc
          ! header%evla = fpar%acqui%evla(this%ievt)
          ! header%evlo = fpar%acqui%evlo(this%ievt)
          ! header%evdp = fpar%acqui%evdp(this%ievt)
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
        if (abs(header%delta - fpar%sim%dt) > 1.0e-3 .or. header%npts /= fpar%sim%nstep) then
          thalf = (header%npts - 1) * header%delta / 2.0_dp
          call interpolate_syn_dp(datarray, -thalf, dble(header%delta), header%npts, &
                                  -dble(T0), dble(fpar%sim%dt), fpar%sim%nstep)
        endif
      end subroutine read_stf
  end subroutine semd2sac

  subroutine preprocess(this, ievt)
    use measure_adj_mod, only: measure_adj_fwat
    class(TeleData), intent(inout) :: this
    integer, intent(in) :: ievt
    integer :: irec, irec_local, icomp
    real(kind=dp), dimension(:,:,:), allocatable :: seismo_dat,&
         seismo_syn, seismo_dat_glob, seismo_syn_glob, adj_src
    real(kind=dp), dimension(:,:), allocatable :: seismo_stf, seismo_stf_glob, stf_array
    real(kind=dp), dimension(:), allocatable :: seismo_inp, stf_local
    real(kind=cr) :: avgamp

    ! read syn data
    call this%read(ievt)

    ! read receiver information
    call this%od%read_stations(this%ievt)

    ! read observed data
    call this%od%read_obs_data()

    ! initialize misfits
    call this%wchi%init()

    call this%calc_fktimes()

    call get_band_name(fpar%sim%SHORT_P(1), fpar%sim%LONG_P(1), this%band_name)

    call synchronize_all()
    
    call log%write('Preprocessing', .true.)
    if (nrec_local > 0) then
      seismo_dat = zeros_dp(fpar%sim%nstep, fpar%sim%NRCOMP, nrec_local)
      seismo_syn = zeros_dp(fpar%sim%nstep, fpar%sim%NRCOMP, nrec_local)
      seismo_stf = zeros_dp(fpar%sim%nstep, nrec_local)
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        do icomp = 1, fpar%sim%NRCOMP
          ! seismo_inp = this%od%data(:, icomp, irec)
          ! call interpolate_syn_dp(seismo_inp, dble(this%od%tarr(irec)), dble(this%od%dt),&
                                  ! this%od%npts, -dble(T0), dble(fpar%sim%dt), fpar%sim%nstep)
          call this%interp_data_to_syn(this%od%data(:, icomp, irec), dble(this%od%tbeg(irec)),&
                                  dble(this%od%tarr(irec)), dble(this%ttp(irec)), dble(this%od%dt), seismo_inp)
          call detrend(seismo_inp)
          call demean(seismo_inp)
          call bandpass_dp(seismo_inp, fpar%sim%nstep, dble(fpar%sim%dt),&
                           1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          seismo_dat(:, icomp, irec_local) = seismo_inp(1:fpar%sim%nstep)

          seismo_syn(:, icomp, irec_local) = this%data(:, icomp, irec)
          call detrend(seismo_syn(:, icomp, irec_local))
          call demean(seismo_syn(:, icomp, irec_local))
          call bandpass_dp(seismo_syn(:, icomp, irec_local), fpar%sim%nstep, dble(fpar%sim%dt),&
                           1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          if (icomp == 1) then
            call this%deconv_for_stf(seismo_dat(:, 1, irec_local), seismo_syn(:, 1, irec_local),&
                                     this%ttp(irec), stf_local)
            seismo_stf(:, irec_local) = stf_local
          endif
        enddo
      enddo
    end if
    call synchronize_all()

    call this%assemble_3d(seismo_dat, seismo_dat_glob, fpar%sim%NRCOMP)
    call this%assemble_3d(seismo_syn, seismo_syn_glob, fpar%sim%NRCOMP)
    call this%assemble_2d(seismo_stf, seismo_stf_glob)
    
    if (worldrank == 0) then
      call log%write('Calculating STF', .true.)
      call this%seis_pca(seismo_dat_glob, seismo_syn_glob, seismo_stf_glob, stf_array)
      avgamp = average_amp_scale(seismo_dat_glob, 1)
      if (IS_OUTPUT_PREPROC) then
        block
          type(sachead) :: header
          integer :: i
          call sacio_newhead(header, real(DT), NSTEP, -real(T0))
          do i = 1, fpar%sim%NRCOMP
            call sacio_writesac(trim(OUTPUT_FILES)//'/STF_'//trim(fpar%acqui%evtid_names(this%ievt))//&
                                '_'//trim(fpar%sim%RCOMPS(i))//'.sac', header, stf_array(:, i), ier)
            if (ier /= 0) call exit_MPI(0, 'Error writing STF file')
          enddo
        end block
      endif
    else
      allocate(stf_array(fpar%sim%nstep, fpar%sim%NRCOMP))
    endif
    call synchronize_all()
    call bcast_all_singlecr(avgamp)
    call bcast_all_dp(stf_array, fpar%sim%nstep*fpar%sim%NRCOMP)

    if (worldrank == 0) deallocate(seismo_dat_glob, seismo_syn_glob, seismo_stf_glob)

    ! convolution with STF
    if (nrec_local > 0) then
      adj_src = zeros_dp(fpar%sim%nstep, fpar%sim%NRCOMP, nrec_local)
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        do icomp = 1, fpar%sim%NRCOMP
          block 
            type(sachead) :: header
            real(kind=dp) :: tstart, tend
            character(len=MAX_STRING_LEN) :: sacfile, file_prefix
            real(kind=dp), dimension(:), allocatable :: seismo_syn_local
            real(kind=dp), dimension(NCHI) :: window_chi
            real(kind=dp), dimension(NDIM_MA) :: adj_src_local
            real(kind=dp) :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn
            integer :: out_imeas

            call myconvolution_dp(seismo_syn(:, icomp, irec_local), stf_array(:, icomp), seismo_syn_local, 0)
            seismo_syn(:, icomp, irec_local) = seismo_syn_local * fpar%sim%dt
            if (IS_OUTPUT_PREPROC) then
              call sacio_newhead(header, real(DT), NSTEP, -real(T0))
              sacfile = trim(OUTPUT_FILES)//'/wsyn'//trim(this%od%netwk(irec))//&
                        '.'//trim(this%od%stnm(irec))//&
                        '.'//trim(fpar%sim%RCOMPS(icomp))//&
                        '.sac.'//trim(this%band_name)
              call sacio_writesac(sacfile, header, seismo_syn(:, icomp, irec_local), ier)
              sacfile = trim(OUTPUT_FILES)//'/wdat'//trim(this%od%netwk(irec))//&
                        '.'//trim(this%od%stnm(irec))//&
                        '.'//trim(fpar%sim%RCOMPS(icomp))//&
                        '.sac.'//trim(this%band_name)
              call sacio_writesac(sacfile, header, seismo_dat(:, icomp, irec_local), ier)
            endif

            tstart = this%ttp(irec) - fpar%sim%time_win(1)
            tend = this%ttp(irec) + fpar%sim%time_win(2)
            call measure_adj_fwat(seismo_dat(:, icomp, irec_local), seismo_syn(:, icomp, irec_local),&
                                  tstart, tend, dble(-t0),dble(DT), NSTEP, this%od%netwk(irec), this%od%stnm(irec),&
                                  fpar%sim%CH_CODE, window_chi, tr_chi,am_chi, T_pmax_dat,T_pmax_syn,&
                                  adj_src_local, file_prefix, out_imeas, this%band_name)
            ! adj_src(:, icomp, irec_local) = adj_src_local(1:NSTEP)
          end block
        enddo
      enddo
    end if
    call synchronize_all()
  end subroutine preprocess

  subroutine seis_pca(this, seismo_dat, seismo_syn, seismo_stf_glob, stf_array)
    use spanlib, only: sl_pca
    class(TeleData), intent(inout) :: this
    real(kind=dp), dimension(:,:,:), intent(in) :: seismo_dat, seismo_syn
    real(kind=dp), dimension(:,:), intent(in) :: seismo_stf_glob
    real(kind=dp), dimension(:,:,:), allocatable :: recp_syn
    real(kind=dp), dimension(:), allocatable :: tmpl
    real(kind=cr), dimension(:,:), allocatable :: xeof, pc, ff
    real(kind=cr), dimension(:), allocatable :: eig, avgarr
    real(kind=dp), dimension(:,:), allocatable, intent(out) :: stf_array
    real(kind=cr) :: valid_avg, valid_sum, avgamp
    integer :: irec, icomp, valid_count, i

    stf_array = zeros_dp(fpar%sim%nstep, fpar%sim%NRCOMP)
    allocate(xeof(nrec, nrec))
    allocate(pc(fpar%sim%nstep, nrec))
    allocate(eig(nrec))
    ff = transpose(real(seismo_stf_glob(:, :)))
    call sl_pca(ff, nrec, xeof, pc, eig)

    ! build reconstructed seismograms
    recp_syn = zeros_dp(fpar%sim%nstep, fpar%sim%NRCOMP, nrec)
    do irec = 1, nrec
      do icomp = 1, fpar%sim%NRCOMP
        call myconvolution_dp(seismo_syn(:, icomp, irec), dble(pc(:, 1)), tmpl, 0)
        recp_syn(:, icomp, irec) = tmpl * fpar%sim%dt
      enddo
    enddo

    ! calculate amplitude correction
    do icomp = 1, fpar%sim%NRCOMP
      avgarr = zeros(nrec)
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
      do i = 1, nrec
        if (abs(avgarr(i) - valid_avg) <= 0.2 .and. avgarr(i) /= -1000.0) then
          valid_sum = valid_sum + avgarr(i)
          valid_count = valid_count + 1
        endif
      enddo

      if (valid_count == 0 .or. valid_count <= 0.1 * real(nrec)) then
        call exit_MPI(0, 'Error: Too few valid amplitude values')
      endif
      avgamp = valid_sum / valid_count
      stf_array(:, icomp) = dble(pc(:, 1) * avgamp)
    enddo
    deallocate(xeof, pc, eig, recp_syn, ff, tmpl, avgarr)

  end subroutine seis_pca

  subroutine deconv_for_stf(this, data_num, data_den, tref, stf)
    class(TeleData), intent(inout) :: this
    real(kind=dp), dimension(:), intent(in) :: data_num, data_den
    real(kind=dp), dimension(:), allocatable, intent(inout) :: stf
    real(kind=cr), dimension(:), allocatable :: stf_dp
    real(kind=cr), intent(in) :: tref
    real(kind=cr) :: tb, te
    integer :: nstep_cut
    real(kind=dp), dimension(:), allocatable :: data_num_win, data_den_win

    ! cut data to from fpar%sim%time_win(1) to fpar%sim%time_win(2)
    data_num_win = data_num(1:fpar%sim%nstep)
    data_den_win = data_den(1:fpar%sim%nstep)
    tb = tref - fpar%sim%time_win(1)
    te = tref + fpar%sim%time_win(2)
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
    stf = dble(stf_dp)
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
    call free_shm_array(this%ttp_win)
    call free_shm_array(this%dat_win)
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

end module tele_data