module telecc_data
  use config
  use measure_adj_mod, only: meas_adj_conv_diff
  use ma_constants
  use common_lib, only: get_band_name, rotate_R_to_NE_dp, dwascii, &
                        get_gauss_fac
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
  use specfem_par, only: T0, DT, NSTEP, nrec, nrec_local,OUTPUT_FILES, &
                        number_receiver_global, ispec_selected_rec, islice_selected_rec
  use tele_data, only: TeleData

  implicit none

  integer, private :: ier, ncomp

  type, extends(TeleData) :: TeleCCData

  contains
  procedure :: preprocess
  procedure, private :: measure_adj, pre_proc
  end type TeleCCData

contains

  subroutine preprocess(this, ievt)
    class(TeleCCData), intent(inout) :: this
    integer, intent(in) :: ievt
    integer :: irec_local, irec
    character(len=MAX_STRING_LEN) :: msg

    fpar%sim%NRCOMP = 1
    ncomp = size(fpar%sim%RCOMPS)

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

    call this%pre_proc()

    call this%measure_adj()

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

  subroutine pre_proc(this)
    class(TeleCCData), intent(inout) :: this
    integer :: irec_local, irec, nstep_cut, icomp
    real(kind=dp), dimension(:), allocatable :: seismo_inp
    real(kind=dp) :: tstart, tend

    if (nrec_local > 0) then
      this%seismo_dat = zeros_dp(NSTEP, ncomp, nrec_local)
      this%seismo_syn = zeros_dp(NSTEP, ncomp, nrec_local)
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        tstart = this%ttp(irec) - fpar%sim%time_win(1)
        tend = this%ttp(irec) + fpar%sim%time_win(2)
        nstep_cut = int((tend - tstart) / dble(DT)) + 1
        do icomp = 1, size(fpar%sim%RCOMPS)
          call this%interp_data_to_syn(this%od%data(:, icomp, irec), dble(this%od%tbeg(irec)),&
                                  dble(this%od%tarr(irec)), dble(this%ttp(irec)), dble(this%od%dt), seismo_inp)
          call detrend(seismo_inp)
          call demean(seismo_inp)
          call bandpass_dp(seismo_inp, NSTEP, dble(DT),&
                           1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          call interpolate_syn_dp(seismo_inp, -dble(T0), dble(DT), NSTEP, &
                                  tstart, dble(dt), nstep_cut)
          this%seismo_dat(:, icomp, irec_local) = seismo_inp

          seismo_inp = this%data(:, icomp, irec)
          call detrend(seismo_inp)
          call demean(seismo_inp)
          call bandpass_dp(seismo_inp, NSTEP, dble(DT),&
                           1/fpar%sim%LONG_P(1), 1/fpar%sim%SHORT_P(1), IORD)
          call interpolate_syn_dp(seismo_inp, -dble(T0), dble(DT), NSTEP, &
                                  tstart, dble(dt), nstep_cut)
          this%seismo_syn(:, icomp, irec_local) = seismo_inp
        enddo
      enddo
    endif
    call synchronize_all()
  
  end subroutine pre_proc

  subroutine measure_adj(this)
    class(TeleCCData), intent(inout) :: this
    integer :: irec_local, irec
    real(kind=dp), dimension(:), allocatable :: adj_r_tw, adj_z_tw
    real(kind=dp), dimension(:,:), allocatable :: adj_src
    real(kind=dp), dimension(NCHI) :: window_chi_local
    real(kind=dp) :: f0
    type(sachead) :: header

    if (nrec_local > 0) then
      allocate(this%window_chi(nrec_local, NCHI, fpar%sim%NRCOMP))
      allocate(this%tr_chi(nrec_local, fpar%sim%NRCOMP))
      allocate(this%am_chi(nrec_local, fpar%sim%NRCOMP))
      allocate(this%T_pmax_dat(nrec_local, fpar%sim%NRCOMP))
      allocate(this%T_pmax_syn(nrec_local, fpar%sim%NRCOMP))
      allocate(this%sta(nrec_local))
      allocate(this%net(nrec_local))
      allocate(this%tstart(nrec_local))
      allocate(this%tend(nrec_local))
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        if (IS_OUTPUT_PREPROC) then
          block
            character(len=MAX_STRING_LEN) :: sacfile
            real(kind=dp), dimension(:), allocatable :: conv
            integer :: nb
            nb = floor(fpar%sim%TIME_WIN(1)/DT)
            call myconvolution_dp(this%seismo_syn(:, 2, irec_local),&
                                  this%seismo_dat(:, 1, irec_local),&
                                  conv, 1)
            call sacio_newhead(header, real(DT), NSTEP, -real(T0))
            header%knetwk = trim(this%od%netwk(irec))
            header%kstnm = trim(this%od%stnm(irec))
            header%kcmpnm = trim(fpar%sim%CH_CODE)//'Z'
            sacfile = trim(OUTPUT_FILES)//'/wconv.'//trim(this%od%netwk(irec))//&
                      '.'//trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
                      'Z.sac.'//trim(this%band_name)
            call sacio_writesac(sacfile, header, conv(nb:NSTEP+nb-1), ier)
            header%kcmpnm = trim(fpar%sim%CH_CODE)//'R'
            call myconvolution_dp(this%seismo_syn(:, 1, irec_local),&
                                  this%seismo_dat(:, 2, irec_local),&
                                  conv, 1)
            sacfile = trim(OUTPUT_FILES)//'/wconv.'//trim(this%od%netwk(irec))//&
                      '.'//trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
                      'R.sac.'//trim(this%band_name)
            call sacio_writesac(sacfile, header, conv(nb:NSTEP+nb-1), ier)
          end block
        endif

        f0 = dble(get_gauss_fac(1/fpar%sim%SHORT_P(1)))
        call meas_adj_conv_diff(this%seismo_dat(:, 2, irec_local), this%seismo_dat(:, 1, irec_local),&
                                this%seismo_syn(:, 2, irec_local), this%seismo_dat(:, 1, irec_local),&
                                -dble(fpar%sim%TIME_WIN(1)), dble(fpar%sim%TIME_WIN(1)), dble(this%ttp(irec)),&
                                NSTEP, dble(fpar%sim%RF%tshift), f0, NITER, 0.001_dp, &
                                dble(1/fpar%sim%LONG_P(1)), dble(1/fpar%sim%SHORT_P(1)), &
                                window_chi_local, adj_r_tw, adj_z_tw)

        call this%write_adj(adj_z_tw(1:NSTEP), trim(this%comp_name(1)), irec)
        adj_src = zeros_dp(NSTEP, 2)
        call rotate_R_to_NE_dp(adj_r_tw(1:NSTEP), adj_src(:, 2), adj_src(:, 1), this%baz)
        call this%write_adj(adj_src(:, 1), trim(this%comp_name(2)), irec)
        call this%write_adj(adj_src(:, 2), trim(this%comp_name(3)), irec)

        this%window_chi(irec_local, :, 1) = window_chi_local
        this%tr_chi(irec_local, 1) = fpar%acqui%src_weight(this%ievt)*window_chi_local(15)
        this%am_chi(irec_local, 1) = fpar%acqui%src_weight(this%ievt)*window_chi_local(15)
        this%T_pmax_dat(irec_local, 1) = 0.
        this%T_pmax_syn(irec_local, 1) = 0.
        this%sta(irec_local) = this%od%stnm(irec)
        this%net(irec_local) = this%od%netwk(irec)
        this%tstart(irec_local) = this%ttp(irec) - fpar%sim%time_win(1)
        this%tend(irec_local) = this%ttp(irec) + fpar%sim%time_win(2)
      enddo
    endif
    call synchronize_all()
  end subroutine measure_adj

end module telecc_data