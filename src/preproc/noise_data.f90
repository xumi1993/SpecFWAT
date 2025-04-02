module noise_data
  use config
  use fwat_constants
  use ma_constants
  use syn_data, only: SynData
  use input_params, fpar => fwat_par_global
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, DT, NSTEP, nrec, nrec_local,OUTPUT_FILES, &
                        number_receiver_global, ispec_selected_rec, islice_selected_rec
  use fwat_mpi
  use utils, only: zeros_dp, zeros
  use distaz_lib
  use common_lib, only: get_band_name, get_icomp_syn, rotate_ZRT_to_ZNE
  use signal, only: interpolate_func_dp, dif1, detrend, demean, bandpass_dp
  use sacio
  use logger, only: log


  implicit none
  integer, private :: ier
  type, extends(SynData) :: NoiseData
  contains
    procedure :: semd2sac, preprocess, finalize
    procedure, private :: calc_distaz, measure_adj, write_in_preocess
  end type NoiseData

contains
  subroutine semd2sac(this, ievt)
    class(NoiseData), intent(inout) :: this
    integer, intent(in) :: ievt
    integer :: irec_local, irec, icomp, icomp_syn
    character(len=MAX_STRING_LEN) :: datafile
    type(sachead) :: header

    call this%init(ievt)

    if (worldrank == 0) call system('mkdir -p '//trim(fpar%acqui%in_dat_path(this%ievt)))

    call this%od%read_stations(ievt, .true.)

    call this%calc_distaz()

    call this%read(this%od%baz)

    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        do icomp = 1, fpar%sim%NRCOMP
          icomp_syn = get_icomp_syn(fpar%sim%RCOMPS(icomp))
          datafile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%od%netwk(irec))//'.'&
                      //trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)&
                      //trim(fpar%sim%RCOMPS(icomp))//'.sac'
          call sacio_newhead(header, real(fpar%sim%dt), fpar%sim%nstep, -real(T0))
          header%baz = this%od%baz(irec)
          header%dist = this%od%dist(irec)
          header%evla = fpar%acqui%evla(this%ievt)
          header%evlo = fpar%acqui%evlo(this%ievt)
          header%evdp = fpar%acqui%evdp(this%ievt)
          header%stla = this%od%stla(irec)
          header%stlo = this%od%stlo(irec)
          header%stel = this%od%stel(irec)
          header%knetwk = this%od%netwk(irec)
          header%kstnm = this%od%stnm(irec)
          header%kcmpnm = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))
          call sacio_writesac(datafile, header, this%data(:, icomp_syn, irec), ier)
          if (ier /= 0) call exit_MPI(0, 'Error writing SAC file '//trim(datafile))
        end do 
      end do
    end if
    call synchronize_all()

  end subroutine semd2sac

  subroutine preprocess(this, ievt)
    class(NoiseData), intent(inout) :: this
    integer, intent(in) :: ievt

    call this%init(ievt)

    call this%od%read_stations(ievt, .true.)

    call this%od%read_obs_data()

    call this%read(this%od%baz)
    
    call this%measure_adj()

  end subroutine preprocess

  subroutine measure_adj(this)
    use measure_adj_mod, only: measure_adj_fwat
    use ma_variables, only: OUT_DIR
    class(NoiseData), intent(inout) :: this
    integer :: irec_local, irec, iflt, icomp, icomp_syn, nb, ne, out_imeas
    real(kind=dp), dimension(:,:), allocatable :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn, adj_loc
    real(kind=dp), dimension(:,:,:), allocatable :: window_chi
    real(kind=dp), dimension(NCHI) :: window_chi_loc
    real(kind=dp), dimension(:,:,:,:), allocatable :: adj_src
    real(kind=dp), dimension(:), allocatable :: tstart, tend, seismo_tmp, seismo_dat,seismo_syn
    real(kind=dp), dimension(:), allocatable :: adj_2, adj_3, adj_1
    real(kind=dp), dimension(NDIM_MA) :: adj_syn_local
    real(kind=dp) :: dist_min, max_amp, tr_chi_loc, am_chi_loc, T_pmax_dat_loc, T_pmax_syn_loc
    character(len=MAX_STRING_LEN) :: file_prefix, msg
    character(len=MAX_STR_CHI), dimension(:), allocatable :: sta, net

    OUT_DIR = trim(OUTPUT_FILES)
    if (nrec_local > 0) adj_src = zeros_dp(NSTEP, fpar%sim%NRCOMP, nrec_local, fpar%sim%NUM_FILTER)
    do iflt = 1, fpar%sim%NUM_FILTER
      call get_band_name(fpar%sim%SHORT_P(iflt), fpar%sim%LONG_P(iflt), this%band_name)
      call this%wchi(iflt)%init(this%ievt, this%band_name)
      if (nrec_local > 0) then
        window_chi = zeros_dp(nrec_local, NCHI, fpar%sim%NRCOMP)
        tr_chi = zeros_dp(nrec_local, fpar%sim%NRCOMP)
        am_chi = zeros_dp(nrec_local, fpar%sim%NRCOMP)
        T_pmax_dat = zeros_dp(nrec_local, fpar%sim%NRCOMP)
        T_pmax_syn = zeros_dp(nrec_local, fpar%sim%NRCOMP)
        tstart = zeros_dp(nrec_local)
        tend = zeros_dp(nrec_local)
        if (iflt == 1) then
          allocate(sta(nrec_local))
          allocate(net(nrec_local))
        endif
        do irec_local = 1, nrec_local
          irec = number_receiver_global(irec_local)
          do icomp = 1, fpar%sim%NRCOMP
            ! get data component
            seismo_dat = interpolate_func_dp(this%od%data(:, icomp, irec), dble(this%od%tbeg(irec)), &
                                             dble(this%od%dt), this%od%npts, -dble(t0), dble(DT), NSTEP)
            if (.not. fpar%sim%SUPPRESS_EGF) call dif1(seismo_dat, dble(DT))
            call detrend(seismo_dat)
            call demean(seismo_dat)
            call bandpass_dp(seismo_dat, fpar%sim%nstep, dble(fpar%sim%dt),&
                             1/fpar%sim%LONG_P(iflt), 1/fpar%sim%SHORT_P(iflt), IORD)
            
            ! get synthetic component
            icomp_syn = get_icomp_syn(fpar%sim%RCOMPS(icomp))
            seismo_syn = this%data(:, icomp_syn, irec)
            call detrend(seismo_syn)
            call demean(seismo_syn)
            call bandpass_dp(seismo_syn, fpar%sim%nstep, dble(fpar%sim%dt),&
                             1/fpar%sim%LONG_P(iflt), 1/fpar%sim%SHORT_P(iflt), IORD)
            ! time window
            tstart(irec_local) = this%od%dist(irec)/fpar%sim%GROUPVEL_MAX(iflt)-fpar%sim%LONG_P(iflt)/2.
            tend(irec_local) = this%od%dist(irec)/fpar%sim%GROUPVEL_MIN(iflt)+fpar%sim%LONG_P(iflt)/2.
            tstart(irec_local) = max(tstart(irec_local), -t0)
            tend(irec_local) = min(tend(irec_local), NSTEP*dble(DT)-t0)

            ! reject near offsets
            if (.not. fpar%sim%USE_NEAR_OFFSET) then
              dist_min = (fpar%sim%GROUPVEL_MAX(iflt) + fpar%sim%GROUPVEL_MIN(iflt))*fpar%sim%LONG_P(iflt)/2.
              if (this%od%dist(irec) < dist_min) then
                seismo_syn = 1.0e-18_dp
                seismo_dat = 1.0e-18_dp
              end if
            end if

            ! normalization
            nb = floor((tstart(irec_local)+t0)/dble(DT)+1)
            ne = floor((tend(irec_local)+t0)/dble(DT)+1)
            seismo_dat = seismo_dat / maxval(abs(seismo_dat(nb:ne))) * maxval(abs(seismo_syn(nb:ne)))
            if (is_output_preproc) then
              call this%write_in_preocess(irec, icomp, tstart(irec_local), tend(irec_local), 'syn', seismo_syn)
              call this%write_in_preocess(irec, icomp, tstart(irec_local), tend(irec_local), 'dat', seismo_dat)
            end if
            call measure_adj_fwat(seismo_dat, seismo_syn, tstart(irec_local), tend(irec_local),&
                                  dble(fpar%sim%SHORT_P(iflt)), dble(fpar%sim%LONG_P(iflt)),&
                                  this%od%netwk(irec), this%od%stnm(irec),&
                                  trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp)), &
                                  window_chi_loc, tr_chi_loc,&
                                  am_chi_loc, T_pmax_dat_loc, &
                                  T_pmax_syn_loc, adj_syn_local, &
                                  file_prefix, out_imeas)
            window_chi(irec_local, :, icomp) = window_chi_loc
            adj_src(:, icomp, irec_local, iflt) = adj_syn_local(1:NSTEP)
            tr_chi(irec_local, icomp) = tr_chi_loc * fpar%acqui%src_weight(this%ievt)
            am_chi(irec_local, icomp) = am_chi_loc * fpar%acqui%src_weight(this%ievt)
            T_pmax_dat(irec_local, icomp) = T_pmax_dat_loc
            T_pmax_syn(irec_local, icomp) = T_pmax_syn_loc
            ! write adjoint source
            if (iflt == 1) then
              sta(irec_local) = this%od%stnm(irec)
              net(irec_local) = this%od%netwk(irec)
            endif
          end do
        end do
      end if
      call synchronize_all()
      call this%wchi(iflt)%assemble_window_chi(window_chi, tr_chi, am_chi,&
                                              T_pmax_dat, T_pmax_syn, sta, net,&
                                              tstart, tend)
      call this%wchi(iflt)%write()
      this%total_misfit(iflt) = this%wchi(iflt)%sum_chi(29)
      write(msg, '(a,f12.6)') 'Total misfit: of '//trim(this%band_name)//': ', this%total_misfit(iflt)
      call log%write(msg, .true.) 
    end do
    call synchronize_all()

    call this%get_comp_name_adj()
    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        do icomp = 1, fpar%sim%NRCOMP
          icomp_syn = get_icomp_syn(fpar%sim%RCOMPS(icomp))
          adj_loc = zeros_dp(NSTEP, 3)
          if (fpar%sim%ADJ_SRC_NORM) then
            max_amp = maxval(abs(adj_src(:, icomp, irec_local, :)))
            do iflt = 1, fpar%sim%NUM_FILTER
              adj_loc(:, icomp) = adj_loc(:, icomp) + adj_src(:, icomp, irec_local, iflt)&
                                  / maxval(abs(adj_src(:, icomp, irec_local, iflt))) * max_amp
            end do
          else
            adj_loc(:, icomp) = sum(adj_src(:, icomp, irec_local, :), dim=2)
          end if
        end do
        adj_1 = zeros_dp(NSTEP)
        adj_2 = zeros_dp(NSTEP)
        adj_3 = zeros_dp(NSTEP)
        call rotate_ZRT_to_ZNE(adj_loc(:, 1), adj_loc(:, 2), adj_loc(:, 3), &
                               adj_1, adj_2, adj_3, NSTEP, dble(this%od%baz(irec)))
        call this%write_adj(adj_1, trim(this%comp_name(1)), irec)
        call this%write_adj(adj_2, trim(this%comp_name(2)), irec)
        call this%write_adj(adj_3, trim(this%comp_name(3)), irec)
      end do
    end if
    call synchronize_all()
  end subroutine measure_adj

  subroutine finalize(this)
    class(NoiseData), intent(inout) :: this
    integer :: iflt

    call this%od%finalize()
    do iflt = 1, fpar%sim%NUM_FILTER
      call this%wchi(iflt)%finalize()
    end do
    call free_shm_array(this%dat_win)
  end subroutine finalize

  subroutine write_in_preocess(this, irec, icomp, tb, te, label, trace)
    class(NoiseData), intent(inout) :: this
    integer, intent(in) :: irec, icomp
    real(kind=dp), intent(in) :: tb, te
    real(kind=dp), dimension(:), intent(in) :: trace
    character(len=*), intent(in) :: label
    type(sachead) :: header
    character(len=MAX_STRING_LEN) :: datafile

    datafile = trim(OUTPUT_FILES)//'/'//trim(this%od%netwk(irec))//'.'//&
               trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)//&
               trim(fpar%sim%RCOMPS(icomp))//'.'//trim(label)//'.sac.'//trim(this%band_name)

    call sacio_newhead(header, real(DT), NSTEP, -real(T0))
    header%baz = this%od%baz(irec)
    header%dist = this%od%dist(irec)
    header%evla = fpar%acqui%evla(this%ievt)
    header%evlo = fpar%acqui%evlo(this%ievt)
    header%evdp = fpar%acqui%evdp(this%ievt)
    header%stla = this%od%stla(irec)
    header%stlo = this%od%stlo(irec)
    header%stel = this%od%stel(irec)
    header%knetwk = this%od%netwk(irec)
    header%kstnm = this%od%stnm(irec)
    header%kcmpnm = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))
    header%t1 = real(tb)
    header%t2 = real(te)
    call sacio_writesac(datafile, header, trace, ier)

  end subroutine write_in_preocess

  subroutine calc_distaz(this)
    class(NoiseData), intent(inout) :: this
    integer :: irec_local, irec
    real(kind=dp) :: dist, bazi, azi, delta

    if (worldrank == 0) then
      do irec = 1, this%od%nsta
        if (SUPPRESS_UTM_PROJECTION) then
          call distaz_cart(dble(fpar%acqui%evlo(this%ievt)), dble(fpar%acqui%evla(this%ievt)), &
                           dble(this%od%stlo(irec)), dble(this%od%stla(irec)), azi, bazi, dist)
          dist = dist / 1000.0_dp
        else
          call distaz(dble(this%od%stla(irec)), dble(this%od%stlo(irec)), &
                      dble(fpar%acqui%evla(this%ievt)), dble(fpar%acqui%evlo(this%ievt)),&
                      azi, bazi, delta, dist)
        endif
        this%od%baz(irec) = real(bazi)
        this%od%dist(irec) = real(dist)
      end do
    end if
    call synchronize_all()
    call sync_from_main_rank_cr_1d(this%od%baz, this%od%nsta)

  end subroutine calc_distaz

end module noise_data