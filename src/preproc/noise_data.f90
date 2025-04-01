module noise_data
  use config
  use fwat_constants
  use syn_data, only: SynData
  use input_params, fpar => fwat_par_global
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, nrec, nrec_local,OUTPUT_FILES, &
                        number_receiver_global, ispec_selected_rec, islice_selected_rec
  use fwat_mpi
  use utils, only: zeros_dp, zeros
  use distaz_lib

  implicit none
  type, extends(SynData) :: NoiseData
  contains
    procedure :: semd2sac
    procedure, private :: calc_distaz
  end type NoiseData

contains
  subroutine semd2sac(this, ievt)
    class(NoiseData), intent(in) :: this
    integer, intent(in) :: ievt
    integer :: irec_local, irec, icomp
    character(len=MAX_STRING_LEN) :: datafile

    call this%init(ievt)

    call this%od%read_stations(ievt, .true.)

    call this%calc_distaz()

    call call this%read(this%od%baz)

    if (nrec_local > 0) then
      do irec_local = 1, nrec_local
        irec = number_receiver_global(irec_local)
        do icomp = 1, fpar%sim%NRCOMP
          datafile = trim(fpar%acqui%in_dat_path(this%ievt))//'/'//trim(this%od%netwk(irec))//'.'&
          //trim(this%od%stnm(irec))//'.'//trim(fpar%sim%CH_CODE)&
          //trim(fpar%sim%RCOMPS(icomp))//'.sac'
          call sacio_newhead(header, real(fpar%sim%dt), fpar%sim%nstep, -real(T0))
          header%baz = this%od%baz(irec)
          header%evla = fpar%acqui%evla(this%ievt)
          header%evlo = fpar%acqui%evlo(this%ievt)
          header%evdp = fpar%acqui%evdp(this%ievt)
          header%stla = this%od%stla(irec)
          header%stlo = this%od%stlo(irec)
          header%stel = this%od%stel(irec)
          header%knetwk = this%od%netwk(irec)
          header%kstnm = this%od%stnm(irec)
          header%kcmpnm = trim(fpar%sim%CH_CODE)//trim(fpar%sim%RCOMPS(icomp))
          call sacio_writesac(datafile, header, this%data(:, icomp, irec), ier)
          if (ier /= 0) call exit_MPI(0, 'Error writing SAC file '//trim(datafile))
        end do 
      end do
    end if
    call synchronize_all()

  end subroutine semd2sac

  subroutine calc_distaz(this)
    class(NoiseData), intent(inout) :: this
    integer :: irec_local, irec
    real(kind=dp) :: dist, bazi, azi, delta

    if (worldrank == 0) then
      do irec = 1, this%od%nsta
        if (SUPPRESS_UTM_PROJECTION) then
          call distaz_cart(fpar%acqui%evlo(this%ievt), fpar%acqui%evla(this%ievt), &
                          this%od%stlo(irec), this%od%stla(irec), bazi, dist)
        else
          call distaz(fpar%acqui%evla(this%ievt), fpar%acqui%evlo(this%ievt), &
                      this%od%stla(irec), this%od%stlo(irec), azi, bazi, delta, dist)
        endif
        this%od%baz(irec) = real(bazi)
      end do
    end if
    call synchronize_all()
    call sync_from_main_rank_cr_1d(this%od%baz, this%od%nsta)

  end subroutine calc_distaz

end module noise_data