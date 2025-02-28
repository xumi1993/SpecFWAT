module FKTimes_mod
  use fullwave_adjoint_tomo_par
  ! use specfem3D_par, only: al_FK, be_FK, mu_FK, h_FK, phi_FK, theta_FK, xx0,yy0,zz0,CUSTOM_REAL
  use specfem_par, only : CUSTOM_REAL, nrec, FKMODEL_FILE, SUPPRESS_UTM_PROJECTION, ILONGLAT2UTM, myrank
  use specfem_par_coupling
  use utils, only : zeros
  use fwat_input, only : acqui_par
  use constants, only : DEG2RAD

  implicit none

contains

  subroutine free_fk_arrays() 
    if (allocated(al_FK)) deallocate(al_FK)
    if (allocated(be_FK)) deallocate(be_FK)
    if (allocated(mu_FK)) deallocate(mu_FK)
    if (allocated(h_FK)) deallocate(h_FK)

  end subroutine free_fk_arrays

  subroutine ReadFKModelInput_fwat()

    real(kind=CUSTOM_REAL) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box
    integer                :: ioerr
    character(len=MAX_STRING_LEN)     :: keyword, keyvalue, line
    character(len=MAX_STRING_LEN)     :: keyword_tmp, incident_wave

    real(kind=CUSTOM_REAL) :: rho_layer, vp_layer, vs_layer, ztop_layer
    real(kind=CUSTOM_REAL) :: Radius_box, wave_length_at_bottom
    real(kind=CUSTOM_REAL), dimension(:), allocatable  :: rho_fk_input, ztop_fk_input, zbot_fk_input, vp_fk_input, vs_fk_input
    integer,  dimension(:), allocatable  :: ilayer_fk_input
    integer  :: ilayer,ier, fid
    logical  :: position_of_wavefront_not_read

    !!--------------------------------------------------------------
    ! # model description :
    ! NLAYER   n # number of layers
    ! LAYER 1  rho, vp ,vs, ztop
    ! LAYER 2  rho, vp, vs, ztop
    ! ...
    ! LAYER n  rho, vp, vs, ztop  # homogenoeus half-space
    ! #
    ! # incident wave description:
    ! INCIDENT_WAVE  "p" or "sv"
    ! BACK_AZITUTH    bazi
    ! INCIDENCE       inc
    ! ORIGIN_WAVEFRONT xx0, yy0, zz0
    ! ORIGIN_TIME      tt0
    ! FREQUENCY_MAX    ff0
    ! TIME_WINDOW      tmax_fk
    !!----------------------------------------------------------------

    ! only master process reads
    call FindBoundaryBox(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)

    if (myrank == 0) then

      !! set default values
      tt0 = 0.
      tmax_fk = 128.
      ff0 = 0.1
      kpsv = 1  ! 1 == P-wave / 2 == SV-wave
      position_of_wavefront_not_read = .true.
      stag = .false.

      !! READING input file
      open(newunit=fid,file=trim(FKMODEL_FILE), status="old", action="read")
      if (ioerr /= 0) then
        write(*,*) " ERROR READING FK INPUT FILE "
        write(*,*) " FILE NOT FOUND ", trim(FKMODEL_FILE)
        stop
      endif
      do
        read(fid, fmt='(a)',iostat=ioerr) line
        !call remove_begin_blanks(line)
        if (ioerr < 0) exit
        if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
        read(line,*) keyword, keyvalue
        select case(trim(keyword))
        case('NLAYER')
            read(line, *) keyword_tmp, nlayer
            ! allocate(al_FK(nlayer), be_FK(nlayer), mu_FK(nlayer), h_FK(nlayer),stat=ier)
            ! if (ier /= 0) call exit_MPI_without_rank('error allocating array 2226')
            ! allocate(rho_fk_input(nlayer),stat=ier)
            ! if (ier /= 0) call exit_MPI_without_rank('error allocating array 2227')
            ! allocate(vp_fk_input(nlayer),stat=ier)
            ! if (ier /= 0) call exit_MPI_without_rank('error allocating array 2228')
            ! allocate(vs_fk_input(nlayer),stat=ier)
            ! if (ier /= 0) call exit_MPI_without_rank('error allocating array 2229')
            ! allocate(ztop_fk_input(nlayer+1),stat=ier)
            ! if (ier /= 0) call exit_MPI_without_rank('error allocating array 2230')
            ! allocate(ilayer_fk_input(nlayer+1),stat=ier)
            ! if (ier /= 0) call exit_MPI_without_rank('error allocating array 2231')

            al_FK = zeros(nlayer)
            be_FK = zeros(nlayer)
            mu_FK = zeros(nlayer)
            h_FK = zeros(nlayer)
            rho_fk_input = zeros(nlayer)
            vp_fk_input = zeros(nlayer)
            vs_fk_input = zeros(nlayer)
            ztop_fk_input = zeros(nlayer)
            zbot_fk_input = zeros(nlayer)
            ilayer_fk_input = zeros(nlayer+1)
            ilayer_fk_input(:)=-1

        case('LAYER')
            read(line, *) keyword_tmp, ilayer, rho_layer, vp_layer, vs_layer, ztop_layer
            ilayer_fk_input(ilayer)=ilayer
            rho_fk_input(ilayer)=rho_layer
            vp_fk_input(ilayer)=vp_layer
            vs_fk_input(ilayer)=vs_layer
            ztop_fk_input(ilayer)=ztop_layer

        case('INCIDENT_WAVE')
            read(line,*)  keyword_tmp, incident_wave

            select case(trim(incident_wave))
                case ('p', 'P')
                  kpsv=1
                case('sv','SV')
                  kpsv=2
                case default
                  kpsv=1
                end select

        case('BACK_AZIMUTH')
            read(line,*)  keyword_tmp, phi_FK
            phi_FK = - phi_FK - 90.

        case('AZIMUTH')
            read(line,*)  keyword_tmp, phi_FK
            phi_FK = 90. - phi_FK

        case('TAKE_OFF')
            read(line,*)  keyword_tmp, theta_FK

        case('ORIGIN_WAVEFRONT')
            read(line,*)  keyword_tmp, xx0, yy0, zz0
            position_of_wavefront_not_read=.false.

          case('ORIGIN_TIME')
            read(line,*)  keyword_tmp, tt0

        case('FREQUENCY_MAX')
            read(line,*)  keyword_tmp, ff0

        case('TIME_WINDOW')
            read(line,*)  keyword_tmp, tmax_fk

        end select
      !!------------------------------------------------------------------------------------------------------
      enddo
      close(fid)

      if (allocated(ilayer_fk_input)) then

        ilayer_fk_input(nlayer+1) = ilayer_fk_input(nlayer)
        ! ztop_fk_input(nlayer+1)=ztop_fk_input(nlayer)
        Z_ref_for_FK=ztop_fk_input(nlayer)
        do ilayer=1, nlayer-1
            zbot_fk_input(ilayer) = ztop_fk_input(ilayer+1)
        enddo
        zbot_fk_input(nlayer) = Z_ref_for_FK
        do ilayer=1, nlayer
            al_FK(ilayer) = vp_fk_input(ilayer)
            be_FK(ilayer) = vs_fk_input(ilayer)
            mu_FK(ilayer) = rho_fk_input(ilayer) * vs_fk_input(ilayer)**2
            h_FK(ilayer) =  ztop_fk_input(ilayer) - zbot_fk_input(ilayer)

            if (ilayer_fk_input(ilayer) == -1) then
              write(*,*) " ERROR READING FK INPUT FILE "
              write(*,*) " MISSING LAYER ", ilayer
              stop
            endif
        enddo

        deallocate(ilayer_fk_input, rho_fk_input, ztop_fk_input, vp_fk_input, vs_fk_input)

      else

        write(*,*) " ERROR READING FK INPUT FILE "
        write(*,*) " NOT BE ABLE TO READ MODEL PROPERTIES "
        stop

      endif

      !! compute position of wave front
      if (position_of_wavefront_not_read) then
        xx0=0.5*(Xmin_box + Xmax_box)
        yy0=0.5*(Ymin_box + Ymax_box)
        Radius_box = sqrt( (Xmin_box - xx0)**2 + (Ymin_box - yy0)**2)

        if (kpsv == 1) then
            wave_length_at_bottom = al_FK(nlayer) / ff0
        else if (kpsv == 2) then
            wave_length_at_bottom = be_FK(nlayer) / ff0
        endif

        zz0 = Zmin_box - Radius_box * sin ( abs (theta_FK) * (acos(-1.d0) / 180.d0)  ) -  &
              3.*wave_length_at_bottom * cos ( abs (theta_FK) * (acos(-1.d0) / 180.d0)  ) -  &
              Z_ref_for_FK

      endif

      ! converts to rad
      phi_FK   = phi_FK * DEG2RAD
      theta_FK = theta_FK * DEG2RAD
    endif
    call synchronize_all()

    ! send FK parameters to others MPI slices
    call bcast_all_singlei(kpsv)
    call bcast_all_singlei(nlayer)

    if (myrank > 0) then
      allocate(al_FK(nlayer),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2206')
      allocate(be_FK(nlayer),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2207')
      allocate(mu_FK(nlayer),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2208')
      allocate(h_FK(nlayer),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2209')
    endif

    call bcast_all_cr(al_FK, nlayer)
    call bcast_all_cr(be_FK, nlayer)
    call bcast_all_cr(mu_FK, nlayer)
    call bcast_all_cr(h_FK, nlayer)
    call bcast_all_singlecr(phi_FK)
    call bcast_all_singlecr(theta_FK)
    call bcast_all_singlecr(ff0)
    call bcast_all_singlecr(xx0)
    call bcast_all_singlecr(yy0)
    call bcast_all_singlecr(zz0)
    call bcast_all_singlecr(tt0)
    call bcast_all_singlecr(Z_REF_for_FK)
    call bcast_all_singlecr(tmax_fk)
    call bcast_all_singlel(stag)

  end subroutine ReadFKModelInput_fwat

  subroutine traveltime(xx, yy, zz, tdelay)
    double precision, intent(in) :: xx, yy, zz
    double precision, intent(out) :: tdelay
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: h, v_fk_input
    double precision :: p, z0, zi, eta
    integer :: ilayer, j

    if (kpsv == 1) then
      v_fk_input = al_FK
      p = sin(theta_FK) / v_fk_input(nlayer)
    else if (kpsv == 2) then
      v_fk_input = be_FK
      p = sin(theta_FK) / v_fk_input(nlayer)
    endif
    tdelay = p * (xx - xx0) * cos(phi_FK) + p * (yy - yy0) * sin(phi_FK)

    z0 = zz0 - Z_ref_for_FK
    zi = zz - Z_ref_for_FK
    h = zeros(nlayer)
    ilayer = nlayer
    do j = nlayer - 1, 1, -1
      if (zi <= sum(h_FK(j:nlayer))) then
        ilayer = j
        exit
      end if
    end do
    h(ilayer+1:nlayer) = h_FK(ilayer+1:nlayer)
    h(ilayer) = zi - sum(h_FK(ilayer+1:nlayer))
    h(nlayer) = 0 - z0
    if (h(ilayer) < 0) then
      print *, 'Error setting layer thickness'
      stop
    end if
    do j = nlayer, ilayer, -1
      eta = sqrt(1 / v_fk_input(j)**2 - p**2)
      tdelay = tdelay + eta * h(j)
    end do

  end subroutine traveltime

  subroutine generate_fk_times(ievt, stla, stlo, stel, tdelay)

    integer, intent(in) :: ievt
    double precision, dimension(nrec), intent(in) :: stla, stlo, stel
    double precision, dimension(nrec), intent(out) :: tdelay
    double precision :: xx, yy, zz
    integer :: irec

    FKMODEL_FILE = 'src_rec/FKmodel_'//trim(acqui_par%evtid_names(ievt))

    call ReadFKModelInput_fwat()
    
    do irec = 1, nrec
      if (.not. SUPPRESS_UTM_PROJECTION) then
        call utm_geo(stlo(irec),stla(irec),xx,yy,ILONGLAT2UTM)
      else
        xx = stlo(irec)
        yy = stla(irec)
      endif
      zz = stel(irec)
      call traveltime(xx, yy, zz, tdelay(irec))
    end do
    call free_fk_arrays()

  end subroutine generate_fk_times

end module FKTimes_mod