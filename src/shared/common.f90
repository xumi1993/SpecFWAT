module common_lib
  use config
  use fwat_constants
  implicit none

contains
  subroutine get_simu_type()
    select case (dat_type)
      case (SIMU_TYPE_NOISE)
        simu_type = SIMU_TYPE_NOISE
      case (SIMU_TYPE_LEQ)
        simu_type = SIMU_TYPE_LEQ
      case (SIMU_TYPE_TELE)
        simu_type = SIMU_TYPE_TELE
      case ('rf')
        simu_type = SIMU_TYPE_TELE
      case ('telecc')
        simu_type = SIMU_TYPE_TELE
      case default
        if (worldrank == 0) call exit_MPI(0, 'Unknown data type')
    end select
  end subroutine get_simu_type

  function read_iter_num()
    integer :: read_iter_num
    
    read(model_name(2:), '(I2.2)') read_iter_num
  end function read_iter_num

  subroutine get_dat_type()
    use input_params, fpar => fwat_par_global
    select case (simu_type)
    case (SIMU_TYPE_NOISE)
      dat_type = SIMU_TYPE_NOISE
    case (SIMU_TYPE_LEQ)
      dat_type = SIMU_TYPE_LEQ
    case (SIMU_TYPE_TELE)
      call get_tele_type(fpar%sim%TELE_TYPE)
    case default
      if (worldrank == 0) call exit_MPI(0, 'Unknown simulation type')
    end select
  end subroutine get_dat_type

  subroutine get_tele_type(tele_type)
    use input_params, fpar => fwat_par_global
    integer, intent(in) :: tele_type
    select case (tele_type)
      case (1)
        dat_type = 'tele'
      case (2)
        dat_type = 'rf'
      case (3)
        dat_type = 'telecc'
      case default
        if (worldrank == 0) call exit_MPI(0, 'Unknown teleseismic type')
    end select
  end subroutine get_tele_type

  subroutine get_kernel_names()
    use input_params, fpar => fwat_par_global
      if (fpar%update%model_type == 1) then
      nkernel = size(KERNEL_ISO)
      kernel_names = KERNEL_ISO
      parameter_names = MODEL_ISO
    elseif (fpar%update%model_type == 2) then
      nkernel = size(KERNEL_AZI_ANI)
      kernel_names = KERNEL_AZI_ANI
      parameter_names = MODEL_AZI_ANI
    else
      call exit_MPI(0, 'Unknown model type')
    endif
  end subroutine get_kernel_names

  integer function find_string(string_list, search_str)
    character(len=MAX_STRING_LEN) :: search_str
    character(len=MAX_STRING_LEN), dimension(:) :: string_list
    integer :: left, right, mid

    left = 1
    right = size(string_list)

    find_string = -1

    do while (left <= right)
      mid = (left + right) / 2
      if (trim(adjustl(string_list(mid))) == trim(adjustl(search_str))) then
        find_string = mid
        return
      elseif (trim(adjustl(string_list(mid))) < trim(adjustl(search_str))) then
        left = mid + 1
      else
        right = mid - 1
      end if
    end do
  end function find_string

  ! Rotation of components
  subroutine rotate_NE_to_RT(vn,ve,vr,vt,bazi)
  
    real(kind=cr),   intent(in) :: bazi
    real(kind=cr), dimension(:),  intent(in) :: vn, ve
    real(kind=cr), dimension(:), intent(out) :: vr, vt
    real(kind=cr) :: baz
    integer :: it, nt

    baz = deg2rad * bazi

    do it = 1, nt
      vr(it) = -ve(it) * sin(baz) - vn(it) * cos(baz)
      vt(it) = -ve(it) * cos(baz) + vn(it) * sin(baz)
    enddo
  
  end subroutine rotate_NE_to_RT

  ! subroutine rotate_RT_to_NE(vr,vt,vn,ve,bazi)

  !   real(kind=CUSTOM_REAL),   intent(in) :: bazi

  !   real(kind=CUSTOM_REAL), dimension(:),  intent(in) :: vz2, vr, vt
  !   real(kind=CUSTOM_REAL), dimension(:), intent(out) :: vz,  vn, ve

  !   real(kind=CUSTOM_REAL) :: baz
  !   integer :: it, nt
    
  !   nt = size(vn)
  !   baz = deg2rad * bazi

  !   do it = 1, nt
  !     ve(it) = -vr(it) * sin(baz) - vt(it) * cos(baz)
  !     vn(it) = -vr(it) * cos(baz) + vt(it) * sin(baz)
  !     vz(it) = vz2(it)
  !   enddo

  ! end subroutine rotate_ZRT_to_ZNE

  subroutine rotate_NE_to_RT_dp(vn,ve,vr,vt,bazi)
  
    real(kind=cr),   intent(in) :: bazi
    real(kind=dp), dimension(:),  intent(in) :: vn, ve
    real(kind=dp), dimension(:), intent(out) :: vr, vt
    real(kind=cr) :: baz
    integer :: it, nt

    nt = size(vn)
    baz = deg2rad * bazi

    do it = 1, nt
      vr(it) = -ve(it) * sin(baz) - vn(it) * cos(baz)
      vt(it) = -ve(it) * cos(baz) + vn(it) * sin(baz)
    enddo
  
  end subroutine rotate_NE_to_RT_dp

  subroutine rotate_R_to_NE(vr, vn, ve, bazi)
    real(kind=cr), dimension(:), intent(in) :: vr
    real(kind=cr), intent(in) :: bazi
    real(kind=cr), dimension(:), intent(out) :: ve, vn
    real(kind=cr), dimension(:), allocatable :: vt
    real(kind=cr) :: baz
    integer :: nt

    nt = size(vr)
    allocate(vt(nt))
    vt = 0.0_cr
    baz = 360.0 - bazi

    call rotate_NE_to_RT(vr, vt, ve, vn, baz)  
    
  end subroutine rotate_R_to_NE

  subroutine rotate_R_to_NE_dp(vr, vn, ve, bazi)
    real(kind=dp), dimension(:), intent(in) :: vr
    real(kind=cr), intent(in) :: bazi
    real(kind=dp), dimension(:), intent(out) :: ve, vn
    real(kind=dp), dimension(:), allocatable :: vt
    real(kind=cr) :: baz
    integer :: nt

    nt = size(vr)
    allocate(vt(nt))
    vt = 0.0_dp
    baz = 360.0 - bazi

    call rotate_NE_to_RT_dp(vr, vt, ve, vn, baz)  
    
  end subroutine rotate_R_to_NE_dp

  subroutine rotate_T_to_NE_dp(vt, vn, ve, bazi)
    real(kind=dp), dimension(:), intent(in) :: vt
    real(kind=cr), intent(in) :: bazi
    real(kind=dp), dimension(:), intent(out) :: ve, vn
    real(kind=dp), dimension(:), allocatable :: vr
    real(kind=cr) :: baz
    integer :: nt

    nt = size(vt)
    allocate(vr(nt))
    vr = 0.0_dp
    baz = 360.0 - bazi

    call rotate_NE_to_RT_dp(vr, vt, ve, vn, baz)  
    
  end subroutine rotate_T_to_NE_dp

  subroutine get_band_name(SHORT_P, LONG_P, bandname)
    real :: SHORT_P, LONG_P, fl, fh
    character(len=MAX_STRING_LEN), intent(out) :: bandname
    character(len=MAX_STRING_LEN) :: bandstr1, bandstr2

    fh = 1./SHORT_P
    fl = 1./LONG_P
    if (SHORT_P < 1) then
      write(bandstr1, '(a1,i3.3)') 'F',int(fh)
    else
      write(bandstr1, '(a1,i3.3)') 'T',int(SHORT_P)
    endif
    if (LONG_P < 1) then
      write(bandstr2, '(a1,i3.3)') 'F',int(fl)
    else
      write(bandstr2, '(a1,i3.3)') 'T',int(LONG_P)
    endif
    bandname = trim(bandstr1)//'_'//trim(bandstr2)

  end subroutine get_band_name

  subroutine dwascii(fname,data,npt,b,dt)

    character(len=*),intent(in) :: fname
    integer, intent(in) :: npt
    real(kind=dp), intent(in) :: data(*)
    real(kind=dp), intent(in) :: b,dt
    integer :: ios,i

    open(9,file=trim(fname),iostat=ios,status='unknown')
    if (ios /= 0) then
       print *,'Error opening ascii file to write: ',trim(fname)
       call exit_MPI(worldrank, 'Error opening ascii file to write')
    endif
    do i = 1,npt
       write(9,*) b + dt*(i-1), data(i)
    enddo
    close(9)

  end subroutine dwascii

  function get_gauss_fac(freq_max) result(gauss_fac)
    real(kind=cr), intent(in) :: freq_max
    real(kind=cr) :: gauss_fac

    gauss_fac =  pi * freq_max / sqrt(-log(0.1))
    ! if (gauss_fac < 1.5) then
      ! gauss_fac = 1.5
    ! endif
  end function get_gauss_fac

  function get_icomp_syn(comp_code) result(icomp)
    character(len=1), intent(in) :: comp_code
    integer :: icomp

    if (comp_code == 'Z') then
      icomp = 1
    elseif (comp_code == 'R') then
      icomp = 2
    elseif (comp_code == 'T') then
      icomp = 3
    else
      if (worldrank == 0) call exit_MPI(0, 'Unknown component')
    endif

  end function get_icomp_syn

  subroutine rotate_ZRT_to_ZNE(vz2,vr,vt,vz,vn,ve,nt,bazi)
  
    integer,                  intent(in) :: nt
    real(kind=dp),   intent(in) :: bazi

    real(kind=dp), dimension(nt),  intent(in) :: vz2, vr, vt
    real(kind=dp), dimension(nt), intent(out) :: vz,  vn, ve

    real(kind=dp) :: baz
    integer :: it

    baz = deg2rad * bazi

    do it = 1, nt
        ve(it) = -vr(it) * sin(baz) - vt(it) * cos(baz)
        vn(it) = -vr(it) * cos(baz) + vt(it) * sin(baz)
        vz(it) = vz2(it)
    enddo
  
  end subroutine rotate_ZRT_to_ZNE

  subroutine mkdir(dirname)
    character(len=*), intent(in) :: dirname
    integer :: ios

    if (worldrank == 0) call execute_command_line('mkdir -p ' // trim(dirname), wait=.true., exitstat=ios)
    call bcast_all_singlei(ios)
    if (ios /= 0) then
      call exit_MPI(0, 'Error creating directory')
    endif
  end subroutine mkdir

end module common_lib