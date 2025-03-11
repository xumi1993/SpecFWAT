module common_lib
  use config
  use fwat_constants
  implicit none

contains
  subroutine get_simu_type()
    if (dat_type == 'noise') then
      simu_type = SIMU_TYPE_NOISE
    else if (index(dat_type, 'tele') /= 0 .or. dat_type == 'rf') then
      simu_type = SIMU_TYPE_TELE
    else
      if (worldrank == 0) call exit_MPI(0, 'Unknown data type')
    endif
  end subroutine get_simu_type

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

    nt = size(vn)
    baz = deg2rad * bazi

    do it = 1, nt
      vr(it) = -ve(it) * sin(baz) - vn(it) * cos(baz)
      vt(it) = -ve(it) * cos(baz) + vn(it) * sin(baz)
    enddo
  
  end subroutine rotate_NE_to_RT

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

end module common_lib