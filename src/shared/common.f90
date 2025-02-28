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
end module common_lib