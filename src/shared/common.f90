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
end module common_lib