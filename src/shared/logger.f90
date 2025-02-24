module logger
  use config

  implicit none
  integer, parameter :: log_unit=888

  type :: logger_type
    character(len=MAX_STRING_LEN) :: log_file
    contains
    procedure :: init => logger_init, write => logger_write, finalize => logger_finalize
  end type logger_type

contains

  subroutine logger_init(this, log_file)
    class(logger_type), intent(inout) :: this
    character(len=MAX_STRING_LEN), intent(in) :: log_file

    if(worldrank == 0) open(unit=log_unit, file=log_file, status='replace', action='write')

  end subroutine logger_init

  subroutine logger_write(this, message, is_time)
    class(logger_type), intent(inout) :: this
    character(len=MAX_STRING_LEN), intent(in) :: message
    logical, optional, intent(in) :: is_time
    logical :: is_time_loc

    if (present(is_time)) then
      is_time_loc = is_time
    else
      is_time_loc = .false.
    end if

    if(worldrank == 0) then
      if (is_time_loc) then
        write(log_unit, '("["A"] ", A)') nowtime(), message
      else
        write(log_unit, '(A)') message
      end if
      call flush(log_unit)
    endif
  end subroutine logger_write

  subroutine logger_finalize(this)
    class(logger_type), intent(inout) :: this

    if(worldrank == 0) close(log_unit)

  end subroutine logger_finalize

  function nowtime() result(timestamp)
    character(len=MAX_STRING_LEN) :: timestamp
    integer, dimension(8) :: time_values

    call date_and_time(VALUES=time_values)
    write(timestamp,'(i4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2)') time_values(1),time_values(2), &
          time_values(3),time_values(5),time_values(6),time_values(7)
  end function

end module logger