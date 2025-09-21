module memory_monitor
  use iso_c_binding
  use fwat_constants
  use config, only: worldrank
  implicit none

  ! C接口声明用于getrusage
  interface
    function c_getrusage(who, usage) bind(C, name="getrusage")
      import :: c_int, c_ptr
      integer(c_int), value :: who
      type(c_ptr), value :: usage
      integer(c_int) :: c_getrusage
    end function c_getrusage
  end interface

  ! rusage结构体的Fortran等价物
  type, bind(C) :: rusage_type
    integer(c_long) :: ru_utime_sec      ! user CPU time used (seconds)
    integer(c_long) :: ru_utime_usec     ! user CPU time used (microseconds)  
    integer(c_long) :: ru_stime_sec      ! system CPU time used (seconds)
    integer(c_long) :: ru_stime_usec     ! system CPU time used (microseconds)
    integer(c_long) :: ru_maxrss         ! maximum resident set size (kB on Linux)
    integer(c_long) :: ru_ixrss          ! integral shared memory size
    integer(c_long) :: ru_idrss          ! integral unshared data size
    integer(c_long) :: ru_isrss          ! integral unshared stack size
    integer(c_long) :: ru_minflt         ! page reclaims (soft page faults)
    integer(c_long) :: ru_majflt         ! page faults (hard page faults)
    integer(c_long) :: ru_nswap          ! swaps
    integer(c_long) :: ru_inblock        ! block input operations
    integer(c_long) :: ru_oublock        ! block output operations
    integer(c_long) :: ru_msgsnd         ! IPC messages sent
    integer(c_long) :: ru_msgrcv         ! IPC messages received
    integer(c_long) :: ru_nsignals       ! signals received
    integer(c_long) :: ru_nvcsw          ! voluntary context switches
    integer(c_long) :: ru_nivcsw         ! involuntary context switches
  end type rusage_type

  ! getrusage的who参数常量
  integer(c_int), parameter :: RUSAGE_SELF = 0
  integer(c_int), parameter :: RUSAGE_CHILDREN = -1

contains

  subroutine print_memory_usage(prefix)
    character(len=*), intent(in), optional :: prefix
    character(len=256) :: mem_prefix
    type(rusage_type), target :: usage
    integer(c_int) :: ret
    real(kind=dp) :: rss_mb, utime, stime

    ! 只有rank 0输出内存信息
    if (worldrank /= 0) return

    ! 设置前缀
    if (present(prefix)) then
      mem_prefix = trim(prefix) // ' - '
    else
      mem_prefix = ''
    endif

    ! 调用getrusage获取资源使用情况
    ret = c_getrusage(RUSAGE_SELF, c_loc(usage))
    
    if (ret /= 0) then
      print '(a)', trim(mem_prefix) // 'Warning: getrusage() failed'
      return
    endif

    ! 转换为MB (Linux上ru_maxrss以kB为单位)
    rss_mb = real(usage%ru_maxrss, kind=dp) / 1024.0_dp
    
    ! 计算CPU时间 (秒)
    utime = real(usage%ru_utime_sec, kind=dp) + real(usage%ru_utime_usec, kind=dp) / 1000000.0_dp
    stime = real(usage%ru_stime_sec, kind=dp) + real(usage%ru_stime_usec, kind=dp) / 1000000.0_dp

    ! 输出内存和CPU使用信息
    print '(a,f8.2,a,i0,a)', trim(mem_prefix) // 'Max RSS Memory: ', &
          rss_mb, ' MB (', usage%ru_maxrss, ' kB)'
    print '(a,f8.3,a,f8.3,a)', trim(mem_prefix) // 'CPU Time - User: ', &
          utime, ' s, System: ', stime, ' s'
    print '(a,i0,a,i0)', trim(mem_prefix) // 'Page Faults - Minor: ', &
          usage%ru_minflt, ', Major: ', usage%ru_majflt

  end subroutine print_memory_usage

  function get_max_rss_mb() result(rss_mb)
    real(kind=dp) :: rss_mb
    type(rusage_type), target :: usage
    integer(c_int) :: ret

    rss_mb = -1.0_dp  ! 错误时返回-1

    ret = c_getrusage(RUSAGE_SELF, c_loc(usage))
    if (ret == 0) then
      rss_mb = real(usage%ru_maxrss, kind=dp) / 1024.0_dp
    endif

  end function get_max_rss_mb

  function get_max_rss_kb() result(rss_kb)
    integer(c_long) :: rss_kb
    type(rusage_type), target :: usage
    integer(c_int) :: ret

    rss_kb = -1  ! 错误时返回-1

    ret = c_getrusage(RUSAGE_SELF, c_loc(usage))
    if (ret == 0) then
      rss_kb = usage%ru_maxrss
    endif

  end function get_max_rss_kb

  subroutine print_memory_summary(prefix, show_cpu, show_pagefaults)
    character(len=*), intent(in), optional :: prefix
    logical, intent(in), optional :: show_cpu, show_pagefaults
    character(len=256) :: mem_prefix
    type(rusage_type), target :: usage
    integer(c_int) :: ret
    real(kind=dp) :: rss_mb, utime, stime
    logical :: do_cpu, do_pf

    ! 只有rank 0输出内存信息
    if (worldrank /= 0) return

    ! 设置选项
    do_cpu = .true.
    do_pf = .false.
    if (present(show_cpu)) do_cpu = show_cpu
    if (present(show_pagefaults)) do_pf = show_pagefaults

    ! 设置前缀
    if (present(prefix)) then
      mem_prefix = trim(prefix) // ' - '
    else
      mem_prefix = ''
    endif

    ! 调用getrusage获取资源使用情况
    ret = c_getrusage(RUSAGE_SELF, c_loc(usage))
    
    if (ret /= 0) then
      print '(a)', trim(mem_prefix) // 'Warning: getrusage() failed'
      return
    endif

    ! 转换为MB
    rss_mb = real(usage%ru_maxrss, kind=dp) / 1024.0_dp
    
    ! 输出内存信息
    print '(a,f8.2,a)', trim(mem_prefix) // 'Max RSS Memory: ', rss_mb, ' MB'

    ! 可选输出CPU时间
    if (do_cpu) then
      utime = real(usage%ru_utime_sec, kind=dp) + real(usage%ru_utime_usec, kind=dp) / 1000000.0_dp
      stime = real(usage%ru_stime_sec, kind=dp) + real(usage%ru_stime_usec, kind=dp) / 1000000.0_dp
      print '(a,f8.3,a,f8.3,a)', trim(mem_prefix) // 'CPU Time - User: ', &
            utime, ' s, System: ', stime, ' s'
    endif

    ! 可选输出页面错误信息
    if (do_pf) then
      print '(a,i0,a,i0)', trim(mem_prefix) // 'Page Faults - Minor: ', &
            usage%ru_minflt, ', Major: ', usage%ru_majflt
    endif

  end subroutine print_memory_summary

  subroutine get_resource_usage(rss_mb, user_time, system_time, minor_faults, major_faults)
    real(kind=dp), intent(out), optional :: rss_mb, user_time, system_time
    integer(c_long), intent(out), optional :: minor_faults, major_faults
    type(rusage_type), target :: usage
    integer(c_int) :: ret

    ! 调用getrusage获取资源使用情况
    ret = c_getrusage(RUSAGE_SELF, c_loc(usage))
    
    if (ret /= 0) then
      ! 出错时设置默认值
      if (present(rss_mb)) rss_mb = -1.0_dp
      if (present(user_time)) user_time = -1.0_dp
      if (present(system_time)) system_time = -1.0_dp
      if (present(minor_faults)) minor_faults = -1
      if (present(major_faults)) major_faults = -1
      return
    endif

    ! 返回请求的值
    if (present(rss_mb)) then
      rss_mb = real(usage%ru_maxrss, kind=dp) / 1024.0_dp
    endif
    
    if (present(user_time)) then
      user_time = real(usage%ru_utime_sec, kind=dp) + real(usage%ru_utime_usec, kind=dp) / 1000000.0_dp
    endif
    
    if (present(system_time)) then
      system_time = real(usage%ru_stime_sec, kind=dp) + real(usage%ru_stime_usec, kind=dp) / 1000000.0_dp
    endif
    
    if (present(minor_faults)) then
      minor_faults = usage%ru_minflt
    endif
    
    if (present(major_faults)) then
      major_faults = usage%ru_majflt
    endif

  end subroutine get_resource_usage

end module memory_monitor