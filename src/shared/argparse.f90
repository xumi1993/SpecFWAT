module argparse
  use config
  use fwat_mpi

  implicit none

  integer, private :: ier

contains
  subroutine parse_args_fwd_meas_adj(ievt)
    integer, parameter :: max_num_args = 8
    character(len=MAX_STRING_LEN), dimension(max_num_args) :: argv
    integer, intent(out) :: ievt
    integer :: i, iarg, argc
    character(len=MAX_STRING_LEN) :: arg, usage

    usage = 'Usage: fwat_fwd_measure_adj -m <model> -s <simu_type> -r <run_mode> [-e <event_index>]'

    argc = command_argument_count()
    do i = 1, argc
      call get_command_argument(i, argv(i))
    enddo
    if (argc > max_num_args .or. argc < max_num_args - 2) then
      if (worldrank == 0) print *, trim(usage)
      call exit_MPI(0, 'ERROR: Too more or too less arguments')
    endif

    ! parse arguments
    run_mode = 0
    ievt = 0

    do i = 1, argc
      arg = argv(i)
      if (arg == '-m' .or. arg == '--model') then
        iarg = i + 1
        if (iarg > argc) then
          if (worldrank == 0) print *, trim(usage)
          call exit_MPI(0, 'ERROR: Model name not set')
        endif
        model_name = argv(iarg)
      elseif (arg == '-s' .or. arg == '--simu-type') then
        iarg = i + 1
        if (iarg > argc) then
          if (worldrank == 0) print *, trim(usage)
          call exit_MPI(0, 'ERROR: data-type not set')
        endif
        simu_type = argv(iarg)
      elseif (arg == '-h' .or. arg == '--help') then
        if (worldrank == 0) print *, trim(usage)
        call finalize_MPI()
        call exit(0)
      elseif (arg == '-r' .or. arg == '--run-mode') then
        iarg = i + 1
        if (iarg > argc) then
          if (worldrank == 0) print *, trim(usage)
          call exit_MPI(0, 'ERROR: run-mode not set')
        endif
        read(argv(iarg), *) run_mode
      elseif (arg == '-e' .or. arg == '--event') then
        iarg = i + 1
        if (iarg > argc) then
          if (worldrank == 0) print *, trim(usage)
          call exit_MPI(0, 'ERROR: event index not set')
        endif
        read(argv(iarg), *) ievt
        single_run = .true.
      endif
    enddo

    if (run_mode == 0 .or. len_trim(model_name) == 0 .or. len_trim(simu_type) == 0) then
      if (worldrank == 0) print *, trim(usage)
      call exit_MPI(0, 'ERROR: Invalid arguments')
    endif
  end subroutine parse_args_fwd_meas_adj

  subroutine parse_args_post_process()
    integer, parameter :: max_num_args = 4, min_num_args = 2
    character(len=MAX_STRING_LEN), dimension(max_num_args) :: argv
    integer :: i, iarg, argc
    character(len=MAX_STRING_LEN) :: usage

    usage = 'Usage: fwat_post_proc -m <model> [-r 1|2] [-h]'

    argc = command_argument_count()
    do i = 1, argc
      call get_command_argument(i, argv(i))
    enddo

    if (argc /= min_num_args .and. argc /= max_num_args) then
      if (worldrank == 0) print *, trim(usage)
      call exit_MPI(0, 'ERROR: Too more or too less arguments')
    endif

    run_mode = 1  ! default run mode
    ! parse arguments
    do i = 1, argc
      if (argv(i) == '-m' .or. argv(i) == '--model') then
        iarg = i + 1
        if (iarg > argc) then
          if (worldrank == 0) print *, usage
          call exit_MPI(0, 'ERROR: Model name not set')
        endif
        model_name = argv(iarg)
      elseif (argv(i) == '-h' .or. argv(i) == '--help') then
        if (worldrank == 0) print *, trim(usage)
        call finalize_MPI()
        stop
      elseif (argv(i) == '-r' .or. argv(i) == '--run-mode') then
        iarg = i + 1
        if (iarg > argc) then
          if (worldrank == 0) print *, usage
          call exit_MPI(0, 'ERROR: run-mode not set')
        endif
        read(argv(iarg), *) run_mode
      endif
    enddo

  end subroutine parse_args_post_process

  subroutine parse_args_optimize()
    integer, parameter :: max_num_args = 2
    character(len=MAX_STRING_LEN), dimension(max_num_args) :: argv
    integer :: i, iarg, argc
    character(len=MAX_STRING_LEN) :: usage

    usage = 'Usage: fwat_optimize -m <model>'

    argc = command_argument_count()
    do i = 1, argc
      call get_command_argument(i, argv(i))
    enddo

    if (argc /= max_num_args) then
      if (worldrank == 0) print *, trim(usage)
      call exit_MPI(0, 'ERROR: Too more arguments')
    endif

    ! parse arguments
    do i = 1, argc
      if (argv(i) == '-m' .or. argv(i) == '--model') then
        iarg = i + 1
        if (iarg > argc) then
          if (worldrank == 0) print *, usage
          call exit_MPI(0, 'ERROR: Model name not set')
        endif
        model_name = argv(iarg)
      endif
    enddo
  end subroutine parse_args_optimize

  subroutine parse_args_mesh_databases()
    integer, parameter :: max_num_args = 2
    character(len=MAX_STRING_LEN), dimension(max_num_args) :: argv
    integer :: i, iarg, argc
    character(len=MAX_STRING_LEN) :: usage

    usage = 'Usage: fwat_mesh_databases -s <simu_type>'

    argc = command_argument_count()
    do i = 1, argc
      call get_command_argument(i, argv(i))
    enddo

    if (argc /= max_num_args) then
      if (worldrank == 0) print *, trim(usage)
      call exit_MPI(0, 'ERROR: Too more arguments')
    endif

    ! parse arguments
    do i = 1, argc
      if (argv(i) == '-s' .or. argv(i) == '--simu_type') then
        iarg = i + 1
        if (iarg > argc) then
          if (worldrank == 0) print *, usage
          call exit_MPI(0, 'ERROR: simu_type not set')
        endif
        simu_type = argv(iarg)
      endif
    enddo

  end subroutine parse_args_mesh_databases
end module argparse