module argparse
  use config
  use fwat_mpi

  implicit none

contains
  ! FK inherits simulation settings from FWAT; only the optional source-list index is selected here.
  subroutine parse_args_fk()
    integer :: iarg, argc, ios
    character(len=MAX_STRING_LEN) :: arg, value
    character(len=*), parameter :: usage = 'Usage: xfwat_fk [-e|--event <event_index>]'

    event_index = 0
    single_run = .false.
    argc = command_argument_count()
    iarg = 1
    do while (iarg <= argc)
      call get_command_argument(iarg, arg)
      select case (arg)
      case ('-h', '--help')
        if (worldrank == 0) print *, usage
        call finalize_mpi()
        stop
      case ('-e', '--event')
        if (single_run) call exit_MPI(worldrank, 'Repeated event option. '//usage)
        iarg = iarg + 1
        if (iarg > argc) call exit_MPI(worldrank, 'Missing event index. '//usage)
        call get_command_argument(iarg, value)
        if (len_trim(value) == 0 .or. verify(trim(value), '0123456789') /= 0) &
          call exit_MPI(worldrank, 'Event index must be a positive integer. '//usage)
        read(value, *, iostat=ios) event_index
        if (ios /= 0) call exit_MPI(worldrank, 'Invalid event index. '//usage)
        if (event_index < 1) call exit_MPI(worldrank, 'Event index must be positive. '//usage)
        single_run = .true.
      case default
        call exit_MPI(worldrank, 'Unknown option: '//trim(arg)//'. '//usage)
      end select
      iarg = iarg + 1
    enddo
  end subroutine parse_args_fk

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
    integer, parameter :: max_num_args = 5, min_num_args = 2
    character(len=MAX_STRING_LEN), dimension(max_num_args) :: argv
    integer :: i, iarg, argc
    character(len=MAX_STRING_LEN) :: usage

    usage = 'Usage: fwat_post_proc -m <model> [-r 1|2] [-h] -g'

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
          if (worldrank == 0) print *, trim(usage)
          call exit_MPI(0, 'ERROR: run-mode not set')
        endif
        read(argv(iarg), *) run_mode
      elseif (argv(i) == '-g' .or. argv(i) == '--use-gll') then
        use_gll = .true.
      endif
    enddo

  end subroutine parse_args_post_process

  subroutine parse_args_optimize()
    integer :: iarg, argc
    character(len=MAX_STRING_LEN) :: arg
    logical :: has_model
    character(len=*), parameter :: usage = 'Usage: fwat_optimize -m <model> [-g|--use-gll]'
    argc = command_argument_count()
    has_model = .false.
    use_gll = .false.
    iarg = 1
    do while (iarg <= argc)
      call get_command_argument(iarg, arg)
      select case (arg)
      case ('-m', '--model')
        if (has_model) call exit_MPI(worldrank, 'Repeated model option. '//usage)
        iarg = iarg + 1
        if (iarg > argc) call exit_MPI(worldrank, 'Model name not set. '//usage)
        call get_command_argument(iarg, model_name)
        if (len_trim(model_name) == 0 .or. model_name(1:1) == '-') &
          call exit_MPI(worldrank, 'Model name not set. '//usage)
        has_model = .true.
      case ('-g', '--use-gll')
        use_gll = .true.
      case ('-h', '--help')
        if (worldrank == 0) print *, usage
        call finalize_MPI()
        stop
      case default
        call exit_MPI(worldrank, 'Unknown option: '//trim(arg)//'. '//usage)
      end select
      iarg = iarg + 1
    enddo
    if (.not. has_model) call exit_MPI(worldrank, 'Model name not set. '//usage)
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
          if (worldrank == 0) print *, trim(usage)
          call exit_MPI(0, 'ERROR: simu_type not set')
        endif
        simu_type = argv(iarg)
      endif
    enddo

  end subroutine parse_args_mesh_databases

  ! Parse xspecfwat options, returning the first model index and iteration count.
  ! Also sets simu_type in config; MPI must be initialized for help/error exits.
  subroutine parse_invert_args(first, iteration_count)
    integer, intent(out) :: first, iteration_count
    integer :: iarg, argc, ios
    character(len=MAX_STRING_LEN) :: arg, value
    logical :: has_model, has_type, has_count
    character(len=*), parameter :: usage = 'Usage: xspecfwat -m M00 -s noise|tele|leq [-n iterations]'
    argc = command_argument_count()
    has_model = .false.
    has_type = .false.
    has_count = .false.
    iteration_count = 1
    first = -1
    iarg = 1
    do while (iarg <= argc)
      call get_command_argument(iarg, arg)
      if (arg == '-h' .or. arg == '--help') then
        if (worldrank == 0) print *, usage
        call finalize_mpi()
        stop
      endif
      if (iarg == argc) call exit_MPI(worldrank, 'Missing option value. '//usage)
      call get_command_argument(iarg+1, value)
      select case (arg)
      case ('-m', '--model')
        if (has_model) call exit_MPI(worldrank, 'Repeated model option')
        has_model = .true.
        if (len_trim(value) /= 3 .or. value(1:1) /= 'M' .or. verify(value(2:3), '0123456789') /= 0) &
          call exit_MPI(worldrank, 'Model must have the form M00 through M98')
        read(value(2:3), '(I2)', iostat=ios) first
        if (ios /= 0) call exit_MPI(worldrank, 'Invalid model index')
      case ('-s', '--simu-type')
        if (has_type) call exit_MPI(worldrank, 'Repeated simulation type')
        has_type = .true.
        simu_type = value
        if (.not. any(INV_TYPE_NAMES == simu_type)) call exit_MPI(worldrank, 'Unknown simulation type. '//usage)
      case ('-n', '--iterations')
        if (has_count) call exit_MPI(worldrank, 'Repeated iteration count')
        has_count = .true.
        if (len_trim(value) == 0 .or. verify(trim(value), '0123456789') /= 0) &
          call exit_MPI(worldrank, 'Iteration count must be a positive integer')
        read(value, *, iostat=ios) iteration_count
        if (ios /= 0) call exit_MPI(worldrank, 'Invalid iteration count')
      case default
        call exit_MPI(worldrank, 'Unknown option: '//trim(arg)//'. '//usage)
      end select
      iarg = iarg+2
    enddo
    if (.not. has_model .or. .not. has_type) call exit_MPI(worldrank, usage)
    if (first < 0 .or. first > 98 .or. iteration_count < 1 .or. iteration_count > 99-first) &
      call exit_MPI(worldrank, 'Requested iterations must produce models no later than M99')
  end subroutine parse_invert_args

end module argparse
