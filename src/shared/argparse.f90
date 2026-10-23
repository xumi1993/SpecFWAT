module argparse
  use config
  use fwat_mpi

  implicit none

  private
  public :: parse_args_fk, parse_args_fwd_meas_adj, parse_args_post_process, &
            parse_args_optimize, parse_args_mesh_databases, parse_invert_args

contains
  ! FK inherits simulation settings from FWAT; only the optional source-list index is selected here.
  subroutine parse_args_fk()
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: argv
    character(len=MAX_STRING_LEN) :: usage
    integer :: i, argc

    usage = 'Usage: xfwat_fk [-e <event_index>] [-h]'

    call read_args(usage, argv, argc)

    single_run = .false.
    event_index = 0

    i = 1
    do while (i <= argc)
      select case (trim(argv(i)))
      case ('-e', '--event')
        call get_int_value(argv, i, usage, event_index, single_run)
      case default
        call arg_error(usage, 'ERROR: Unknown argument: '//trim(argv(i)))
      end select
      i = i + 1
    enddo

    if (single_run .and. event_index < 1) call arg_error(usage, 'ERROR: event index must be >= 1')
  end subroutine parse_args_fk

  subroutine parse_args_fwd_meas_adj(ievt)
    integer, intent(out) :: ievt
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: argv
    character(len=MAX_STRING_LEN) :: usage
    integer :: i, argc
    logical :: has_model, has_simu_type, has_run_mode

    usage = 'Usage: xfwat_fwd_measure_adj -m <model> -s <simu_type> -r <run_mode 1-5> [-e <event_index>] [-h]'

    call read_args(usage, argv, argc)

    has_model = .false.
    has_simu_type = .false.
    has_run_mode = .false.
    single_run = .false.
    ievt = 0

    i = 1
    do while (i <= argc)
      select case (trim(argv(i)))
      case ('-m', '--model')
        call get_str_value(argv, i, usage, model_name, has_model)
      case ('-s', '--simu-type', '--simu_type')
        call get_str_value(argv, i, usage, simu_type, has_simu_type)
      case ('-r', '--run-mode')
        call get_int_value(argv, i, usage, run_mode, has_run_mode)
      case ('-e', '--event')
        call get_int_value(argv, i, usage, ievt, single_run)
      case default
        call arg_error(usage, 'ERROR: Unknown argument: '//trim(argv(i)))
      end select
      i = i + 1
    enddo

    if (.not. has_model) call arg_error(usage, 'ERROR: Model name not set')
    if (.not. has_simu_type) call arg_error(usage, 'ERROR: simu-type not set')
    if (.not. has_run_mode) call arg_error(usage, 'ERROR: run-mode not set')
    if (run_mode < FORWARD_ONLY .or. run_mode > ADJOINT_ONLY) &
      call arg_error(usage, 'ERROR: run-mode must be in 1-5')
    if (single_run .and. ievt < 1) call arg_error(usage, 'ERROR: event index must be >= 1')
  end subroutine parse_args_fwd_meas_adj

  subroutine parse_args_post_process()
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: argv
    character(len=MAX_STRING_LEN) :: usage
    integer :: i, argc
    logical :: has_model, has_run_mode

    usage = 'Usage: xfwat_post_proc -m <model> [-r 1|2] [-g] [-h]'

    call read_args(usage, argv, argc)

    has_model = .false.
    has_run_mode = .false.
    run_mode = 1  ! default run mode
    use_gll = .false.

    i = 1
    do while (i <= argc)
      select case (trim(argv(i)))
      case ('-m', '--model')
        call get_str_value(argv, i, usage, model_name, has_model)
      case ('-r', '--run-mode')
        call get_int_value(argv, i, usage, run_mode, has_run_mode)
      case ('-g', '--use-gll')
        call set_flag(argv, i, usage, use_gll)
      case default
        call arg_error(usage, 'ERROR: Unknown argument: '//trim(argv(i)))
      end select
      i = i + 1
    enddo

    if (.not. has_model) call arg_error(usage, 'ERROR: Model name not set')
    if (run_mode /= 1 .and. run_mode /= 2) call arg_error(usage, 'ERROR: run-mode must be 1 or 2')
  end subroutine parse_args_post_process

  subroutine parse_args_optimize()
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: argv
    character(len=MAX_STRING_LEN) :: usage
    integer :: i, argc
    logical :: has_model

    usage = 'Usage: xfwat_optimize -m <model> [-g] [-h]'

    call read_args(usage, argv, argc)

    has_model = .false.
    use_gll = .false.

    i = 1
    do while (i <= argc)
      select case (trim(argv(i)))
      case ('-m', '--model')
        call get_str_value(argv, i, usage, model_name, has_model)
      case ('-g', '--use-gll')
        call set_flag(argv, i, usage, use_gll)
      case default
        call arg_error(usage, 'ERROR: Unknown argument: '//trim(argv(i)))
      end select
      i = i + 1
    enddo

    if (.not. has_model) call arg_error(usage, 'ERROR: Model name not set')
  end subroutine parse_args_optimize

  subroutine parse_args_mesh_databases()
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: argv
    character(len=MAX_STRING_LEN) :: usage
    integer :: i, argc
    logical :: has_simu_type

    usage = 'Usage: xfwat_mesh_databases -s <simu_type> [-h]'

    call read_args(usage, argv, argc)

    has_simu_type = .false.

    i = 1
    do while (i <= argc)
      select case (trim(argv(i)))
      case ('-s', '--simu-type', '--simu_type')
        call get_str_value(argv, i, usage, simu_type, has_simu_type)
      case default
        call arg_error(usage, 'ERROR: Unknown argument: '//trim(argv(i)))
      end select
      i = i + 1
    enddo

    if (.not. has_simu_type) call arg_error(usage, 'ERROR: simu-type not set')
  end subroutine parse_args_mesh_databases

  ! Parse xspecfwat options, returning the first model index and iteration count.
  ! Also sets simu_type in config; MPI must be initialized for help/error exits.
  subroutine parse_invert_args(first, iteration_count)
    integer, intent(out) :: first, iteration_count
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: argv
    character(len=MAX_STRING_LEN) :: usage, model
    integer :: i, argc
    logical :: has_model, has_simu_type, has_count

    usage = 'Usage: xspecfwat -m <M00-M98> -s noise|tele|leq [-n <iterations>] [-h]'

    call read_args(usage, argv, argc)

    has_model = .false.
    has_simu_type = .false.
    has_count = .false.
    first = -1
    iteration_count = 1

    i = 1
    do while (i <= argc)
      select case (trim(argv(i)))
      case ('-m', '--model')
        call get_str_value(argv, i, usage, model, has_model)
      case ('-s', '--simu-type', '--simu_type')
        call get_str_value(argv, i, usage, simu_type, has_simu_type)
      case ('-n', '--iterations')
        call get_int_value(argv, i, usage, iteration_count, has_count)
      case default
        call arg_error(usage, 'ERROR: Unknown argument: '//trim(argv(i)))
      end select
      i = i + 1
    enddo

    if (.not. has_model) call arg_error(usage, 'ERROR: Model name not set')
    if (.not. has_simu_type) call arg_error(usage, 'ERROR: simu-type not set')
    if (len_trim(model) /= 3 .or. model(1:1) /= 'M' .or. verify(model(2:3), '0123456789') /= 0) &
      call arg_error(usage, 'ERROR: Model must have the form M00 through M98')
    read(model(2:3), '(I2)') first
    if (first > 98) call arg_error(usage, 'ERROR: Model must have the form M00 through M98')
    if (.not. any(INV_TYPE_NAMES == simu_type)) &
      call arg_error(usage, 'ERROR: Unknown simu-type: '//trim(simu_type))
    if (iteration_count < 1 .or. iteration_count > 99 - first) &
      call arg_error(usage, 'ERROR: Requested iterations must produce models no later than M99')
  end subroutine parse_invert_args

  ! read all command line arguments; print usage and exit if -h/--help is given
  subroutine read_args(usage, argv, argc)
    character(len=*), intent(in) :: usage
    character(len=MAX_STRING_LEN), dimension(:), allocatable, intent(out) :: argv
    integer, intent(out) :: argc
    integer :: i, ier

    argc = command_argument_count()
    allocate(argv(argc))
    do i = 1, argc
      call get_command_argument(i, argv(i), status=ier)
      if (ier /= 0) call arg_error(usage, 'ERROR: Cannot read argument (too long?)')
    enddo

    if (any(argv == '-h') .or. any(argv == '--help')) then
      if (worldrank == 0) print '(a)', trim(usage)
      call finalize_mpi()
      stop
    endif
  end subroutine read_args

  ! take argv(i+1) as the value of option argv(i), and advance i to it;
  ! seen flags the option as given, so a repeated option is rejected
  subroutine get_str_value(argv, i, usage, val, seen)
    character(len=*), dimension(:), intent(in) :: argv
    integer, intent(inout) :: i
    character(len=*), intent(in) :: usage
    character(len=*), intent(out) :: val
    logical, intent(inout) :: seen

    if (seen) call arg_error(usage, 'ERROR: Duplicate argument: '//trim(argv(i)))
    if (i >= size(argv)) call arg_error(usage, 'ERROR: Missing value for '//trim(argv(i)))
    if (len_trim(argv(i+1)) == 0 .or. argv(i+1)(1:1) == '-') &
      call arg_error(usage, 'ERROR: Missing value for '//trim(argv(i)))
    i = i + 1
    val = argv(i)
    seen = .true.
  end subroutine get_str_value

  ! same as get_str_value, but the value must be a non-negative integer
  subroutine get_int_value(argv, i, usage, ival, seen)
    character(len=*), dimension(:), intent(in) :: argv
    integer, intent(inout) :: i
    character(len=*), intent(in) :: usage
    integer, intent(out) :: ival
    logical, intent(inout) :: seen
    character(len=MAX_STRING_LEN) :: str
    integer :: ier

    call get_str_value(argv, i, usage, str, seen)
    ier = verify(trim(str), '0123456789')
    if (ier == 0) read(str, *, iostat=ier) ival
    if (ier /= 0) call arg_error(usage, 'ERROR: Invalid integer for '//trim(argv(i-1))//': '//trim(str))
  end subroutine get_int_value

  ! switch option without value
  subroutine set_flag(argv, i, usage, flag)
    character(len=*), dimension(:), intent(in) :: argv
    integer, intent(in) :: i
    character(len=*), intent(in) :: usage
    logical, intent(inout) :: flag

    if (flag) call arg_error(usage, 'ERROR: Duplicate argument: '//trim(argv(i)))
    flag = .true.
  end subroutine set_flag

  subroutine arg_error(usage, msg)
    character(len=*), intent(in) :: usage, msg

    if (worldrank == 0) print '(a)', trim(usage)
    call exit_MPI(worldrank, msg)
  end subroutine arg_error

end module argparse
