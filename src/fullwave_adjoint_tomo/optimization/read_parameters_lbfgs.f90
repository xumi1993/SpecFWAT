! This program is used to read parameters for add_model_iso_lbfgs
! Author: Kai Wang, wangkaim8@gmail.com
! University of Toronto, ON, Canada
! Last modified: Dec 24, 2018


subroutine read_parameters_lbfgs()
use tomography_par
  implicit none
  integer :: ier
  character(len=MAX_STRING_LEN) :: s_step_fac,arg
  character(len=MAX_STRING_LEN) :: s_iter_start,s_iter_current

  ! subjective step length to multiply to the gradient
  ! e.g. step_fac = 0.03

  call get_command_argument(1,s_step_fac)
  call get_command_argument(2,s_iter_start)
  call get_command_argument(3,s_iter_current)

  if (command_argument_count() /= 7) then
    call usage()
  endif

  ! read in parameter information
  read(s_step_fac,*) step_fac
  read(s_iter_start,*) iter_start
  read(s_iter_current,*) iter_current

  if (myrank == 0) print *, 'starting iteration for this period band:',iter_start
  if (myrank == 0) print *, 'current iteration:',iter_current ! safety check
  if (abs(step_fac) < 1.e-15) then
    print *,'Error: step factor ',step_fac,' is too small and will lead to no update...'
    call exit_MPI(myrank,'Error step factor too small')
  endif
  if (iter_start < 0 .or. iter_current <0) then
    print *,'Error: iter number is negative ...'
    call exit_MPI(myrank,'Error iter number should be positive')
  endif

  ! input directory for old model files
  call get_command_argument(4,arg)
  if (len_trim(arg) > 0) then
    if (arg(len_trim(arg):len_trim(arg)) /= '/') then
      INPUT_MODEL_DIR = trim(arg) // '/'
    else
      INPUT_MODEL_DIR = trim(arg)
    endif
  endif

  ! output directory for new model files
  call get_command_argument(5,arg)
  if (len_trim(arg) > 0) then
    if (arg(len_trim(arg):len_trim(arg)) /= '/') then
      OUTPUT_MODEL_DIR = trim(arg) // '/'
    else
      OUTPUT_MODEL_DIR = trim(arg)
    endif

    ! sets output directory for statistics to same directory
    if (PRINT_STATISTICS_FILES) then
      OUTPUT_STATISTICS_DIR = trim(OUTPUT_MODEL_DIR)
    endif
  endif

  ! input directory which holds old (summed) kernel files
  call get_command_argument(6,arg)
  if (len_trim(arg) > 0) then
    if (arg(len_trim(arg):len_trim(arg)) /= '/') then
      KERNEL_OLD_DIR = trim(arg) // '/'
    else
      KERNEL_OLD_DIR = trim(arg)
    endif
  endif

  ! input directory which holds new (summed) kernel files
  call get_command_argument(7,arg)
  if (len_trim(arg) > 0) then
    if (arg(len_trim(arg):len_trim(arg)) /= '/') then
      INPUT_KERNELS_DIR = trim(arg) // '/'
    else
      INPUT_KERNELS_DIR = trim(arg)
    endif
  endif

  ! statistics
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_step_fac',status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(OUTPUT_STATISTICS_DIR)//'statistics_step_fac'
      print *,'Please make sure that directory '//trim(OUTPUT_STATISTICS_DIR)//' exists...'
      print *
      stop 'Error opening statistics file'
    endif
    write(IOUT,'(1e24.12)') step_fac
    close(IOUT)
  endif


contains

  subroutine usage()

  implicit none

  if (myrank == 0) then
    print *,'Usage: add_model_lbfgs step_factor iter_star iter_current OLD-MODEL-DIR NEW-MODEL-DIR OLD-KERNELS-DIR NEW-KERNELS-DIR'
    print *
    print *,'with'
    print *,'  step_factor   - factor to scale gradient (e.g. 0.03 for 3 percent update)'
    print *,'  iter_start    - Start number of iteration (e.g. 1)'
    print *,'  iter_current  - Current number of iteration (e.g. 8)'
    print *,'  OLD-MODEL-DIR/  - (optional) directory which will hold old model files (e.g. vp_new.bin,..)'
    print *,'  NEW-MODEL-DIR/  - (optional) directory which will hold new model files (e.g. vp_new.bin,..)'
    print *,'  OLD-KERNELS-DIR/ - (optional) directory which holds old summed kernels (e.g. alpha_kernel.bin,..)'
    print *,'  NEW-KERNELS-DIR/ - (optional) directory which holds new summed kernels (e.g. alpha_kernel.bin,..)'
    print *
    print *
  endif
  call synchronize_all()
  call exit_MPI(myrank,'Error usage: add_model step_factor iter_star iter_current OLD-MODEL-DIR NEW-MODEL-DIR &
      &OLD-KERNELS-DIR NEW-KERNELS-DIR')

  end subroutine usage
end subroutine read_parameters_lbfgs

