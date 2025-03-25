program fwat_post_proc
  use fwat_mpi
  use config
  use post_processing
  use input_params, fpar => fwat_par_global
  use common_lib, only: get_simu_type
  use argparse, only: parse_args_post_process

  implicit none

  type(PostFlow) :: fpp
  integer :: itype

  call init_mpi()
  call init_mpi_fwat()

  call parse_args_post_process()

  call fpar%read(FWAT_PAR_FILE)

  if (count(fpar%postproc%INV_TYPE) > 1) then
    is_joint = .true.
  else
    is_joint = .false.
  endif

  call fpp%init()

  do itype = 1, NUM_INV_TYPE
    if (fpar%postproc%INV_TYPE(itype)) then
      ! set simu type

      ! generate kernels for this type
      call fpp%generate_for_type(itype)
    
      ! sum kernels for this type
      call fpp%sum_kernel()

      if (.not. is_joint .and. fpar%postproc%IS_PRECOND) call fpp%sum_precond()

      ! smooth kernels
      call fpp%smooth_kernel()

      ! write kernels
      call fpp%write()
    end if
  end do

  if (is_joint) then
    call fpp%sum_joint_kernel()
  endif
  call log%write('*******************************************', .false.)
  call log%write('********** POST-PROCESSING DONE ***********', .false.)
  call log%write('*******************************************', .false.)

  call finalize_mpi()
    
end program fwat_post_proc