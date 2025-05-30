program fwat_post_proc
  use fwat_mpi
  use config
  use post_processing
  use input_params, fpar => fwat_par_global
  use argparse, only: parse_args_post_process

  implicit none

  type(PostFlow) :: fpp
  integer :: itype

  call init_mpi()
  call init_mpi_fwat()

  call parse_args_post_process()

  call fpar%read(FWAT_PAR_FILE)
  call read_parameter_file(.true.)

  call fpp%init(.true.)
  do itype = 1, NUM_INV_TYPE
    if (fpar%postproc%INV_TYPE(itype)) then

      ! generate kernels for this type
      call fpp%init_for_type(itype)
    
      ! sum kernels for this type
      if (run_mode == 1) then
        call fpp%sum_kernel()
      elseif (run_mode == 2) then
        call fpp%read_sum_kernel()
      endif
      
      if (fpar%postproc%IS_PRECOND) then
        call fpp%apply_precond()
      else
        if (.not. (is_joint .and. itype == 1)) call fpp%sum_precond()
      endif

      call fpp%pde_smooth()

      call fpp%taper_kernel_grid()

      call fpp%write_gradient_grid()

      ! remove event kernels
      call fpp%remove_ekernel()

      call fpar%acqui%finalize()

      call fpp%finalize()

      call synchronize_all()
    end if
  end do

  if (is_joint) then
    call sum_joint_kernel_grid()
  endif

  call log%write('*******************************************', .false.)
  call log%write('********** POST-PROCESSING DONE ***********', .false.)
  call log%write('*******************************************', .false.)

  call finalize_mpi()
    
end program fwat_post_proc