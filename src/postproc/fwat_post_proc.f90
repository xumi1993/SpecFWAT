program fwat_post_proc
  use fwat_mpi
  use config
  use post_processing
  use input_params, fpar => fwat_par_global
  use argparse, only: parse_args_post_process
  use param_check, only: check_model_name, check_postproc_params, check_nevents

  implicit none

  type(PostFlow) :: fpp
  integer :: itype, iter

  call init_mpi()
  call init_mpi_fwat()

  call parse_args_post_process()

  call fpar%read(FWAT_PAR_FILE)
  call read_parameter_file(.true.)

  call check_model_name(iter)

  call fpp%init()
  do itype = 1, NUM_INV_TYPE
    if (fpar%postproc%INV_TYPE(itype)) then

      ! generate kernels for this type
      call fpp%init_for_type(itype)

      ! init_for_type selected the data type and read its source list, so the
      ! smoothing and tapering parameters of this type can be checked here.
      call check_postproc_params(iter)
      call check_nevents()

      ! sum kernels for this type
      if (run_mode == 1) then
        call fpp%sum_kernel()
      elseif (run_mode == 2) then
        call fpp%read_sum_kernel()
      endif
      
      if (fpar%postproc%IS_HESS_PRECOND) then
        call fpp%apply_precond()
      else
        if ((.not. (is_joint .and. (itype == 1 .or. itype == 3))) .and. &
            run_mode == 1) call fpp%sum_precond()
      endif

      call fpp%pde_smooth()
      
      if (use_gll) then
        call fpp%taper_kernel_gll()
        call fpp%write_gradient_gll()
      else
        call fpp%taper_kernel_grid()
        call fpp%write_gradient_grid()
      endif

      ! remove event kernels
      call remove_ekernel()

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
  call log%finalize()
  
  call synchronize_all()

  call finalize_mpi()
    
end program fwat_post_proc