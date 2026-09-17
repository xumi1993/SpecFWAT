program fwat_post_proc
  use fwat_mpi
  use config
  use post_processing, only: remove_ekernel
  uer post_rtm, only: PostRTM
  use input_params, fpar => fwat_par_global
  use argparse, only: parse_args_post_process
  use param_check, only: check_model_name, check_postproc_params, check_nevents

  implicit none

  type(PostRTM) :: fpp

  call init_mpi()
  call init_mpi_fwat()

  model_name = "RTM"

  call fpar%read(FWAT_PAR_FILE)
  call read_parameter_file(.true.)

  call fpp%init()

  call fpp%sum_kernel()
      
  call fpp%apply_precond()
    
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