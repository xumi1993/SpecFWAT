program fwat_post_proc
  use fwat_mpi
  use config
  use post_processing
  use input_params, fpar => fwat_par_global
  use common_lib, only: get_simu_type

  implicit none

  type(PostFlow) :: fpp
  integer :: itype

  call init_mpi()
  call init_mpi_fwat()

  if (count(fpar%postproc%INV_TYPE) > 1) then
    is_joint = .true.
  else
    is_joint = .false.
  endif

  call fpp%init()

  do itype = 1, NUM_INV_TYPE
    if (fpar%postproc%INV_TYPE(itype)) then
      ! set simu type
      simu_type = fpp%simu_types(itype)
      call fpar%select_simu_type()

      ! generate kernels for this type
      call fpp%generate_for_type()

      ! sum kernels for this type
      call fpp%sum_kernel()

      if (.not. is_joint .or. fpar%postproc%IS_PRECOND) call fpp%sum_precond()

      ! smooth kernels
      call fpp%smooth_kernel()

      ! write kernels
      call fpp%write()
    end if
  end do

  if (is_joint) then
    call fpp%sum_joint_kernel()
  endif

  call finalize_mpi()
    
end program fwat_post_proc