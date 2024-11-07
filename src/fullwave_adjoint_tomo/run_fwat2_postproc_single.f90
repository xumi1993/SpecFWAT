!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
!

subroutine run_fwat2_postproc_single(model, simu_type)

  use fullwave_adjoint_tomo_par
  use fwat_input
  use fwat_utils
  use tomography_par, only: USE_ALPHA_BETA_RHO,USE_ISO_KERNELS, step_fac, iter_start,iter_current, &
           INPUT_MODEL_DIR,OUTPUT_MODEL_DIR,KERNEL_OLD_DIR,INPUT_KERNELS_DIR,PRINT_STATISTICS_FILES
  !!! for reading Database
  use specfem_par
  ! use specfem_par_elastic, only: ELASTIC_SIMULATION,ispec_is_elastic,rho_vp,rho_vs,min_resolved_period
  ! use specfem_par_acoustic, only: ACOUSTIC_SIMULATION,ispec_is_acoustic
  ! use specfem_par_poroelastic, only: POROELASTIC_SIMULATION,ispec_is_poroelastic,rho_vpI,rho_vpII,rho_vsI, &
    ! phistore,tortstore,rhoarraystore
  !!! for model update
  use tomography_model_iso
  use tomography_kernels_iso
  use taper3d
  use constants, only: IMAIN
  use postproc_sub, only: post_proc, sum_joint_kernels, read_database, get_kernel_names, type_name, this_model=>model, &
                          model_prev, model_next, imod_current,imod_up, imod_down, get_model_idx
  implicit none

  ! real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  ! real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  ! real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  ! real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  ! real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob


  character(len=MAX_STRING_LEN)                   :: model, simu_type
  ! character(len=MAX_STRING_LEN)                   :: evtsetb,evtsete,is_smooth 
  character(len=MAX_STRING_LEN)                   :: strstep,strinv
  ! character(len=MAX_STRING_LEN) :: input_dir,output_dir,ekernel_dir_list
  integer                                         :: i,istep
  !real :: t1, t2
  ! logical :: USE_GPU ! TODO: Smoothing using GPU
  
  call world_rank(myrank)
  this_model=model
  call get_model_idx()
  
  ! read(evtsetb(4:),'(I3)') setb 
  ! read(evtsete(4:),'(I3)') sete
  ! ekernel_dir_list='optimize/ekernel_dir.lst'
  ! output_dir='optimize/SUM_KERNELS_'//trim(model)
  !**************** Build directories for the storage ****************************

  if (myrank == 0) then
    open(unit=OUT_FWAT_LOG,file='output_fwat2_log_'//trim(model)//'.'//trim(simu_type)//'.txt')
    write(OUT_FWAT_LOG,*) 'Running XFWAT2_POSTPROC_OPT !!!'
  !  write(OUT_FWAT_LOG,*) 'model,evtsetb,evtsete: ',trim(model),' ', trim(evtsetb),'-',trim(evtsete)
    write(OUT_FWAT_LOG,*) 'INV_TYPE: ',trim(simu_type)
    do i = 1, 3
      if (tomo_par%INV_TYPE(i)) then
        write(OUT_FWAT_LOG,*) 'DATA TYPE: ',trim(tomo_par%INV_TYPE_NAME(i))
        write(OUT_FWAT_LOG,*) '-- SIGMA_H, SIGMA_V: ', tomo_par%NOISE_SIGMA_H,tomo_par%NOISE_SIGMA_V
        write(OUT_FWAT_LOG,*) '-- WEIGHT: ', tomo_par%JOINT_WEIGHT(i)
      endif
    enddo
    if (strinv == 'joint') write(OUT_FWAT_LOG,*) 'NORM_TYPE: '//trim(tomo_par%NORM_TYPE)
    write(OUT_FWAT_LOG,*) 'USE_SPH_SMOOTH: ',tomo_par%USE_SPH_SMOOTH
    write(OUT_FWAT_LOG,*) 'model_prev,model_next: ',trim(model_prev),' ',trim(model_next)
    write(OUT_FWAT_LOG,*) 'OPT_METHOD: ',trim(tomo_par%OPT_METHOD) 
    write(OUT_FWAT_LOG,*) 'DO_LS: ',tomo_par%DO_LS 
    write(OUT_FWAT_LOG,*) 'MAX_SLEN: ',tomo_par%MAX_SLEN
    write(OUT_FWAT_LOG,*) 'STEP_LENS: ',tomo_par%STEP_LENS 
    write(OUT_FWAT_LOG,*) 'VPVS_RATIO_RANGE: ',tomo_par%VPVS_RATIO_RANGE
    write(OUT_FWAT_LOG,*) '*******************************************************'
    flush(OUT_FWAT_LOG)
    call system('mkdir -p optimize')
    if(imod_current==0) then
      call system('mkdir -p optimize/MODEL_M00 && cp model_initial/* optimize/MODEL_M00/')
    endif
    call system('mkdir -p optimize/MODEL_'//trim(model_next))
  endif 
  call synchronize_all()
! -----------------------------------------------------------------
! Sum and smooth kernels
  if (is_read_database) call read_database()
  call get_kernel_names()
  call synchronize_all()
  type_name = simu_type
  call post_proc()
  ! do i = 1,NUM_INV_TYPE
  !   type_name = tomo_par%INV_TYPE_NAME(i)
  !   if (type_name == simu_type .and. tomo_par%INV_TYPE(i)) then
  !     call post_proc()
  !   endif
  ! enddo
  ! if (count(tomo_par%INV_TYPE)>1) then
  !   call sum_joint_kernels()
  !   if ((.not. VERBOSE_MODE) .and. imod_current /= tomo_par%ITER_START .and. myrank == 0) then
  !     do i = 1,NUM_INV_TYPE
  !       if (type_name == simu_type .and. tomo_par%INV_TYPE(i)) then
  !         call system('rm -rf ./optimize/SUM_KERNELS_'//trim(model)//'_'//trim(tomo_par%INV_TYPE_NAME(i)))
  !       endif
  !     enddo
  !   endif
  ! endif
  call synchronize_all()
  
  if(myrank==0)  write(OUT_FWAT_LOG,*) '*******************************************************'
  if(myrank==0)  write(OUT_FWAT_LOG,*) 'Finished FWAT stage2 here!!!'
  if(myrank==0) close(OUT_FWAT_LOG)
  if(myrank==0) close(IMAIN)


end subroutine run_fwat2_postproc_single


