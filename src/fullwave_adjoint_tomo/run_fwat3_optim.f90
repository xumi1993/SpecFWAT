subroutine run_optim()

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
  use postproc_sub

  implicit none

  double precision                                                :: misfit, chi, chi0, this_misfit, misfit0(NUM_INV_TYPE)
  integer                                                         :: nchan, i, j, sit
  integer                                                         :: imod_current,imod_up, imod_down
  integer,                               parameter                :: LOG_UNIT = 898
  character(len=MAX_STRING_LEN)                                   :: model_ls, evtset, msg, strstep

! ------------------------------------------------------------------
! optimize (SD, CG or LBFGS)
  call world_rank(myrank)
  read(model(2:),'(I2.2)') imod_current
  imod_up=imod_current+1
  imod_down=imod_current-1
  if(myrank==0) then
    open(unit=LOG_UNIT,file='output_fwat3_log_'//trim(model)//'.txt')
    write(LOG_UNIT,*) '*******************************************************'
    call write_timestamp_log(LOG_UNIT, 'This is optimize ...')
    flush(LOG_UNIT)
  endif
  INPUT_MODEL_DIR='optimize/MODEL_'//trim(model)//'/'
  if (imod_down<tomo_par%ITER_START) then
    KERNEL_OLD_DIR='none'
  else
    KERNEL_OLD_DIR='optimize/SUM_KERNELS_'//trim(model_prev)//'/'
  endif
  INPUT_KERNELS_DIR='optimize/SUM_KERNELS_'//trim(model)//'/'
  PRINT_STATISTICS_FILES=.false.
  INPUT_DATABASES_DIR=trim(LOCAL_PATH)//'/'
  iter_start=tomo_par%ITER_START
  iter_current=imod_current
  if (tomo_par%DO_LS) then
    ! read current misfit
    step_fac = tomo_par%MAX_SLEN
    misfit = 0.0
    misfit0 = 0.0
    do i = 1, NUM_INV_TYPE
      if (tomo_par%INV_TYPE(i)) then
        call read_misfit(tomo_par%INV_TYPE_NAME(i), model, chi, nchan)
        call read_misfit(tomo_par%INV_TYPE_NAME(i), 'M00', chi0, nchan)
        misfit0(i) = chi0
        misfit = misfit + chi / misfit0(i)
      endif
    enddo
    call synchronize_all()
      
    do sit = 1, tomo_par%MAX_SUB_ITERS
      write(strstep,'(F5.3)') step_fac
      if (myrank == 0) then
        write(msg, '("Starting ",I0,"th sub-iteration with step_length = ",F5.3)') sit, step_fac
        call write_timestamp_log(LOG_UNIT, trim(msg))
        flush(LOG_UNIT)
      endif

      ! get L-BFGS direction
      OUTPUT_MODEL_DIR=trim(LOCAL_PATH)//'/' 
      if (myrank == 0) then
        call write_timestamp_log(LOG_UNIT, 'Write tmp model to: '//trim(OUTPUT_MODEL_DIR))
        flush(LOG_UNIT)
      endif
      call model_update_opt()
      call synchronize_all()
      call generate_database_fwat()
      call synchronize_all()

      ! sub-iterations
      this_misfit = 0.0
      do i = 1, NUM_INV_TYPE
        model_ls =  trim(model)//'_step'//trim(strstep)
        if (tomo_par%INV_TYPE(i)) then
          type_name = tomo_par%INV_TYPE_NAME(i)  
          call select_set_range()
          do j = isetb, isete
            write(evtset,'("set",I0)') j
            call run_linesearch(model_ls, evtset, type_name)
          enddo
          call read_misfit(type_name, model_ls, chi, nchan)
          this_misfit = this_misfit + chi / misfit0(i)
        endif
      enddo
      call synchronize_all()

      ! check misfit
      if (this_misfit < misfit) then
        write(msg, '("Norm misfit reduced from ",F20.6," to ",F20.6)') misfit, this_misfit
        call write_timestamp_log(LOG_UNIT, msg)
        call write_timestamp_log(LOG_UNIT, 'Accept model')
        flush(LOG_UNIT)
        exit
      else
        write(msg, '("Norm misfit increased from ",F20.6," to ",F20.6)') misfit, this_misfit
        call write_timestamp_log(LOG_UNIT, msg)
        flush(LOG_UNIT)
      endif
    enddo
  else
    step_fac=tomo_par%MAX_SLEN
    if (trim(LOCAL_PATH) /= 'optimize/MODEL_'//trim(model) .and. &
        trim(LOCAL_PATH) /= 'optimize/MODEL_'//trim(model)//'/' .and. &
        trim(LOCAL_PATH) /= './optimize/MODEL_'//trim(model) .and. &
        trim(LOCAL_PATH) /= './optimize/MODEL_'//trim(model)//'/') then
      OUTPUT_MODEL_DIR=trim(LOCAL_PATH)//'/' ! run second time to be called for next iteration
      call model_update_opt()
    endif
  endif
  OUTPUT_MODEL_DIR='optimize/MODEL_'//trim(model_next)//'/'
  if (myrank == 0) call write_timestamp_log(LOG_UNIT, 'Write model to: '//trim(OUTPUT_MODEL_DIR))
  call model_update_opt() ! run first time to save model 
! -----------------------------------------------------------------
  if (myrank==0) write(LOG_UNIT,*) '*******************************************************'
  if (myrank==0) write(LOG_UNIT,*) 'Finish stage3 optimization'
  if (myrank==0) close(LOG_UNIT)
  ! if(myrank==0) close(IMAIN)


end subroutine run_optim


subroutine run_optim_man()
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
  use postproc_sub

  implicit none

  integer                                         :: i,istep,imod_current,imod_up, imod_down
  character(len=MAX_STRING_LEN)                   :: strstep

  read(model(2:),'(I2.2)') imod_current
  imod_up=imod_current+1
  imod_down=imod_current-1
  ! optimize (SD, CG or LBFGS)
  ! if(myrank==0) then
    ! call write_timestamp_log(OUT_FWAT_LOG, 'This is optimize ...')
    ! flush(OUT_FWAT_LOG)
  ! endif
  INPUT_MODEL_DIR='optimize/MODEL_'//trim(model)//'/'
  if (imod_down<tomo_par%ITER_START) then
    KERNEL_OLD_DIR='none'
  else
    KERNEL_OLD_DIR='optimize/SUM_KERNELS_'//trim(model_prev)//'/'
  endif
  INPUT_KERNELS_DIR='optimize/SUM_KERNELS_'//trim(model)//'/'
  PRINT_STATISTICS_FILES=.false.
  INPUT_DATABASES_DIR=trim(LOCAL_PATH)//'/'
  iter_start=tomo_par%ITER_START
  iter_current=imod_current
  if (tomo_par%DO_LS) then
    do istep=1,tomo_par%NUM_STEP
      step_fac=tomo_par%STEP_LENS(istep)
      write(strstep,'(F5.3)') step_fac
      OUTPUT_MODEL_DIR='optimize/MODEL_'//trim(model)//'_step'//trim(strstep)//'/'
      call system('mkdir -p '//trim(OUTPUT_MODEL_DIR))
      ! if (myrank == 0) call write_timestamp_log(OUT_FWAT_LOG, 'Write model to: '//trim(OUTPUT_MODEL_DIR))
      call model_update_opt()
    enddo
  else
    step_fac=tomo_par%MAX_SLEN
    OUTPUT_MODEL_DIR='optimize/MODEL_'//trim(model_next)//'/'
    ! if (myrank == 0) call write_timestamp_log(OUT_FWAT_LOG, 'Write model to: '//trim(OUTPUT_MODEL_DIR))
    call model_update_opt() ! run first time to save model 
    if (trim(LOCAL_PATH) /= 'optimize/MODEL_'//trim(model) .and. &
        trim(LOCAL_PATH) /= 'optimize/MODEL_'//trim(model)//'/' .and. &
        trim(LOCAL_PATH) /= './optimize/MODEL_'//trim(model) .and. &
        trim(LOCAL_PATH) /= './optimize/MODEL_'//trim(model)//'/') then
      OUTPUT_MODEL_DIR=trim(LOCAL_PATH)//'/' ! run second time to be called for next iteration
      call model_update_opt()
    endif
  endif 
end subroutine run_optim_man