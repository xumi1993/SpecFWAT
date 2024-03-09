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

 subroutine run_sele_ctrlgrp(model,evtsetb,evtsete,nctrl)

  use fullwave_adjoint_tomo_par
  use fwat_input
  use tomography_par, only: USE_ALPHA_BETA_RHO,USE_ISO_KERNELS
  !!! for model update
  use tomography_model_iso
  use tomography_kernels_iso
  implicit none


  character(len=MAX_STRING_LEN)            :: model 
  character(len=MAX_STRING_LEN)            :: evtsetb,evtsete 
  character(len=MAX_STRING_LEN)            :: strset
  character(len=MAX_STRING_LEN)            :: ekernel_dir_list, output_dir
  integer                                  :: iset, setb, sete,ier,nctrl

  character(len=MAX_STRING_LEN)            :: evtset_file 
  character(len=MAX_STRING_LEN)            :: evtnm 
  character(len=MAX_STRING_LEN)            :: fwd_dir 


    
  read(evtsetb(4:),'(I3)') setb 
  read(evtsete(4:),'(I3)') sete
  ekernel_dir_list='optimize/ekernel_dir.lst'
  output_dir='optimize/SUM_KERNELS_'//trim(model)
  !**************** Build directories for the storage ****************************
  call read_fwat_par_file()
  if (myrank == 0) then
     open(unit=OUT_FWAT_LOG,file='output_dynmbat_log_'//trim(model)//'.'//trim(evtsetb)&
                                 //'-'//trim(evtsete)//'.txt')
     write(OUT_FWAT_LOG,*) 'Running XDYNMBAT SELE CTRL GRP !!!'
     write(OUT_FWAT_LOG,*) 'model,evtsetb,evtsete: ',trim(model),' ', trim(evtsetb),'-',trim(evtsete)
     write(OUT_FWAT_LOG,*) 'SAVE_OUTPUT_EACH_EVENT: ',tomo_par%SAVE_OUTPUT_EACH_EVENT
     call system('mkdir -p optimize')
    
     open(unit=400,file=trim(ekernel_dir_list),iostat=ier)
     do iset=setb,sete
        write(strset,'(I0)') iset
        evtset_file='src_rec/sources_'//'set'//trim(strset)//'.dat'
        open(unit=401,file=trim(evtset_file),status='old',iostat=ier)
        if (ier /=0) then
          print *,'Error could not open source subset file: ',trim(evtset_file) 
          call exit_mpi(myrank,'Error opening source subset file')
        endif
        do 
          read(401,*,iostat=ier) evtnm 
          if (ier /=0) exit
          fwd_dir='solver/'//trim(model)//'.set'//trim(strset)//'/'//trim(evtnm)
          write(400,*) trim(fwd_dir)//'/EKERNEL'
        enddo
        close(401)
     enddo
     close(400)
     write(OUT_FWAT_LOG,*) '*******************************************************'
    
  endif 
  call synchronize_all()

! -----------------------------------------------------------------
! sum 
  if(myrank==0) write(OUT_FWAT_LOG,*) 'This is compute_ctrl_grp ...'
  call compute_ctrl_grp(ekernel_dir_list,output_dir,model,nctrl)
! -----------------------------------------------------------------
  
  if(myrank==0)  write(OUT_FWAT_LOG,*) '*******************************************************'
  if(myrank==0)  write(OUT_FWAT_LOG,*) 'Finished DYNMBAT select control group here!!!'
  if(myrank==0) close(OUT_FWAT_LOG)


end subroutine run_sele_ctrlgrp
