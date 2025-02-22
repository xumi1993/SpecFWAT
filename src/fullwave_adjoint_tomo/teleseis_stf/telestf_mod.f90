module telestf_mod
   use fullwave_adjoint_tomo_par, only: NRCOMP, RCOMPS, VERBOSE_MODE, Niter
   use fwat_utils, only: bandpass_iir
   use interpolation_mod, only: interpolate_syn_alloc, myconvolution, time_deconv
   use ma_constants, only: NDIM
   use spanlib
   use specfem_par, only: OUTPUT_FILES, MAX_STRING_LEN
   use sacio
   implicit none

contains

subroutine time_iterdeconv(stnm,seismo_dat,seismo_syn,dat_tw,syn_tw,ff,&
                          ttp,tb,te,irec,nrec,npts,beg,del,fstart0,fend0)
   implicit none
   type(sachead) :: head
   integer :: nrec,irec
   integer :: npts, nerr
   real(kind=4) , dimension(npts)  :: seismo_dat,seismo_syn
   real(kind=4) , dimension(nrec,npts,NRCOMP)  :: dat_tw,syn_tw ! windowed dat and syn
   real(kind=4) , dimension(nrec,npts)  :: ff
  !!! for reading sac data
   real(kind=4) :: beg, del
   character(len=MAX_STRING_LEN), dimension(nrec)  :: stnm
   real(kind=4)  :: ttp,tb,te
   !!! for time_iterdecon
   double precision                                  :: t0_inp,t1_inp,dt_inp
   double precision, dimension(NDIM)                 :: dat_inp,syn_inp,rfn_bp
   real(kind=4) , dimension(:), allocatable :: num,den,rfn,tmpl
   integer :: NDIM_CUT
   !!! for PCA
   double precision, intent(in) :: fstart0,fend0

  !===========================================================
  ! preprocessing and cut
  dt_inp=dble(del) !!!! Here, I set dt=0.01 to interpolate data and syn
!   t0_inp=ttp+beg-tb
!   t1_inp=ttp+beg+te
  t0_inp=ttp-tb
  t1_inp=ttp+te
  NDIM_CUT=(t1_inp-t0_inp)/dt_inp+1
  !write(*,*)'irec,dt_inp,t0_inp,t1_inp,NDIM_CUT=',irec,dt_inp,t0_inp,t1_inp,NDIM_CUT
  dat_inp(:)=0.
  syn_inp(:)=0.
  dat_inp(1:npts)=dble(seismo_dat(1:npts))
  syn_inp(1:npts)=dble(seismo_syn(1:npts))
  ! cut to P window -5 to 45 ps
  call interpolate_syn_alloc(dat_inp,dble(beg),dble(del),npts,t0_inp,dt_inp,NDIM_CUT)
  call interpolate_syn_alloc(syn_inp,dble(beg),dble(del),npts,t0_inp,dt_inp,NDIM_CUT)
  call interpolate_syn_alloc(dat_inp,t0_inp,dt_inp,NDIM_CUT,dble(beg),dble(del),npts)
  call interpolate_syn_alloc(syn_inp,t0_inp,dt_inp,NDIM_CUT,dble(beg),dble(del),npts)
  dat_tw(irec,1:npts,1)=real(dat_inp(1:npts))
  syn_tw(irec,1:npts,1)=real(syn_inp(1:npts))
  NDIM_CUT=npts
  allocate(num(NDIM_CUT))
  allocate(den(NDIM_CUT))
  num=real(dat_inp(1:NDIM_CUT))
  den=real(syn_inp(1:NDIM_CUT))
  if (VERBOSE_MODE) then
     call sacio_newhead(head, del, npts, beg)
     call sacio_writesac(trim(stnm(irec))//'.num', head, dble(num), nerr)
     call sacio_writesac(trim(stnm(irec))//'.den', head, dble(den), nerr)
   !   call wsac1(trim(stnm(irec))//'.num',num,npts,beg,del,nerr)  
   !   call wsac1(trim(stnm(irec))//'.den',den,npts,beg,del,nerr) 
  endif
  ! call time_deconv
  allocate(rfn(NDIM_CUT))
  if (maxval(abs(den))==0) then
     rfn=0.
  else
     call time_deconv(num,den,del,NDIM_CUT,Niter,rfn)
  endif
  allocate(tmpl(NDIM_CUT))
  if (maxval(abs(den)).ne.0) then
     call myconvolution(den,rfn,NDIM_CUT,NDIM_CUT,tmpl,0)
     tmpl=tmpl * del
  else
     tmpl=0
  endif
  if (VERBOSE_MODE) then
     call sacio_newhead(head, del, NDIM_CUT, beg)
     call sacio_writesac(trim(stnm(irec))//'.pre', head, dble(tmpl), nerr)
   !   call wsac1(trim(stnm(irec))//'.pre',tmpl,NDIM_CUT,beg,del,nerr)  
  endif

  ! filter deconv traces
  rfn_bp(:)=0.
  rfn_bp(1:NDIM_CUT)=rfn(1:NDIM_CUT)
  call bandpass_iir(rfn_bp,NDIM_CUT,dt_inp,fstart0,fend0)
  if (VERBOSE_MODE) then
     call sacio_writesac(trim(stnm(irec))//'.dec', head, rfn_bp(1:NDIM_CUT), nerr)
   !  call wsac1(trim(stnm(irec))//'.dec',real(rfn_bp),NDIM_CUT,beg,del,nerr)  
  endif
  ff(irec,1:NDIM_CUT)=rfn_bp(1:NDIM_CUT)
  deallocate(num,den,rfn,tmpl)
end subroutine time_iterdeconv

subroutine seis_pca(stnm,dat_tw,syn_tw,ff,nrec,npts,beg,del,stf_array)
   implicit none
   type(sachead) :: head
   integer :: npts, nerr
   integer :: i,j,nrec,irec,icomp
   real(kind=4) , dimension(nrec,npts,NRCOMP)  :: dat_tw,syn_tw
   real(kind=4) , dimension(nrec,npts)  :: ff
   real(kind=4) , dimension(npts,NRCOMP)  :: stf_array
  !!! for reading sac data
   real(kind=4) :: beg, del
   character(len=256) kname
   character(len=256), dimension(nrec)  :: stnm
   !!! for PCA
   real(kind=4) , dimension(:,:), allocatable :: xeof    ! m x m
   real(kind=4) , dimension(:,:), allocatable :: pc      ! n x m
   real(kind=4) , dimension(:,:), allocatable :: pcT      ! m x n
   real(kind=4) , dimension(:,:), allocatable :: pdm      ! n x m
   real(kind=4) , dimension(:), allocatable :: ev      ! m 
   integer :: nkeep 
   real :: sumev
   !!! seis_ampl_scale
   real(kind=4) , dimension(:,:,:), allocatable :: recp_syn
   real(kind=4) , dimension(:), allocatable :: tmpl
   !!! for geting average amplitude factor
   real(kind=4) sumamp,avgamp0,avgamp(NRCOMP)
   real(kind=4) , dimension(:), allocatable :: avgarr

  !===========================================================
  ! apply PCA 
   write(*,*)'run PCA ...'
   nkeep=nrec
   allocate(xeof(nkeep,nkeep))
   allocate(pc(npts,nkeep))
   allocate(pcT(nkeep,npts))
   allocate(pdm(npts,nkeep))
   allocate(ev(nkeep))
   call sl_pca(ff, nkeep, xeof, pc, ev) 
   !!! OUTPUT to files
   open(9,file=trim(OUTPUT_FILES)//'/eigenvalue.dat')
   sumev=0.
   do i=1,nkeep
     write(9,*)ev(i)
     sumev=sumev+ev(i)
   enddo
   close(9)
   open(10,file=trim(OUTPUT_FILES)//'/contribution.dat')
   do i=1,nkeep
      write(10,*) 100.*ev(i)/sumev
   enddo
   close(10)
   ! Calculate the primary data mode reconstucted from the first PC
   pcT=transpose(pc)
   pdm=0. 
   !do i=1,1
    pdm=matmul(pc(:,1:1),transpose(xeof(:,1:1)))
   !enddo
   if (VERBOSE_MODE) then
      do i=1,nkeep
         write(kname,'(A3,I3.3,A4)')'PCs',i,'.sac'
         call sacio_newhead(head, del, npts, beg)
         call sacio_writesac(trim(OUTPUT_FILES)//'/'//trim(kname), head, dble(pc(1:npts,i)), nerr)
         ! call wsac1(trim(OUTPUT_FILES)//'/'//trim(kname),pc(1:npts,i),npts,beg,del,nerr)
         if(nerr .NE. 0) then
            write(*,*)'Error writing SAC File: ', kname
            call exit(-1)
         endif
         write(kname,'(A4,I3.3,A4)')'PDMs',i,'.sac'
         call sacio_writesac(trim(OUTPUT_FILES)//'/'//trim(kname), head, dble(pdm(1:npts,i)), nerr)
         ! call wsac1(trim(OUTPUT_FILES)//'/'//trim(kname),pdm(1:npts,i),npts,beg,del,nerr)

         write(kname,'(A4,I3.3,A4)')'EOFs',i,'.sac'
         call sacio_newhead(head, del, nkeep, beg)
         call sacio_writesac(trim(OUTPUT_FILES)//'/'//trim(kname), head, dble(xeof(1:nkeep,i)), nerr)
         ! call wsac1(trim(OUTPUT_FILES)//'/'//trim(kname),xeof(1:nkeep,i),nkeep,beg,del,nerr)
      enddo  
   endif
   !===========================================================
   ! reconstruct syntheic using primary average STF   
   write(*,*)'reconstruct syntheic using primary average STF ...'
   allocate(recp_syn(nrec,npts,NRCOMP))
   allocate(tmpl(npts))
   call sacio_newhead(head, del, npts, beg)
   do irec=1,nrec
      do icomp=1,NRCOMP
         call myconvolution(syn_tw(irec,1:npts,icomp),pc(1:npts,1),npts,npts,tmpl,0) 
         tmpl=tmpl * del
         recp_syn(irec,1:npts,icomp)=tmpl(1:npts)
         if (VERBOSE_MODE) then
            call sacio_writesac(trim(stnm(irec))//'.recp'//trim(RCOMPS(icomp)), head, dble(tmpl), nerr)
            ! call wsac1(trim(stnm(irec))//'.recp'//trim(RCOMPS(icomp)),tmpl,npts,beg,del,nerr)  
         endif
      enddo
   enddo
   ! call  seis_amp_scale
   write(*,*)'calculate ampl scale factors ...'
   do icomp=1,NRCOMP
      allocate(avgarr(nrec))
      write(*,*)'Inital amp factors:',trim(RCOMPS(icomp))
      open(10,file=trim(OUTPUT_FILES)//'/ampl_fac_'//trim(RCOMPS(icomp))//'.dat')
      sumamp=0.
      j=0
      do i=1,nrec
         if (maxval(abs(recp_syn(i,:,icomp))).ne.0) then
            avgarr(i)=dot_product(dat_tw(i,:,icomp),recp_syn(i,:,icomp))/dot_product(recp_syn(i,:,icomp),recp_syn(i,:,icomp))    
            write(*,*) trim(stnm(i)),avgarr(i)
            sumamp=sumamp+avgarr(i)
            j=j+1
         else
            avgarr(i)=-1000.
         endif
         write(10,*) trim(stnm(i)),avgarr(i)
      enddo
      close(10)
      write(*,*) "There are ", j,' non zero rec found'
      if (j==0) then 
         write(*,*) 'skip icomp: ',trim(RCOMPS(icomp))
         deallocate(avgarr)
         cycle
      endif
      avgamp0=sumamp/j
      write(*,*)'Initial averge amp: ',avgamp0 
      write(*,*)'selected amp factors (|A-A0|<=0.2):'
      sumamp=0.
      j=0
      do i=1,nrec
      !if (abs(avgarr(i)-avgamp0)<=0.2*abs(avgamp0) .and. avgarr(i).ne.-1000.) then
         if (abs(avgarr(i)-avgamp0)<=0.2 .and. avgarr(i).ne.-1000.) then
            write(*,*) avgarr(i)
            sumamp=sumamp+avgarr(i)
            j=j+1
         endif
      enddo
      if(j==0.or.j<=0.1*real(nrec)) then
         write(*,*) 'Error: too little data satisfy |A-A0|<=0.2'
         call exit(-1)
      endif
      avgamp(icomp)=sumamp/j
      write(*,*)'Final averge amp: ',avgamp(icomp)
      stf_array(:,icomp)=pc(1:npts,1)*avgamp(icomp)
      deallocate(avgarr)
      ! get new scaled STFs for Z and R
      call sacio_writesac(trim(OUTPUT_FILES)//'/stf_pca001_'//trim(RCOMPS(icomp))//'.sac', head, dble(stf_array(:,icomp)), nerr)
      ! call wsac1(trim(OUTPUT_FILES)//'/stf_pca001_'//trim(RCOMPS(icomp))//'.sac',stf_array(:,icomp),npts,beg,del,nerr)  
   enddo ! end loop icomp
      
   ! reconstruct new synthetics using new scaled STFs
   write(*,*)'reconstruct syntheic using scaled average STF ...'
   do irec=1,nrec
      do icomp=1,NRCOMP
         call myconvolution(syn_tw(irec,1:npts,icomp),pc(1:npts,1)*avgamp(icomp),npts,npts,tmpl,0)
         tmpl=tmpl * del
         recp_syn(irec,1:npts,icomp)=tmpl(1:npts)
         if (VERBOSE_MODE) then
            call sacio_writesac(trim(stnm(irec))//'.rec'//trim(RCOMPS(icomp)), head, dble(tmpl), nerr)
            ! call wsac1(trim(stnm(irec))//'.rec'//trim(RCOMPS(icomp)),tmpl,npts,beg,del,nerr)
         endif
      enddo
   enddo 
   ! calculate adjoint sources and residuals 

   deallocate(xeof)
   deallocate(pc)
   deallocate(pcT)
   deallocate(pdm)
   deallocate(ev) 


end subroutine seis_pca

end module telestf_mod
