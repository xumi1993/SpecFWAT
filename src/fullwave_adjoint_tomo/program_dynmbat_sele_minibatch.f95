! Kai Wang, Macquarie University
! Aug 06, 2021
!========================================================
! This program use the Mitchell’s best-candidate algorithm
! to find a dynamic mini-batch The first source in the 
! initial mini-batch is chosen randomly, the source furthest
! away from the already selected sources is added to the mini-batch.

! ==== Dynamic mini-batch ====
! Reference: 
! van Herwaarden, D.-P., Boehm, C., Afanasiev, M., Thrastarson, S., Krischer, L.,
! Trampert, J. & Fichtner, A., 2020. Accelerating full-waveform inver- sion using
! dynamic mini-batches, Geophys. J. Int., 221, 1427–1438.

! ==== Mitchell's best-candidate algorithm ====
! Reference: Mitchell, 1991. Spectrally Optimal Sampling for Distribution Ray
!            Tracing.
! Method   : Begin by choosing the first sample at random in a region. To add
!  the (n + 1)th sample, generate m*n uniformly distributed random candidate points
!  (where m is a constant parameter).
!  For each of these random points, compute the distance to the closest of the n
!  points already in the pattern. Then chose the candidate point with the largest
!  closest-point distance, and add it to the pattern. By scaling up the number of random
!  candidate points, in proportion to n, we maintain a constant ratio m of candidates
!  to pattern points in the process. 
!========================================================
program dynmbat_sele_minibatch
   implicit none
   character(len=50)                             :: model,model_prev,model_next 
   character(len=50)                             :: infname,outfname,nk_str
   integer                                       :: i,j,k,ier,nsrc,nk,nsrcb,nsrcr
   integer                                       :: narg
   integer                                       :: imod_current,imod_up, imod_down
   ! source list 
   character(len=50)                             :: kname
   real(kind=4)                                  :: lat,lon,ele,bur  
   character(len=50), dimension(:), allocatable  :: evnm,evnmb,evnmr
   real(kind=4), dimension(:), allocatable       :: evla,evlo,evdp,evbur
   real(kind=4)                                  :: dummy
   integer, dimension(:), allocatable            :: bk_prev,rk_prev
   ! sampling
   real(kind=4)                                  :: r1 ! random number (0,1] for choosing the first candidate
   real(kind=4), dimension(:), allocatable       :: rr ! m*n candidates, (0,1]
   integer, dimension(:), allocatable            :: rk ! m*n candidates
   integer, dimension(:), allocatable            :: bk ! mini-batch series
   !real(kind=4)                                  :: az, baz, dist_deg, dist_km
   real(kind=4)                                  ::dist_min, dist_max
   integer                                       :: is_bk, m, rk0
   double precision                              :: eqlt,eqln,stlt,stln
   ! double precision                              :: delta,deltdg,deltkm,azeqst,azesdg,azsteq,azsedg
   double precision                              :: delta,deltkm,azeqst,azsteq
      
   narg=command_argument_count()
   if ( narg .ne. 2 ) then
     print *,'Usage: sampling_minibatch model nsrc'
     print *,'model     - current model name (iteration) ' 
     print *,'nsrc      - number of sources in the mini-batch ' 
     stop
   endif
   call get_command_argument(1,model) 
   call get_command_argument(2,nk_str) 
   read(nk_str,*) nk
   allocate(bk(nk))

   read(model(2:),'(I2.2)') imod_current
   imod_up=imod_current+1
   imod_down=imod_current-1
   write(model_next,'(A1,I2.2)') 'M',imod_up
   if (imod_down <0) then
      model_prev='none'
   else
      write(model_prev,'(A1,I2.2)') 'M',imod_down
   endif

   ! read source list
   infname='src_rec/sources.dat'
   open(9,file=trim(infname),status='old')
   nsrc=0
   do
     read(9,*,iostat=ier) kname,lat,lon,ele,bur
     !write(*,*) trim(kname),lat,lon,ele,bur
     if (ier < 0 ) exit
     nsrc=nsrc+1
   enddo
   close(9)

   print *,'infname,nsrc,nk=',trim(infname),nsrc,nk
   allocate(evnm(nsrc))
   allocate(evla(nsrc))
   allocate(evlo(nsrc))
   allocate(evdp(nsrc))
   allocate(evbur(nsrc))

   open(9,file=trim(infname),status='old')
   i=1
   do
     read(9,*,iostat=ier) evnm(i),evla(i),evlo(i),evdp(i),evbur(i)
     if (ier < 0 ) exit
     !write(*,*) trim(evnm(i)),evla(i),evlo(i),evdp(i),evbur(i)
     i=i+1
   enddo
   close(9)
  
   if (imod_current==0) then
      ! choose the first source randomly
      call random_seed()
      call random_number(r1)
      bk(1)=floor(r1*(nsrc-1))+1
      print *,'Choose the first source randomly'
      print *,'  Random number: ', r1
      print *,'  The first source number: ', bk(1) 
      rk0=2 ! 
   else ! M01, M02, ...
      !infname='src_rec/sources_batch.'//trim(model_prev)//'.dat'
      call system('cat src_rec/sources_batch.M??.dat |sort -n |uniq >src_rec/sources_batch.used.dat')
      infname='src_rec/sources_batch.used.dat'
      open(9,file=trim(infname),status='old')
       nsrcb=0
      do
        read(9,*,iostat=ier) kname,lat,lon,ele,bur
        if (ier < 0 ) exit
        nsrcb=nsrcb+1
      enddo
      close(9)

      print *,'infname,nsrcb,nk=',trim(infname),nsrcb,nk
      allocate(evnmb(nsrcb))
      allocate(bk_prev(nsrcb))

      open(9,file=trim(infname),status='old')
      do i=1,nsrcb
        read(9,*,iostat=ier) evnmb(i),dummy,dummy,dummy,dummy
        do j=1,nsrc
           if (trim(evnmb(i))==trim(evnm(j))) then
              bk_prev(i)=j
           endif
        enddo
        write(*,*) trim(evnmb(i)),bk_prev(i)
      enddo
      close(9)

      infname='src_rec/sources_ctrlgrp.'//trim(model_prev)//'.dat'
      open(9,file=trim(infname),status='old')
       nsrcr=0
      do
        read(9,*,iostat=ier) kname,lat,lon,ele,bur
        if (ier < 0 ) exit
        nsrcr=nsrcr+1
      enddo
      close(9)

      print *,'infname,nsrcr,nk=',trim(infname),nsrcr,nk
      if (nk <2*nsrcr) then
         stop 'Batch size is less then 2*ctrlgrp size'
      endif
      allocate(evnmr(nsrcr))
      allocate(rk_prev(nsrcr))

      open(9,file=trim(infname),status='old')
      do i=1,nsrcr
        read(9,*,iostat=ier) evnmr(i),dummy,dummy,dummy,dummy
        do j=1,nsrc
           if (trim(evnmr(i))==trim(evnm(j))) then
              rk_prev(i)=j
           endif
        enddo
        write(*,*) trim(evnmr(i)),rk_prev(i)
      enddo
      close(9)
      ! set rki_prev to the initial bk
      do i=1,nsrcr
         bk(i)=rk_prev(i)
      enddo
      rk0=nsrcr+1
   endif ! if imod_current==0
   

   m=10 ! constant ratio to scale up of number of random candidate points, m*n
   do k=rk0, nk 
     allocate(rr(m*(k-1)+1))
     allocate(rk(m*(k-1)+1))
     call random_seed()
     call random_number(rr)
     rk=floor(rr*(nsrc-1))+1
     !print *,'k,rk=',k,rk
     dist_max=0.
     do i=1, m*(k-1)+1
       ! if i has been selected
       is_bk=0
       do j=1,k-1 ! loop current bk
          if (rk(i)==bk(j)) then
            is_bk=1
            exit
          endif
       enddo
       if (imod_current>0) then
         do j=1,nsrcb ! loop previous bk
            if (rk(i)==bk_prev(j)) then
              is_bk=1
              exit
            endif
         enddo
       endif
       
       if (is_bk==0) then
         !write(*,*) 'Not selected before, try to find a new one ...'
         dist_min=10000000000.
         do j=1,k-1
            !write(*,*) 'Test=',j,evla(bk(j)),evlo(bk(j)),evla(rk(i)),evlo(rk(i))
            eqlt=evla(bk(j))
            eqln=evlo(bk(j))
            stlt=evla(rk(i))
            stln=evlo(rk(i))
            ! call delaz(eqlt,eqln,stlt,stln,delta,deltdg,deltkm,azeqst,azesdg,azsteq,azsedg,0)
            ! ===Dec 2022===: Mijan modify distaz here. Tests are required.
            call distaz(stlt, stln, eqlt, eqln, azeqst, azsteq, delta, deltkm)
            !write(*,*) 'DELAZ: deltdg,deltkm=',azesdg,azsedg,deltdg,deltkm
            ! WK: distaz sometime gets float error when the lat/lon of EQ and ST are very close
            !call distaz(evla(bk(j)),evlo(bk(j)),evla(rk(i)),evlo(rk(i)),az,baz,dist_deg,dist_km)  
            !write(*,*) 'distaz: dist_deg,dist_km=',az,baz,dist_deg,dist_km
            !if (dist_km<dist_min) dist_min=dist_km
            if (deltkm<dist_min) dist_min=deltkm
         enddo
         if (dist_min>dist_max) then
            dist_max=dist_min
            bk(k)=rk(i)
         endif
       endif
     enddo
     print *,'k,bk(k),dist_max= ',k,bk(k),dist_max
     deallocate(rr)
     deallocate(rk)
   enddo 

   ! write mini-batch
   outfname='src_rec/sources_batch.'//trim(model)//'.dat'
   open(9,file=trim(outfname))
   do k=1,nk
     write(9,"(A16,4F22.4)") evnm(bk(k)),evla(bk(k)),evlo(bk(k)),evdp(bk(k)),evbur(bk(k))
   enddo
   close(9)
    
   deallocate(evnm)
   if (imod_current>0) then
      deallocate(evnmb,evnmr,bk_prev,rk_prev) 
   endif
   deallocate(evla) 
   deallocate(evlo) 
   deallocate(evdp) 
   deallocate(evbur) 


end program dynmbat_sele_minibatch

