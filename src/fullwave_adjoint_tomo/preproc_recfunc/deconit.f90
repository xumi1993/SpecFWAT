module decon_mod
    use fftpack, only: fft_cls
    
contains 
subroutine deconit(utr, wtr, nt, dt, tshift, f0, &
                   maxit, minderr, ipart, rfi)
    ! -----------------------------------------
    ! Time iterative deconvolution for receiver 
    ! function calculation. Code refer to that in
    ! CPS330.
    ! 
    ! Dec 31, 2021 Mijian Xu @ Nanyang Technological University
    ! ------------------------------------------- 

    implicit none
    type(fft_cls)                                               :: fftcls
    real, intent(in)                                            :: dt,tshift,f0,minderr
    integer, intent(in)                                         :: nt,ipart,maxit
    double precision, dimension(nt), intent(in)                 :: utr, wtr
    double precision, dimension(nt), intent(inout)              :: rfi
    integer                                                     :: nft,i,maxlag,it
    integer, dimension(1)                                       :: i1
    complex, dimension(:), allocatable                          :: wf, gauss
    double precision, dimension(:), allocatable                 :: uflt,wflt, rflt, rw, p0, pflt
    double precision                                            :: powerU, sumsq_i, d_error, amp, sumsq
    double precision, dimension(:), allocatable                 :: rms

    nft=nt
    call npow2(nft)
    if (allocated(gauss)) deallocate(gauss)
    allocate(gauss(nft))
    if (allocated(uflt)) deallocate(uflt)
    allocate(uflt(nft))
    if (allocated(wflt)) deallocate(wflt)
    allocate(wflt(nft))
    if (allocated(rflt)) deallocate(rflt)
    allocate(rflt(nft))
    if (allocated(wf)) deallocate(wf)
    allocate(wf(nft))
    if (allocated(rw)) deallocate(rw)
    allocate(rw(nft))
    if (allocated(p0)) deallocate(p0)
    allocate(p0(nft))
    if (allocated(pflt)) deallocate(pflt)
    allocate(pflt(nft))
    if (allocated(rms)) deallocate(rms)
    allocate(rms(maxit))
    gauss(:) = 0.
    uflt(:) = 0.
    wflt(:) = 0.
    rflt(:) = 0.
    wf(:) = 0.
    rw(:) = 0.
    pflt(:) = 0.
    p0(:) = 0.
    rms(:) = 0.

    ! call FFT%Create(nft)
    call gaussfilter(dt, nft, f0, gauss)
    
    uflt(1:nt) = utr
    wflt(1:nt) = wtr

    wf = fftcls%fft(wflt, nft)
    
    call filter(uflt, nft, dt, gauss)
    call filter(wflt, nft, dt, gauss)

    ! wf = FFT%Forward(real(wflt))
    rflt = uflt
    powerU = sum(uflt ** 2)

    sumsq_i = 1.0
    d_error = 100 * powerU + minderr
    if (ipart == 0) then
        maxlag = 0.5 * nft
    else
        maxlag = nft
    endif

    i=0
    do while(abs(d_error) > minderr .and. i < maxit)
        i = i+1
        call correl(rflt, wflt, nft, rw)
        rw = rw / sum(wflt ** 2)
        i1 = maxloc(abs(rw(1:maxlag)))
        amp = rw(i1(1)) / dt
        p0(i1(1)) = p0(i1(1)) + amp
        pflt = p0
        call filter(pflt, nft, dt, gauss)
        call filter(pflt, nft, dt, wf)
        rflt = uflt - pflt
        sumsq = sum(rflt ** 2) / powerU
        rms(i) = sumsq
        d_error = 100 * (sumsq_i - sumsq)
        sumsq_i = sumsq
    enddo
    it = i
    write(*, *) 'Iter=',it,'; RMS=',sumsq
    ! call dwsac1(trim(fname)//'z.sac',wflt,nft,dble(0.),dble(dt))
    call filter(p0, nft, dt, gauss)
    call phaseshift(p0, nft, dt, tshift)
    ! call fft%Destory()
    rfi = p0(1:nt)

    deallocate(gauss)
    deallocate(uflt)
    deallocate(wflt)
    deallocate(rflt)
    deallocate(wf)
    deallocate(rw)
    deallocate(p0)
    deallocate(pflt)
    deallocate(rms)

end subroutine deconit

subroutine filter(x, nft, dt, gauss)
    implicit  none
    type(fft_cls) :: fftcls
    integer :: nft
    real :: dt
    ! double precision, dimension(nft), intent(in):: gauss
    double precision, dimension(nft), intent(inout):: x
    complex, dimension(nft) :: z1, gauss
    ! z1 = cmplx(x, 0.0)

    z1 = fftcls%fft(x, nft)
    z1 = z1 * gauss *dt
    ! z1 = z1/nft
    x = fftcls%ifft(z1, nft)

end subroutine filter


subroutine gaussfilter(dt, nft, f0, gauss)
    implicit none
    integer nft, nft21, i
    real :: dt, f0, df
    complex, dimension(nft), intent(out) :: gauss
    double precision, dimension(nft) :: gauss1

    gauss1(:) = 0.
    gauss =  cmplx(0., 0.)
    df = 1.0 / (nft * dt)
    nft21 = 0.5 * nft + 1
    do i = 1,nft21
        gauss1(i) = exp(-0.25 * (2 * 3.14159265 * df * (i-1) / f0) ** 2)/dt
        ! print *, gauss(i) 
    enddo
    do i=nft21+1,nft
        gauss1(i) = gauss1(2*nft21-i)
    enddo
    gauss = cmplx(gauss1)

end subroutine gaussfilter


subroutine correl(r, w, nft, rw)
    implicit none
    type(fft_cls) :: fftcls
    ! Type(CLS_FFT) :: fftcls 
    integer, intent(in) :: nft
    complex, dimension(nft) :: r1, w1, dat
    double precision, dimension(nft), intent(in):: r, w
    double precision, dimension(nft), intent(out):: rw
    ! real, dimension()

    r1 = fftcls%fft(r, nft)
    w1 = fftcls%fft(w, nft)

    dat = r1* conjg(w1)
    ! dat = dat / nft
    rw = fftcls%ifft(dat, nft)
       
end subroutine correl


subroutine phaseshift(x, nft, dt, tshift)
    implicit none
    type(fft_cls) :: fftcls
    integer :: nft, shift_i, i
    real :: dt, tshift, p
    double precision, dimension(nft), intent(inout) :: x
    complex, dimension(nft) :: z1

    z1 = fftcls%fft(x, nft)
    shift_i=int((tshift/dt)+0.5)
    do i=1,nft
        p = 2*3.14159265*(i-1)*shift_i/nft
        z1(i) = z1(i) * cmplx(cos(p), -sin(p))
    enddo
    ! z1 = z1/nft
    x = fftcls%ifft(z1, nft)

end subroutine phaseshift

subroutine npow2(npts)
!-----
!      Given npts, determine the N=2**m such that N >= npts
!      return the new ntps
!-----
        integer nsamp, npts
        nsamp = npts
        npts = 1
 1000   continue
        npts = 2*npts
        if(npts.lt.nsamp)goto 1000
        return
        
end subroutine npow2
end module decon_mod
