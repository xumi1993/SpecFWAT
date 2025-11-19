module decon_mod
  use fftpack, only: fft_cls
    
contains
  subroutine deconvolve(utr, wtr, dt, tshift, f0, &
                       maxit, minderr, ipart, rfi, use_gpu)
    ! -----------------------------------------
    ! Time iterative deconvolution for receiver 
    ! function calculation. Code refer to that in
    ! CPS330.
    ! 
    ! Nov 19, 2025 Mijian Xu @ University of Toronto
    ! ------------------------------------------- 

    implicit none
    double precision, dimension(:), intent(in)                  :: utr, wtr
    real, intent(in)                                            :: dt,tshift,f0,minderr
    integer, intent(in)                                         :: ipart,maxit
    double precision, dimension(:), allocatable, intent(out)    :: rfi
    logical, intent(in)                                         :: use_gpu

    if (use_gpu) then
      call deconit_gpu(utr, wtr, dt, tshift, f0, maxit, minderr, ipart, rfi)
    else 
      call deconit(utr, wtr, dt, tshift, f0, maxit, minderr, ipart, rfi)
    endif

  end subroutine deconvolve

  subroutine deconit_gpu(utr, wtr, dt, tshift, f0, &
                     maxit, minderr, ipart, rfi)
    implicit none
    real, intent(in)                                            :: dt,tshift,f0,minderr
    integer, intent(in)                                         :: ipart,maxit
    double precision, dimension(:), intent(in)                  :: utr, wtr
    double precision, dimension(:), allocatable, intent(out)    :: rfi
    integer                                                     :: nft, nt

    interface
      subroutine deconit_cuda(utr, wtr, nt, nft, dt, tshift, f0, &
                              maxit, minderr, ipart, rfi) bind(C, name="deconit_cuda")
        use iso_c_binding
        real(c_double), dimension(*), intent(in) :: utr, wtr
        integer(c_int), value :: nt, nft
        real(c_float), value :: dt, tshift, f0
        integer(c_int), value :: maxit
        real(c_float), value :: minderr
        integer(c_int), value :: ipart
        real(c_double), dimension(*), intent(out) :: rfi
      end subroutine deconit_cuda
    end interface

    nt = size(utr)
    nft=nt
    call npow2(nft)
    
    if (allocated(rfi)) deallocate(rfi)
    allocate(rfi(nt))
    
    call deconit_cuda(utr, wtr, nt, nft, dt, tshift, f0, maxit, minderr, ipart, rfi)

  end subroutine deconit_gpu

  subroutine deconit(utr, wtr, dt, tshift, f0, &
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
    integer, intent(in)                                         :: ipart,maxit
    double precision, dimension(:), intent(in)                  :: utr, wtr
    double precision, dimension(:), allocatable, intent(out)    :: rfi
    integer                                                     :: nft,i,maxlag,it, nt
    integer, dimension(1)                                       :: i1
    complex, dimension(:), allocatable                          :: wf, gauss
    double precision, dimension(:), allocatable                 :: uflt,wflt, rflt, rw, p0, pflt,rms
    double precision                                            :: powerU, sumsq_i, d_error, amp, sumsq

    nt = size(utr)
    nft=nt
    call npow2(nft)
    allocate(gauss(nft), uflt(nft), wflt(nft), rflt(nft), wf(nft),&
             rw(nft), p0(nft), pflt(nft), rms(maxit))
    gauss(:) = 0.
    uflt(:) = 0.
    wflt(:) = 0.
    rflt(:) = 0.
    wf(:) = 0.
    rw(:) = 0.
    pflt(:) = 0.
    p0(:) = 0.
    rms(:) = 0.

    call gaussfilter(dt, nft, f0, gauss)
    
    uflt(1:nt) = utr
    wflt(1:nt) = wtr

    wf = fftcls%fft(wflt, nft)
    
    call filter(uflt, nft, dt, gauss)
    call filter(wflt, nft, dt, gauss)

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
    ! write(*, *) 'Iter=',it,'; RMS=',sumsq
    call filter(p0, nft, dt, gauss)
    call phaseshift(p0, nft, dt, tshift)
    ! call fft%Destory()
    rfi = p0(1:nt)

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
        gauss1(i) = exp(-0.25 * (2 * 3.1415926535 * df * (i-1) / f0) ** 2)/dt
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
    integer :: nsamp, npts
    nsamp = npts
    npts = 1
    
    do
        npts = 2 * npts
        if (npts >= nsamp) exit
    end do
    
    return
        
  end subroutine npow2
end module decon_mod
