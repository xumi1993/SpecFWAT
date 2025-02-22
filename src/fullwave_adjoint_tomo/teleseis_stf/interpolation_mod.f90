module interpolation_mod

  !use constants, only: CUSTOM_REAL

  implicit none

  !**************************************************
  !* CONTENT OF MY MODULE PRECISION_MOD
  !*** Single precision
  integer, parameter :: si = selected_int_kind(8)
  integer, parameter :: sp = selected_real_kind(4)

  !*** Double precision
  integer, parameter :: di = selected_int_kind(16)
  integer, parameter :: dp = selected_real_kind(8)

  !*** Custom precision
  integer, parameter :: cp = selected_real_kind(4)
  !integer, parameter :: cp = CUSTOM_REAL

  !*** Special precision
  integer, parameter :: hp = selected_real_kind(8)
  !***************************************************

  !***************************************************
  !* CONTENT OF MY MODULE CONSTANTS_MOD
  !*** Constants for inititalization
  real(kind=cp), parameter :: mytiny = 1.e-9_cp
  real(kind=cp), parameter :: myhuge = 1.e30_cp

  real(kind=cp),    parameter :: myverysmall = 1.e-25_cp
  real(kind=dp),    parameter :: mysmall     = 0.000001_dp
  integer(kind=si), parameter :: myhugeint   = 100000000

  real(kind=cp), parameter :: zero = 0._cp
  real(kind=cp), parameter :: one  = 1._cp
  real(kind=cp), parameter :: half = 0.5_cp

  real(kind=cp), parameter :: onethird  = 1._cp/3._cp
  real(kind=cp), parameter :: twothird  = 2._cp/3._cp
  real(kind=cp), parameter :: fourthird = 4._cp/3._cp

  real(kind=cp), parameter :: onefourth = 1._cp/4._cp
  real(kind=cp), parameter :: onesixth  = 1._cp/6._cp

  !*** Constants
  real(kind=cp), parameter :: pi  = 3.141592653589793_cp
  real(kind=cp), parameter :: twopi = 2._cp * pi

  !*** Constants for conversions
  real(kind=cp), parameter :: deg2rad = 0.017453292519943
  real(kind=cp), parameter :: rad2deg = 57.295779513082323

  !*** Check stability in modeling
  real(kind=cp), parameter :: stability_criterion = 1.e28_cp
  !**********************************************************


contains


!================================================================================
! Filtering routines
  subroutine bwfilt (x, y, dt, n, irek, norder, f1, f2)

    ! recursive filtering of data with butterworth filter
    ! x: input array
    ! y: output array
    ! dt: time increment
    ! n: number of data points

    ! irek=0: forward filtering only
    ! irek=1: forward and backward filtering

    ! norder: order of butterworth filter
    ! norder=0: only filtering, no determination of coefficients
    ! norder < 0: no starplots of transfer function and impulse response

    ! f1: low cutoff frequency (Hz)
    ! f1=0: low pass filter

    ! f2: high cutoff frequency (Hz)
    ! f2>0.5/dt: high pass filter

    implicit none

    integer(kind=si) :: iunit, npoles,n,lx
    integer(kind=si) :: irek,norder
    real(kind=cp), dimension(n) ::x,y
    real(kind=cp), dimension (10) ::  a, b1, b2
    real(kind=cp) :: dt,f1,f2

    !real(kind(0d0)) :: x(n),y(n)

    iunit = 3

    if (norder /= 0) then
       npoles=abs(norder)
       !determination of filter coefficients
       call bpcoeff(f1,f2,npoles, dt, a,b1, b2)
       if (norder >= 0) then
          !plot of transfer function and impuulse response
          lx = 100
          !filtering
       endif
    endif


    if (n /= 0) then
       call rekurs(x,y,n,a,b1,b2,npoles,irek)
    endif
    return
  end subroutine bwfilt

  subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
    ! performs recursive filtering of data in array x of length ndat
    ! filtered output in y
    ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
    ! npoles is the number of poles
    ! iflag=0: forward filtering only
    ! iflag /= 0: forward and backward filtering

    implicit none

    real(kind=cp), dimension(10) :: z,z1,z2 ,a,b1,b2
    real(kind=cp)  ::  x1,x2
    integer(kind=si) :: ndat, npoles,n,i
    integer(kind=si) :: iflag
    real(kind=cp), dimension(ndat) :: x, y

    !forward

    x1 = 0.d0
    x2 = 0.d0

    do i = 1, npoles
       z1(i) = 0.d0
       z2(i) = 0.d0
    enddo

    do n = 1, ndat
       z(1) = a(1)*(x(n)-x2) -b1(1)*z1(1) -b2(1)*z2(1)
       do i = 2, npoles
          z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
       enddo
       x2=x1
       x1=x(n)
       do i = 1, npoles
          z2(i) =z1(i)
          z1(i) =z(i)
       enddo
       y(n) = z(npoles)
    enddo

    if (iflag == 0) then
       return
    endif

    !backward

    x1 =0.d0
    x2 =0.d0

    do i = 1, npoles
       z1(i) = 0.d0
       z2(i) = 0.d0
    enddo

    do n = ndat, 1, -1
       z(1) = a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
       do i =2, npoles
          z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
       enddo
       x2=x1
       x1=y(n)
       do i = 1,npoles
          z2(i)=z1(i)
          z1(i)=z(i)
       enddo
       y(n) = z(npoles)
    enddo
    return
  end subroutine rekurs

  subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
    !determines filtercoefficients for recursive bandpassfilter

    real(kind=cp),dimension(10) :: a,b1,b2
    complex(kind=4) :: s(20), t1,t2,p
    real(kind=cp), parameter :: pi = 3.141592653589793
    real(kind=cp) :: f1,f2,dt,d2,w0,w1,w2,ssum, sprod,fact1,fact2,fact3
    integer(kind=si) :: i,npol2,n,npoles


    if (npoles > 10) then
       stop ' npoles greater than 10: STOP '
    endif

    d2= 2/dt
    w1=d2*tan(2.*pi*f1/d2)
    w2=d2*tan(2.*pi*f2/d2)
    w0=0.5*(w2-w1)

    i=1
    npol2=npoles/2+1
    do n =1,npoles
       p = cexp(cmplx(0.,real(2*n-1+npoles)*pi/real(2*npoles)))
       t1 = p*cmplx(w0,0.)
       t2 = sqrt(t1*t1-cmplx(w1*w2,0.))
       s(i)=t1+t2
       s(i+1)=t1-t2
       i=i+2
    enddo

    do n=1,npoles
       ssum=2*real(s(n))
       sprod=real(s(n)*conjg(s(n)))
       fact1=d2*d2-d2*ssum+sprod
       fact2=2.*(sprod-d2*d2)
       fact3=d2*d2+d2*ssum+sprod
       a(n)=2.*d2*w0/fact1
       b1(n)=fact2/fact1
       b2(n)=fact3/fact1
    enddo
    return
  end subroutine bpcoeff
!--------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! taper in traces form index i0,..., i3
!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine taper_window_W(t_w,i0,i1,i2,i3,nstep,W)

    implicit none

    integer,                                             intent(in)    :: i0, i1, i2, i3, nstep
    real(kind=cp),                              intent(in)    :: W
    real(kind=cp), dimension(:), allocatable,   intent(inout) :: t_w

    integer                                                   :: i
    real(kind=cp)                                             :: omega, phi, pi

    PI = 3.1415926d0
    t_w(1:nstep)=0.

    ! take off
    omega = pi / (i1 - i0)
    phi = pi / 2 - omega * i1
    do i = i0, i1
       t_w(i) = W*(0.5 + 0.5 *sin(omega * i + phi))
    enddo

    ! flying
    do i = i1+1,i2-1
       t_w(i)=W
    enddo

    ! landing
    omega = pi / (i3 - i2)
    phi = pi/2 - omega * i2
    do i= i2,i3
       t_w(i) = W*(0.5 + 0.5 * sin(omega * i + phi))
    enddo

  end subroutine taper_window_W


!================================================================================
! Iterative time domain deconvolution
  subroutine time_deconv(dobs,dcal,dt,nt,nit,src_sum)

    integer(kind=si), intent(in) :: nt, nit
    integer(kind=si)             :: i, ii
    integer(kind=si), parameter  :: part=0

    real(kind=cp), intent(in) :: dt
    real(kind=cp)             :: amp_corr

    real(kind=cp), dimension(nt), intent(in)    :: dobs, dcal
    real(kind=cp), dimension(:), allocatable, intent(inout) :: src_sum

    real(kind=cp), dimension(:), allocatable :: dobs2, autocorr, crosscorr, new_obs, src_one


    !* 0. Alloctae
    if (.not. allocated(autocorr))  allocate(autocorr(nt))
    if (.not. allocated(crosscorr)) allocate(crosscorr(nt))
    if (.not. allocated(src_sum))   allocate(src_sum(nt))
    if (.not. allocated(src_one))   allocate(src_one(nt))
    if (.not. allocated(dobs2))     allocate(dobs2(nt))
    if (.not. allocated(new_obs))   allocate(new_obs(nt))
    src_sum = 0._cp
    src_one = 0._cp

    !* 1. Estimate auto-corr of computed data
    call mycorrelation(dcal,dcal,nt,nt,autocorr,part)
    autocorr = autocorr * dt
    amp_corr = 1._cp / maxval(abs(autocorr))
   !  write(*,*) 'amp_corr',amp_corr

    !* 2. Time iteration to estimate sources
    dobs2 = dobs
    do i = 1, nit

       !* Compute correlation
       call mycorrelation(dobs2,dcal,nt,nt,crosscorr,part)
       crosscorr = crosscorr * dt

       !* Find maximum of correlation
       ii = maxloc(abs(crosscorr),dim=1)

       !* Put local contibution to src_one and src_sum
       src_one     = 0._cp
       src_one(ii) = crosscorr(ii) * amp_corr
       src_sum     = src_sum + src_one

       !* Convolve
       call myconvolution(src_one,dcal,nt,nt,new_obs,part)
       new_obs = new_obs * dt
       dobs2 = dobs2 - new_obs

    enddo

  end subroutine time_deconv

!--------------------------------------------------------------------------------

!================================================================================
! Convolution routine
   subroutine myconvolution(sig1,sig2,n1,n2,conv,part)

     integer(kind=si), intent(in) :: n1, n2, part

     real(kind=cp), dimension(n1), intent(in) :: sig1
     real(kind=cp), dimension(n2), intent(in) :: sig2

!     real(kind=cp), dimension(n1), intent(out) ::conv
      real(kind=cp), dimension(:), allocatable, intent(inout) ::conv

     real(kind=cp), dimension(n1+n2-1) :: convtmp !, intent(out) :: conv

     integer(kind=si) :: i1, i2, ind

     !*** Put to zero
     convtmp = zero

     !*** Convolve
     do i1=1,n1
        do i2=1,n2
           convtmp(i1+i2-1) = convtmp(i1+i2-1) + sig1(i1) * sig2(i2)
        enddo
     enddo

     if (part == 0) then !(middle regular)
        !*** Take good parts (this is wrong...)
        if (modulo(n2,2) == 0) then
           ind = n2/2+1
        else
           ind = ceiling(real(n2/2,kind=cp))
        endif
        if (.not. allocated(conv)) allocate(conv(n2))
        conv(1:n2) = convtmp(ind:ind+n2-1)

     else if (part == 1) then ! full convolution

        if (.not. allocated(conv)) allocate(conv(n1+n2-1))
        conv(:) = convtmp(:)

     else if (part == 2) then !(middle irregular)
        if (.not. allocated(conv)) allocate(conv(n2))
        conv(1:n2) = convtmp(n2:n2+n2-1)
     endif

   end subroutine myconvolution
!--------------------------------------------------------------------------------

!================================================================================
! Correlation routine
   subroutine mycorrelation(sig1,sig2,n1,n2,corr,part)

     integer(kind=si), intent(in) :: n1, n2, part

     real(kind=cp), dimension(n1), intent(in) :: sig1
     real(kind=cp), dimension(n2), intent(in) :: sig2

!     real(kind=cp), dimension(n1), intent(out) :: corr
     real(kind=cp), dimension(:), allocatable, intent(inout) :: corr

     real(kind=cp), dimension(n2) :: flipsig2
     integer(kind=si) :: i

     !*** Choose size of corr
     if (part == 0) then !(middle)
        if (.not. allocated(corr)) allocate(corr(n2))
     else if (part == 1) then !(full)
        if (.not. allocated(corr)) allocate(corr(n1+n2-1))
     endif

     !*** Flip second signal
     do i=1,n2
        flipsig2(i) = sig2(n2-i+1)
     enddo

     !*** Use regular convolution
     call myconvolution(sig1,flipsig2,n1,n2,corr,part)

   end subroutine mycorrelation
!--------------------------------------------------------------------------------
!================================================================================
!         Added by Kai Wang from CPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine zfour(zarr,nn,isign,dt,df) 
!-----
!     THE input is a complex array
!     which has numbers stored in memory as
!     R1, I1, R2, I2, ..., Rnn, Inn
!     where nn must be a power of 2 R and I are the real and imaginary
!     parts of the complex number
!
!     For isign -1 this is a complex time series
!     For isign +1 this is a complex frequency series with
!        index 1 (in fortran corresponding to f=0
!              2                              f=df
!            nn/2 + 1                         f = 1/2dt = Nyquist
!            nn - 1                           f = -2df
!             nn                              f = -df

!-----
!     the cooley-tookey fast fourier transform in usasi basic fortran
!     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
!     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
!     dimensional complex array (i.e., the real and imaginary parts of
!     datc are located immediately adjacent in storage, such as fortran
!     places them) whose length nn is a power of two.  isign
!     is +1 or -1, giving the sign of the transform.  transform values
!     are returned in array datc, replacing the input datc.  the time is
!     proportional to n*log2(n), rather than the usual n**2
!     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
!     b is the number of bits in the floating point fraction.
!
!     the program computes df from dt, dt from df and checks to see
!     if they are consistent. In addition, the transforms are multiplied
!     by dt or df to make the results dimensionally correct
!
!     This is a slightly modified version of the original Brenner routine
!     The changes were to add physical dimensions to the transform
!     and to make it all complex
!-----
        complex zarr(*) 
        integer ii,jj,mm,nn,mmax,istep, isign
        real dt, df
        real theta,sinth,wstpr,wstpi,wr,wi,tempr

        complex ztemp
!-----
!       ensure that the dt and df are defined and
!       consistent
!-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
!-----
!       now begin the transform
!-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
!-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
!-----
!       use trig relations to compute the next sin/cos
!       without actually calling sin() or cos()
!-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
!-----
!       transform is done
!-----
   10   continue 
!-----
!     give the arrays the proper physical dimensions
!-----
        if(isign.lt.0)then
!-----
!             time to frequency domain
!-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
!-----
!             frequency to time domain
!-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
        endif
        return
        end subroutine zfour

        subroutine myzero(x,n)
!-----
!      zero a real array
!-----
!      x   R*4 - array to be zeroed
!      n   I*4 - number of points in array
!-----
        !real x(n)
        integer n
        real(kind=4), dimension(n) :: x
        integer i
        
        do 1 i = 1,n
            x(i) = 0
    1   continue
        return
        end subroutine myzero



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

        subroutine gfilter(x,gwidth,n,nn,dt)
!c-----
!c      convolve a function with a unit-area Gaussian filter
!c-----
!c      x   R*4 - array to be filtered
!c      Gwidth
!c          R*4 - Filter is exp( - ( pi freq/gwidth)**2)
!c      n   I*4 - number of points in time series
!c      nn  I*4 - number of points at next largest power of 2
!c      dt  R*4 - sampling interval
!c-----
        !real x(n)
        integer n
        real(kind=4), dimension(n) :: x
        real gwidth

        !parameter (NSAMP=131072)
        integer(kind=4), parameter :: NSAMP=131072
        common/zval/z1, z2
        complex z1(NSAMP), z2(NSAMP)

        real :: fac,df,freq,dt
        integer :: i, n21,nn

        do 1000 i = 1, nn
            if(i.le.n)then
                z1(i) = cmplx(x(i),0.0)
            else
                z1(i) = cmplx(0.0,0.0)
            endif
 1000   continue
!c-----
!c      get Fourier transform
!c-----
!        call zfour(z1,nn,-1,dt,df)
!c-----
!c      Gaussian filter
!c-----
        n21 = nn / 2 + 1
        do 2000 i=1, n21
            freq = ( i - 1) * df
            fac = 3.1415927*freq/gwidth
            fac = fac * fac
            if(fac.gt.50.0)then
                fac = 0.0
            else
                fac = exp( -fac)
            endif
            z1(i) = z1(i) * fac
            if(i.gt.1)then
                z1(nn + 2 - i) = conjg(z1(i))
            endif
 2000   continue
!c-----
!c      ensure Nyquist frequency element is real
!c-----
!        z1(n21) = cmplx(real(z1(n21)),0.0)
!c-----
!c      get inverse Fourier transform
!c-----
        call zfour(z1,nn,+1,dt,df)
!c-----
!c      reconstitute the series
!c-----
        do 3000 i=1,nn
            x(i) = real(z1(i))
 3000   continue
        return
        end subroutine gfilter

        subroutine phs_shift(x,theshift,n,nn,dt)
!c-----
!c      time shift a signal
!c------
!c      X   R*4 - signal to be shifted
!c      theshift R*4    - time shift in seconds
!c      n   I*4 - length of signal
!c      nn  I*4 - length of signal to nearest power of 2
!c      dt  R*4 - sampling interval
!c-----
        integer n
        real(kind=4), dimension(n) :: x
        real pi, two_pi, theshift, d_omega
        real dt, df, ang, c, s
        integer i, nn, n21
        integer forward, inverse

        !parameter (NSAMP=131072)
        integer(kind=4), parameter :: NSAMP=131072
        common/zval/z1, z2
        complex z1(NSAMP), z2(NSAMP)
!c     
        forward = 1
        inverse = -1
        pi = 3.1415927
        two_pi = 2 * pi
        df=0.0

        do 1000 i=1,nn
            if(i.le.n)then
                z1(i) = cmplx(x(i),0.0)
            else
                z1(i) = cmplx(0.0,0.0)
            endif
 1000   continue
!c-----
!c      get the Fourier transform
!c-----
        call zfour(z1,nn,-1,dt,df)
        d_omega = two_pi * df
!c-----
!c      time shift
!c-----
        n21 = nn / 2 + 1
        do 2000 i=1,n21
            ang = (i-1)*d_omega*theshift
            c = cos(ang)
            s = sin(ang)
            z1(i) = z1(i) * cmplx(c,-s)
            if(i.gt.1)then
                z1(nn + 2 - i) = conjg(z1(i))
            endif
 2000   continue
!c-----
!c      ensure Nyquist frequency element is real
!c-----
        z1(n21) = cmplx(real(z1(n21)),0.0)
!c-----
!c      inverse Fourier transform to the time domain
!c-----
        call zfour(z1,nn,+1,dt,df)
!c-----
!c      redonstitute the time series
!c-----
        do 3000 i=1,nn
            x(i) = real(z1(i))
 3000   continue
        return
        end subroutine phs_shift
             

        subroutine convolve(x,y,n,nn,dt)
!c-----
!c      convolve x and y, replacing the x array
!c-----
!c      x   R*4 - array
!c      y   R*4 - array
!c      n   I*4 - number of points in time series
!c      nn  I*4 - number of points rounded up to next power of 2
!c      dt  R*4 - sampling interval
!c-----
        real  dt, df
        integer n
        real(kind=4), dimension(n) :: x
        real(kind=4), dimension(n) :: y
        integer i, nn, n21

        !parameter (NSAMP=131072)
        integer(kind=4), parameter :: NSAMP=131072
        common/zval/z1, z2
        complex z1(NSAMP), z2(NSAMP)

!c------
!c      convolve  by frequency domain multiplication
!c-----
        do 1000 i=1,nn
            if(i.le.n)then
                z1(i) = cmplx(x(i),0.0)
                z2(i) = cmplx(y(i),0.0)
            else
                z1(i) = cmplx(0.0,0.0)
                z2(i) = cmplx(0.0,0.0)
            endif
 1000   continue
!c-----
!c      get Fourier transforms
!c-----
        df=0.0
        call zfour(z1,nn,-1,dt,df)
        call zfour(z2,nn,-1,dt,df)
!c-----
!c      convolution = F  G
!c-----
        n21 = nn / 2 + 1
        do 2000 i=1,n21
            z1(i) = z2(i) * z1(i)
            if(i.gt.1)then
                z1(nn + 2 - i ) = conjg(z1(i))
            endif
 2000   continue
!c-----
!c      ensure Nyquist frequency element is real
!c-----
        z1(n21) = cmplx(real(z1(n21)),0.0)
!c------
!c      compute inverse Fourier transform
!c-----
        call zfour(z1,nn,+1,dt,df)   
!c-----
!c      save the correlated function
!c-----
        do 20 i = 1,nn
            x(i) = real(z1(i))
   20   continue
        return
        end subroutine convolve
 
!--------------------------------------------------------------------------------

!================================================================================
! Find lag
   subroutine determine_lag(sig1,sig2,n1,n2,lag)

     integer(kind=si) :: part=0

     integer(kind=si), intent(in) :: n1, n2

     real(kind=cp), dimension(n1), intent(in) :: sig1
     real(kind=cp), dimension(n2), intent(in) :: sig2

     real(kind=cp), dimension(:), allocatable :: corr

     integer(kind=si), intent(out) :: lag

     integer(kind=si) :: it, ind
     real(kind=cp)    :: maxcorr

     !*** Take good parts define middle
     if (.not. allocated(corr)) allocate(corr(n1))

     if (modulo(n2,2) == 0) then
        ind = n2/2
     else
        ind = (n2-1)/2 +1
     endif

     !*** Cross-correlation
     call mycorrelation(sig1,sig2,n1,n2,corr,part)

     !*** Find delay
     maxcorr = maxval(abs(corr))
     do it=1,n1
        if (abs(corr(it)) == maxcorr) then
           lag = it-ind
           exit
!        else
!           stop 'there is a problem....'
        endif
     enddo

   end subroutine determine_lag
!--------------------------------------------------------------------------------

!================================================================================
! Compute stalta ratio and give first pick
   subroutine substalta(sig,n,nsta,nlta,crit,stalta,tpick)

     integer(kind=si), intent(in) :: n, nsta, nlta
     real(kind=cp), intent(in) :: crit
     real(kind=cp), dimension(n), intent(in) :: sig

     integer(kind=si), intent(out) :: tpick

     real(kind=cp), dimension(n) :: sta, lta
     real(kind=cp), dimension(n), intent(out) :: stalta
     real(kind=cp), dimension(n+2*nsta) :: tmpsta
     real(kind=cp), dimension(n+2*nlta) :: tmplta

     integer(kind=si) :: i

     !*** Compute the short time average (STA)
     tmpsta(1:nsta) = sig(1)
     tmpsta(nsta+1:nsta+n) = sig(:)
     tmpsta(nsta+n+1:n+2*nsta) = sig(n)
     sta = zero
     do i=1+nsta,n+nsta
        sta(i-nsta) = sum(tmpsta(i-nsta:i+nsta)**2)
     enddo
     sta = 0.5 * sta / nsta

     !*** Compute the long time average (LTA)
     tmplta(1:nlta) = sig(1)
     tmplta(nlta+1:nlta+n) = sig(:)
     tmplta(nlta+n+1:n+2*nlta) = sig(n)
     lta = zero
     do i=1+nlta,n+nlta
        lta(i-nlta) = sum(tmplta(i-nlta:i+nlta)**2)
     enddo
     lta = 0.5 * lta / nlta

     !*** Compute ratio and gives first pick
     stalta = sta / lta
     do i=1,n
        if (stalta(i) >= crit) then
           tpick = i
           exit
        else
           tpick = n
        endif
     enddo

   end subroutine substalta
!--------------------------------------------------------------------------------

!================================================================================
! Determine first arrival
   subroutine pick_first_arrival(sig1,n1,tpick,dt)

     integer(kind=si), intent(in) :: n1
     real(kind=cp)   , intent(in) :: dt
     real(kind=cp), dimension(n1), intent(in) :: sig1

     integer(kind=si), intent(out) :: tpick

     real(kind=cp), dimension(n1) :: stalta

     integer(kind=si) :: nsta, nlta
     real(kind=cp)    :: crit

     nsta = int(ceiling(5. / dt))   !* 5s pour sta
     nlta = int(ceiling(60. / dt))  !* 60s pour lta

     crit = 2.1                     !* a bit more than twice...

     call substalta(sig1,n1,nsta,nlta,crit,stalta,tpick)

   end subroutine pick_first_arrival
!--------------------------------------------------------------------------------

!================================================================================

! Trilinear interpolation

  subroutine trilin_interp(x,y,z,valx,valy,valz,nx,ny,nz,valin,valout)!,indx,indy,indz)

    !*** Size of grid
    integer(kind=si), intent(in) :: nx, ny, nz

    !*** Lookup 1D tables with Cartesian coordinates
    real(kind=cp), dimension(nx), intent(in) :: valx
    real(kind=cp), dimension(ny), intent(in) :: valy
    real(kind=cp), dimension(nz), intent(in) :: valz

    !*** Input grid of data
    real(kind=cp), dimension(nx,ny,nz), intent(in) :: valin

    !*** Coordinate for interpolated value
    real(kind=cp), intent(in) :: x, y, z

    !*** Temporary variables
    real(kind=cp) :: xix1, xix2, yiy1, yiy2, ziz1, ziz2
    real(kind=cp) :: x2x1, x1x2, y2y1, y1y2, z2z1, z1z2
    real(kind=cp) :: facx1y1z1, facx1y1z2, facx1y2z1, facx1y2z2, facx2y1z1, facx2y1z2, facx2y2z1, facx2y2z2

    !*** Output value to be interpolated
    real(kind=cp), intent(out) :: valout
    integer(kind=si) :: indx, indy, indz  !!! Previous guest
    integer(kind=si) :: kx, ky, kz, m

    !*** Lookup table
!    call locate_hunt(valx,nx,x,indx)
!    call locate_hunt(valy,ny,y,indy)
!    call locate_hunt(valz,nz,z,indz)
    call locate_bissection(valx,nx,x,indx)
    call locate_bissection(valy,ny,y,indy)
    call locate_bissection(valz,nz,z,indz)


    m=2
    kx = min(max(indx-(m-1)/2,1),nx+1-m)
    ky = min(max(indy-(m-1)/2,1),ny+1-m)
    kz = min(max(indz-(m-1)/2,1),nz+1-m)


    !*** x_i - x1
    xix1 = x - valx(kx)
    yiy1 = y - valy(ky)
    ziz1 = z - valz(kz)

    !*** x_i - x2
    xix2 = x - valx(kx+1)
    yiy2 = y - valy(ky+1)
    ziz2 = z - valz(kz+1)

    !*** x1 - x2
    x1x2 = 1./ (valx(kx) - valx(kx+1))
    y1y2 = 1./ (valy(ky) - valy(ky+1))
    z1z2 = 1./ (valz(kz) - valz(kz+1))

    !*** x2 - x1
    x2x1 = 1./(valx(kx+1) - valx(kx))
    y2y1 = 1./(valy(ky+1) - valy(ky))
    z2z1 = 1./(valz(kz+1) - valz(kz))

    !*** Factors
    facx1y1z1 = xix2*yiy2*ziz2 * x1x2*y1y2*z1z2 * valin(kx,ky,kz)
    facx1y1z2 = xix2*yiy2*ziz1 * x1x2*y1y2*z2z1 * valin(kx,ky,kz+1)
    facx1y2z1 = xix2*yiy1*ziz2 * x1x2*y2y1*z1z2 * valin(kx,ky+1,kz)
    facx1y2z2 = xix2*yiy1*ziz1 * x1x2*y2y1*z2z1 * valin(kx,ky+1,kz+1)
    facx2y1z1 = xix1*yiy2*ziz2 * x2x1*y1y2*z1z2 * valin(kx+1,ky,kz)
    facx2y1z2 = xix1*yiy2*ziz1 * x2x1*y1y2*z2z1 * valin(kx+1,ky,kz+1)
    facx2y2z1 = xix1*yiy1*ziz2 * x2x1*y2y1*z1z2 * valin(kx+1,ky+1,kz)
    facx2y2z2 = xix1*yiy1*ziz1 * x2x1*y2y1*z2z1 * valin(kx+1,ky+1,kz+1)

    !*** Final value
    valout = facx1y1z1 + facx1y1z2 + facx1y2z1 + facx1y2z2 + facx2y1z1 + facx2y1z2 + facx2y2z1 + facx2y2z2

  end subroutine trilin_interp
!--------------------------------------------------------------------------------



!================================================================================
! Locate with bisection
  subroutine locate_bissection(valx,n,x,ind)

    integer(kind=si), intent(in) :: n
    real(kind=cp),    intent(in) :: x

    real(kind=cp), dimension(n), intent(in) :: valx

    integer(kind=si), intent(out) :: ind

    integer(kind=si) :: jl, ju, jm

    jl = 0
    ju = n+1

    do while ( (ju-jl) > 1 )
       jm = (ju + jl) / 2
       if ( x >= valx(jm) ) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if ( x == valx(1) ) then
       ind = 1
    else if ( x == valx(n) ) then
       ind = n-1
    else
       ind = jl
    endif

  end subroutine locate_bissection
!--------------------------------------------------------------------------------


!================================================================================
! Locate with hunt
  subroutine locate_hunt(valx,n,x,jl)

    integer(kind=si), intent(in) :: n
    real(kind=cp),    intent(in) :: x

    real(kind=cp), dimension(n), intent(in) :: valx

    integer(kind=si), intent(inout) :: jl    ! inout because initial guess

    integer(kind=si) :: ju, jm,  inc

    inc = 1                      ! Set the hunting increment

    if (jl <= 0 .or. jl > n) then ! input guess not useful goto bissection
       jl = 0
       ju = n+1
    else
       if (x >= valx(jl)) then         ! Hunt up
          do while ( x >= valx(jl) )   ! Not done with hunting
             ju  = jl + inc
             if ( ju > n ) then        ! Done hunting... end of table
                ju = n+1
                exit
             else if (x < valx(ju)) then    ! Found bracket
                exit
             else                      ! Not done => double increment
                jl = ju
                inc = inc + inc
             endif
          enddo
       else                            ! Hunt down
          ju = jl
          do while ( x < valx(jl) )    ! not done hunting
             jl = jl - inc
             if (jl <= 0) then         ! Done hunting... end of table
                jl= 0
                exit
             else if (x >= valx(ju)) then   ! Found bracket
                exit
             else
                ju = jl
                inc = inc + inc
             endif
          enddo
       endif
    endif

    !** hunt is done
    ! Dichotomy
    do while ( (ju-jl) > 1 )
       jm = (ju + jl) / 2
       if ( x >= valx(jm) ) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if ( x == valx(1) ) then
       jl = 1
    else if ( x == valx(n) ) then
       jl = n-1
    endif

  end subroutine locate_hunt
!--------------------------------------------------------------------------------


!================================================================================
! Taper 3D
  subroutine taper_3D(ndom,taper,isr,sizetapx,sizetapy,sizetapz)

     integer(kind=si), dimension(3), intent(in) :: ndom

     integer(kind=si) :: i, j, k, isr

     real(kind=cp), dimension(:,:,:), allocatable, intent(inout) :: taper

     real(kind=cp), dimension(:), allocatable :: tapx, tapy, tapz
     real(kind=cp) :: alpha
     real(kind=cp),intent(in) :: sizetapx, sizetapy, sizetapz

     if (.not. allocated(tapx)) allocate(tapx(ndom(1)))
     if (.not. allocated(tapy)) allocate(tapy(ndom(2)))
     if (.not. allocated(tapz)) allocate(tapz(ndom(3)))

     alpha = sizetapx*2./ndom(1)
     tapx = tuckeywin(ndom(1),alpha)

     alpha = sizetapy*2./ndom(2)
     tapy = tuckeywin(ndom(2),alpha)

     alpha = sizetapz*2./ndom(3)
     tapz = tuckeywin(ndom(3),alpha)

     !*** Remove taper at the top
     if (isr == 2) then
        tapz(ndom(3)/2:ndom(3)) = 1.
        tapz(ndom(3))   = 0.1
        tapz(ndom(3)-1) = 0.2
        tapz(ndom(3)-2) = 0.3
        tapz(ndom(3)-3) = 0.4
        tapz(ndom(3)-4) = 0.5
        tapz(ndom(3)-5) = 0.6
        tapz(ndom(3)-6) = 0.7
        tapz(ndom(3)-7) = 0.8
        tapz(ndom(3)-8) = 0.9

     else if (isr == 1) then
                tapz(ndom(3)/2:ndom(3))=1.
     endif

     do k=1,ndom(3)
        do j=1,ndom(2)
           do i=1,ndom(1)
              taper(i,j,k) = tapx(i) * tapy(j) * tapz(k)
           enddo
        enddo
     enddo

     if (allocated(tapx)) deallocate(tapx)
     if (allocated(tapy)) deallocate(tapy)
     if (allocated(tapz)) deallocate(tapz)

  end subroutine taper_3D
!--------------------------------------------------------------------------------


!================================================================================
! Tukey tapering windows
!--------------------------------------------------
! N     : number of samples
! alpha : percentage of signal to taper
! tuk   : tapered window
  function tuckeywin(N,alpha) result(tuk)

    integer(kind=si), intent(in)   :: N
    real(kind=cp), intent(in)      :: alpha

    integer(kind=si) :: i
    real(kind=cp), parameter :: pipi=3.141592653589793
    real(kind=cp), dimension(N) :: tuk

    !*** Central part
    tuk(:) = 1.

    !*** Left part
    do i=0,int(0.5*alpha*(N-1))
       if (i+1 > 0 .and. i+1 <= N) then
          tuk(i+1) = 0.5*(1+cos(pipi*(2.*i/(alpha*(N-1.))-1.)))
       endif
    enddo

    !*** Right part
    do i=int((N-1)*(1-alpha/2.)),N-1
       if (i+1 > 0 .and. i+1 <= N) then
          tuk(i+1) = 0.5*(1+cos(pipi*(2.*i/(alpha*(N-1.))-(2./alpha)+1.)))
       endif
    enddo

  end function tuckeywin
!--------------------------------------------------------------------------------



!================================================================================
! Cardinal sine function
!--------------------------------------------------
! x = argument of sinc function
  real(kind=dp) function mysinc(x)

    real(kind=dp) :: x
    real(kind=dp), parameter :: pipi=3.141592653589793

    if (abs(x) >= 1e-13) then
       mysinc = sin(pipi*x)/(pipi*x)
    else
       mysinc = 1._dp
    endif

  end function mysinc
!--------------------------------------------------------------------------------


  subroutine interpolate_syn_alloc(syn,t1,dt1,npt1,t2,dt2,npt2)

   implicit none
   double precision, dimension(:),intent(inout) :: syn
   integer,intent(in) :: npt1,npt2
   double precision,intent(in) :: t1,dt1,t2,dt2

   double precision, dimension(:), allocatable :: syn1
   double precision :: time, tt
   integer i, ii

   ! initializes trace holding interpolated values
   ! syn1(1:npt2) = 0.
   allocate(syn1(npt2+npt1))
   syn1=0.

   ! loops over number of time steps in complete trace
   do i = 1, npt2

     ! sets time (in s) at this time step:
     ! t2 : start time of trace
     ! dt2: delta_t of a single time step
     time = t2 + (i-1) * dt2

     ! checks if time is within measurement window
     ! t1: start time of measurement window
     ! npt1: number of time steps in measurement window
     ! dt1: delta_t of a single time step in measurement window
     if (time > t1 .and. time < t1 + (npt1-1)*dt1) then

       ! sets index of time steps within this window: is 1 at the beginning of window
       ii = floor((time-t1)/dt1) + 1

       ! time increment within this single time step to match the exact value of time
       tt = time - ((ii-1)*dt1 + t1)

       ! interpolates value of trace for the exact time
       syn1(i) = (syn(ii+1)-syn(ii)) * tt/dt1 + syn(ii)
     endif
   enddo

   ! saves interpolated values to output trace
   syn(1:npt2) = syn1(1:npt2)
   ! LQY: zero out any thing beyond npts
   if (npt1 > npt2) syn(npt2+1:npt1)=0.
   deallocate(syn1)
 end subroutine interpolate_syn_alloc


!!$!================================================================================
!!$! Make lookup tables
!!$  subroutine create_lookup_tables(nx,ny,nz,xo,yo,zo,dx,dy,dz, &
!!$                                  xcrd, ycrd, zcrd)
!!$
!!$    integer(kind=si), intent(in) :: nx, ny, nz
!!$    real(kind=cp),    intent(in) :: xo, yo, zo
!!$    real(kind=cp),    intent(in) :: dx, dy, dz
!!$
!!$    integer(kind=si) :: xmid, ymid, zmid, i
!!$
!!$    real(kind=cp), dimension(nx), intent(out) :: xcrd
!!$    real(kind=cp), dimension(ny), intent(out) :: ycrd
!!$    real(kind=cp), dimension(nz), intent(out) :: zcrd
!!$
!!$    !*** Creates coordinates
!!$    !* For x
!!$    if (modulo(nx,2)==0) then
!!$       xmid = nx / 2
!!$    else
!!$       xmid = (nx / 2) + 1
!!$    endif
!!$
!!$    do i=1,nx
!!$       xcrd(i) = -dble(xmid - i) * dx
!!$    enddo
!!$
!!$    !* For y
!!$    if (modulo(ny,2)==0) then
!!$       ymid = ny / 2
!!$    else
!!$       ymid = (ny / 2) + 1
!!$    endif
!!$
!!$    do i=1,ny
!!$       ycrd(i) = -dble(ymid - i) * dy
!!$    enddo
!!$
!!$    !* For z
!!$    zmid = 1       ! Indice of zorigine = 0
!!$    do i=1,nz
!!$       zcrd(i) = -model_size(3) + dble(i-zmid) * inv_dz
!!$    enddo
!!$
!!$  end subroutine create_inversion_grid

end module interpolation_mod
