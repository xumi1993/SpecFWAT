module signal
  use config
  use ma_constants
  implicit none

contains

  subroutine bandpass_dp(x ,n, delta_t, f1, f2, order)
    ! a double-precision wrapper around sac xapiir()
    ! modified from FLEXWIN subroutines on 26-July-2009
    integer, intent(in) :: n, order
    real(kind=dp), intent(inout),  dimension(:) :: x
    real(kind=dp), intent(in) :: delta_t
    real(kind=cr), intent(in) :: f1,f2
    real(kind=cr), dimension(:), allocatable :: x_sngl

    allocate(x_sngl(n))

    x_sngl(1:n) = sngl(x(1:n))
    !  delta_t_sngl = sngl(delta_t)

    ! old version - uses old SacLib
    ! does band-pass filter
    !call xapiir(x_sngl,n,'BU',sngl(TRBDNDW),sngl(APARM),IORD,'BP',sngl(FSTART),sngl(FEND),delta_t_sngl,PASSES)

    ! new version, uses subroutines in libsac.a
    ! does band-pass filter
    ! BU - butterworth
    ! BP - bandpass
    ! LQY: Shouldn't  delta_t_sngl = sngl(delta_t) still be done? same for f1,f2?
    call xapiir(x_sngl,n,'BU',sngl(TRBDNDW),sngl(APARM),order,'BP',sngl(f1),sngl(f2),sngl(delta_t),PASSES)

    x(1:n) = dble(x_sngl(1:n))

    deallocate(x_sngl)

  end subroutine bandpass_dp

  subroutine interpolate_syn_dp(syn,t1,dt1,npt1,t2,dt2,npt2)
    real(kind=dp), dimension(:),intent(inout) :: syn
    integer, intent(in) :: npt1,npt2
    real(kind=dp),intent(in) :: t1,dt1,t2,dt2
    real(kind=dp), dimension(:), allocatable :: syn1
    real(kind=dp) :: time, tt
    integer i, ii

    ! initializes trace holding interpolated values
    allocate(syn1(npt2))
    syn1 = 0.

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

  end subroutine interpolate_syn_dp

  function interpolate_func_dp(syn, t1, dt1, npt1, t2, dt2, npt2) result(syn_interp)
    real(kind=dp), dimension(:), intent(in) :: syn
    integer, intent(in) :: npt1, npt2
    real(kind=dp), intent(in) :: t1, dt1, t2, dt2
    real(kind=dp), dimension(:), allocatable :: syn_interp
    real(kind=dp), dimension(:), allocatable :: syn1
    real(kind=dp) :: time, tt
    integer i, ii

    ! initializes trace holding interpolated values
    allocate(syn1(npt2))
    syn1 = 0.

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

    ! allocate and save interpolated values to output trace
    allocate(syn_interp(npt2))
    syn_interp(1:npt2) = syn1(1:npt2)

    ! deallocate temporary array
    deallocate(syn1)

  end function interpolate_func_dp

  subroutine dif1(data, dt)
    double precision, intent(in)                   :: dt
    double precision, dimension(:), intent(inout)  :: data
    double precision, dimension(:), allocatable    :: temp_signal
    integer :: i, nt

    nt = size(data)
    temp_signal=data(1:nt)

    do i=2,nt-1
      data(i) = - 0.5 * ( temp_signal(i+1) - temp_signal(i-1) ) / dt
    enddo

    data(1)=0.
    data(nt)=0.
  end subroutine dif1

  subroutine myconvolution_dp(sig1,sig2,conv,part)

    integer, intent(in) :: part
    integer :: n1, n2
    real(kind=dp), dimension(:), intent(in) :: sig1
    real(kind=dp), dimension(:), intent(in) :: sig2
    real(kind=dp), dimension(:), allocatable, intent(out) ::conv
    real(kind=dp), dimension(:), allocatable :: convtmp !, intent(out) :: conv
    integer :: i1, i2, ind
    
    n1 = size(sig1)
    n2 = size(sig2)
    allocate(convtmp(n1+n2-1))
    convtmp = 0.0_dp

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
          ind = ceiling(real(n2/2,kind=cr))
      endif
      allocate(conv(n2))
      conv(1:n2) = convtmp(ind:ind+n2-1)

    else if (part == 1) then ! full convolution

      allocate(conv(n1+n2-1))
      conv(:) = convtmp(:)

    else if (part == 2) then !(middle irregular)
      allocate(conv(n2))
      conv(1:n2) = convtmp(n2:n2+n2-1)
    endif

  end subroutine myconvolution_dp

  subroutine mycorrelation_dp(sig1,sig2,corr,part)

    integer, intent(in) :: part
    integer :: n1, n2
    real(kind=dp), dimension(:), intent(in) :: sig1
    real(kind=dp), dimension(:), intent(in) :: sig2
    real(kind=dp), dimension(:), allocatable, intent(out) :: corr
    real(kind=dp), dimension(:), allocatable :: flipsig2
    integer :: i

    n1 = size(sig1)
    n2 = size(sig2)
    allocate(flipsig2(n2))
    flipsig2 = 0._dp

    !*** Choose size of corr
    if (part == 0) then !(middle)
      allocate(corr(n2))
    else if (part == 1) then !(full)
      allocate(corr(n1+n2-1))
    endif

    !*** Flip second signal
    do i=1,n2
      flipsig2(i) = sig2(n2-i+1)
    enddo

    !*** Use regular convolution
    call myconvolution_dp(sig1,flipsig2,corr,part)

  end subroutine mycorrelation_dp

  ! Convolution routine
  subroutine myconvolution(sig1,sig2,conv,part)

    integer, intent(in) :: part

    real(kind=cr), dimension(:), intent(in) :: sig1
    real(kind=cr), dimension(:), intent(in) :: sig2

!     real(kind=cp), dimension(n1), intent(out) ::conv
    real(kind=cr), dimension(:), allocatable, intent(inout) ::conv

    real(kind=cr), dimension(:), allocatable :: convtmp !, intent(out) :: conv

    integer :: i1, i2, ind, n1, n2

    n1 = size(sig1)
    n2 = size(sig2)
    allocate(convtmp(n1+n2-1))

    !*** Put to zero
    convtmp = 0.0_cr

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
          ind = ceiling(real(n2/2,kind=cr))
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


  subroutine mycorrelation(sig1,sig2,corr,part)

    integer, intent(in) :: part

    real(kind=cr), dimension(:), intent(in) :: sig1
    real(kind=cr), dimension(:), intent(in) :: sig2

!     real(kind=cp), dimension(n1), intent(out) :: corr
    real(kind=cr), dimension(:), allocatable, intent(inout) :: corr

    real(kind=cr), dimension(:), allocatable :: flipsig2
    integer :: i, n1, n2

    n1 = size(sig1)
    n2 = size(sig2)
    allocate(flipsig2(n2))
    flipsig2 = 0._cr

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
    call myconvolution(sig1,flipsig2,corr,part)

  end subroutine mycorrelation

  !================================================================================
! Iterative time domain deconvolution
  subroutine time_deconv(dobs,dcal,dt,nt,nit,src_sum)

    integer, intent(in) :: nt, nit
    integer             :: i, ii
    integer, parameter  :: part=0

    real(kind=cr), intent(in) :: dt
    real(kind=cr)             :: amp_corr

    real(kind=cr), dimension(nt), intent(in)    :: dobs, dcal
    real(kind=cr), dimension(:), allocatable, intent(inout) :: src_sum

    real(kind=cr), dimension(:), allocatable :: dobs2, autocorr, crosscorr, new_obs, src_one


    !* 0. Alloctae
    if (.not. allocated(autocorr))  allocate(autocorr(nt))
    if (.not. allocated(crosscorr)) allocate(crosscorr(nt))
    if (.not. allocated(src_sum))   allocate(src_sum(nt))
    if (.not. allocated(src_one))   allocate(src_one(nt))
    if (.not. allocated(dobs2))     allocate(dobs2(nt))
    if (.not. allocated(new_obs))   allocate(new_obs(nt))
    src_sum = 0._cr
    src_one = 0._cr

    !* 1. Estimate auto-corr of computed data
    call mycorrelation(dcal,dcal,autocorr,part)
    autocorr = autocorr * dt
    amp_corr = 1._cr / maxval(abs(autocorr))
   !  write(*,*) 'amp_corr',amp_corr

    !* 2. Time iteration to estimate sources
    dobs2 = dobs
    do i = 1, nit

       !* Compute correlation
       call mycorrelation(dobs2,dcal,crosscorr,part)
       crosscorr = crosscorr * dt

       !* Find maximum of correlation
       ii = maxloc(abs(crosscorr),dim=1)

       !* Put local contibution to src_one and src_sum
       src_one     = 0._cr
       src_one(ii) = crosscorr(ii) * amp_corr
       src_sum     = src_sum + src_one

       !* Convolve
       call myconvolution(src_one,dcal,new_obs,part)
       new_obs = new_obs * dt
       dobs2 = dobs2 - new_obs

    enddo

  end subroutine time_deconv

  subroutine detrend_legacy(x)
    implicit none
    real(kind=dp), dimension(:) :: x
    real(kind=dp) :: ds1,ds2,dan,davei,davex,dslope,dai
    integer :: i, an, n

    n = size(x)
    an = n
    dan=n
    ds1=0
    ds2=0

    do i=1,n
      ds1 = ds1+ x(i)
      ds2 = ds2 + ds1
    enddo
    davei = 0.5 * (dan+1.0)
    davex = ds1/dan
    dslope = -12.0*(ds2-davei*ds1)/(dan*(dan*dan-1.0))
    do i=1,n
      dai = i-1
      x(i) = x(i)- davex - dslope*(dai-davei)
    enddo

  end subroutine detrend_legacy

  subroutine detrend(x)
    implicit none
    real(kind=dp), dimension(:), intent(inout) :: x
    real(kind=dp) :: x1, x2
    integer :: i, n

    n = size(x)

    x1 = x(1)
    x2 = x(n)

    do i=1,n
      x(i) = x(i) - (x1 + (x2 - x1) * real(i-1, kind=dp) / real(n-1, kind=dp))
    enddo

  end subroutine detrend

  subroutine demean(x)
    implicit none
    double precision, dimension(:), intent(inout) :: x
    double precision :: sum
    integer :: i, n

    n = size(x)
    x = x - sum(x)/n

  end subroutine demean

end module signal