module fftpack
  implicit none
  integer, private, parameter :: dp = kind(1.0d0), sp = kind(1.0)
  real(kind=sp), parameter, private :: PI = 3.14159265358979323846_sp

  type, public :: fft_cls
    contains
      procedure :: fft, fft_dp
      procedure :: ifft
      procedure :: ifft_complex
      procedure :: fftfreq
      procedure :: hilbert
  end type fft_cls

contains

  function fft_dp(this, x, nft) result(fpx)
    class(fft_cls) :: this
    integer, intent(in) :: nft
    real(kind=dp), dimension(nft) :: x
    complex(kind=dp), dimension(nft) :: fpx
    complex(kind=sp), dimension(nft) :: fpr
    integer :: n_input

    n_input = size(x)
    if (n_input > nft) then
      write(*,*) 'Error: fft_dp input size greater than nft'
      error stop
    end if

    fpr = 0.0
    fpr(1:n_input) = cmplx(real(x), 0.0)
    call fft_raw(nft, fpr, -1)
    fpx = cmplx(fpr, kind=dp)
  end function fft_dp

  function fft(this, x, nft) result(fpx)
    class(fft_cls) :: this
    integer, intent(in) :: nft
    real(kind=dp), dimension(nft) :: x
    complex(kind=sp), dimension(nft) :: fpx
    integer :: n_input
    
    n_input = size(x)
    if (n_input > nft) then
      write(*,*) 'Error: fft input size greater than nft'
      error stop
    end if
    
    fpx = 0.0
    fpx(1:n_input) = cmplx(real(x), 0.0)
    call fft_raw(nft, fpx, -1)
  end function fft

  function ifft(this, fpx, nft) result(x)
    class(fft_cls) :: this
    integer :: nft
    double precision, dimension(nft) :: x
    complex, dimension(nft) :: fpx, fp_norm
    fp_norm = fpx/nft
    call fft_raw(nft, fp_norm, +1)
    x = dble(real(fp_norm))
  end function ifft

  function ifft_complex(this, fpx, nft) result(x)
    class(fft_cls) :: this
    integer :: nft
    complex(kind=dp), dimension(nft) :: x
    complex, dimension(nft) :: fpx, fp_norm
    fp_norm = fpx/nft
    call fft_raw(nft, fp_norm, +1)
    x = cmplx(fp_norm, kind=dp)
  end function ifft_complex

  subroutine fft_raw(n, x, ind)
    integer :: n, ind, j, i, kmax, istep,m, k
    complex :: temp, theta
    complex, intent(inout) :: x(n)

    j = 1
    do i = 1, n
      if (i < j) then
        temp = x(j)
        x(j) = x(i)
        x(i) = temp
      endif
      m = n/2
      do while (.true.)
        if (j > m) then
          j = j - m
          m = m / 2
          if (m < 2) exit
        else
          exit
        endif
      enddo
      j = j + m
    enddo
    kmax = 1
    do while (kmax < n)
      istep = kmax * 2
      do k = 1, kmax
        theta = cmplx(0.0, PI * ind * (k - 1) / kmax)
        do i = k, n, istep
          j = i + kmax
          temp = x(j) * cexp(theta)
          x(j) = x(i) - temp
          x(i) = x(i) + temp
        enddo
      enddo
      kmax=istep
    enddo
  end subroutine fft_raw

  function fftfreq(this, n, d) result(freqs)
    class(fft_cls) :: this
    integer, intent(in) :: n         
    double precision, intent(in) :: d        
    double precision, allocatable :: freqs(:) 

    integer :: i, N_half

    allocate(freqs(n))
    N_half = n / 2

    do i = 1, N_half
      freqs(i) = (i - 1) / (n * d)
    end do

    do i = N_half + 1, n
      freqs(i) = - (n - i + 1) / (n * d)
    end do

  end function fftfreq

  function hilbert(this, x, n_opt) result(xa)
    ! FFT-based computation of the analytic signal
    ! The analytic signal is calculated by filtering out the negative frequencies and
    ! doubling the amplitudes of the positive frequencies in the FFT domain
    
    class(fft_cls) :: this
    real(kind=dp), dimension(:), intent(in) :: x     ! Input real signal
    integer, intent(in), optional :: n_opt           ! Number of Fourier components
    complex(kind=dp), dimension(:), allocatable :: xa ! Analytic signal output
    
    ! Local variables
    integer :: n, n_input, nft
    real(kind=dp), dimension(:), allocatable :: x_padded
    complex, dimension(:), allocatable :: xf, h_weights
    
    n_input = size(x)
    
    ! Set the FFT size
    if (present(n_opt)) then
      if (n_opt < n_input) then
        n = n_input
      else
        n = n_opt
      end if
    else
      n = n_input
    end if
    
    nft = 2**exponent(real(n))

    ! Prepare input data with zero padding if necessary
    allocate(x_padded(nft))
    x_padded = 0.0_dp
    x_padded(1:n_input) = x(1:n_input)
    
    ! Compute FFT
    ! allocate(xf(n))
    xf = this%fft(x_padded, nft)
    
    ! Create Hilbert filter weights
    allocate(h_weights(nft))
    h_weights = (0.0, 0.0)
    
    if (mod(n, 2) == 0) then
      ! Even length
      h_weights(1) = (1.0, 0.0)          ! DC component
      h_weights(nft/2 + 1) = (1.0, 0.0)    ! Nyquist frequency
      h_weights(2: nft/2) = (2.0, 0.0)      ! Positive frequencies
      ! Negative frequencies remain zero
    else
      ! Odd length
      h_weights(1) = (1.0, 0.0)          ! DC component
      h_weights(2: (nft + 1)/2) = (2.0, 0.0) ! Positive frequencies
      ! Negative frequencies remain zero
    end if
    
    ! Apply Hilbert filter
    xf = xf * h_weights
    
    ! Compute inverse FFT to get analytic signal
    ! xa = this%ifft_complex(xf, n)
    call fft_raw(nft, xf, +1)
    allocate(xa(n))
    xa = cmplx(xf(1:n)/real(nft), kind=dp)
    deallocate(x_padded, xf, h_weights)
    
  end function hilbert
end module fftpack