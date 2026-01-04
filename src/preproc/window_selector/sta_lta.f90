module sta_lta_mod
  use signal
  implicit none
  contains

  function STA_LTA(env_synt_lp, dt, min_period)
    real(kind=dp), dimension(:), intent(in) :: env_synt_lp
    real(kind=dp), intent(in) :: dt, min_period
    real(kind=dp) :: TOL, Cs, Cl, sta, lta, noise
    integer :: i, n_extend
    real(kind=dp), dimension(:), allocatable :: extended_syn
    real(kind=dp), dimension(:), allocatable :: STA_LTA
    integer :: npts

    ! the sta_lta value peaks at the begining of the trace.
    ! there is no phase incoming at that time, but there is some small trend to positive values.
    ! tweak only sets the sta/lta value when a threshold value in the envelope is reached
    logical,parameter:: TWEAK = .false.
    integer :: thres_noise
    real(kind=dp) :: lta_max, lta_org

  ! set the Cs Cl for STA/LTA ratio using the Bai & Kennett (2001) expressions

    Cs = 10.0_dp**(-dt/min_period)
    Cl = 10.0_dp**(-dt/(12.0_dp*min_period))
    TOL = 1.0e-9_dp

  ! set pre-extension for synthetic data and allocate extended_syn
  !  n_extend=5*12*min_period/dt
    npts = size(env_synt_lp)
    n_extend=1000
    allocate(extended_syn(npts+n_extend))

  ! set noise level
  !  n_sample=12*min_period/dt
  !  noise=sum(env_synt_lp(1:n_sample))/n_sample
    noise=maxval(env_synt_lp)/1e5

    ! determines the onset of the envelope (above some threshold value)
    if( TWEAK ) then
      do i=1,npts
        if( env_synt_lp(i) >= noise*1000. ) then
          thres_noise = i
          exit
        endif
      enddo
    endif

  ! copy the original synthetic into the extended array, right justified
  !  call random_number(extended_syn)
  !  extended_syn=noise*extended_syn
    extended_syn=noise
    extended_syn(n_extend+1:n_extend+npts)=env_synt_lp(1:npts)+noise

    ! if (DEBUG) write(*,*) 'DEBUG : Cs, Cl = ', Cs, Cl
    ! if (DEBUG) write(*,*) 'Number of points used to pre-extend synthetics ', n_extend

  ! initialise sta/lta variables
    allocate(STA_LTA(npts))
    STA_LTA(:)=0.0
    sta=0.0 ; lta = 0.0

  ! warm up the sta and lta
    do i = 1, n_extend
      sta=(Cs*sta+extended_syn(i))
      lta=(Cl*lta+extended_syn(i))
    enddo

    ! determine long term average maximum value
    if( TWEAK ) then
      lta_org = lta
      lta_max = 0.0
      do i = 1, npts
        lta=(Cl*lta+extended_syn(n_extend+i))
        if( lta > lta_max ) lta_max = lta
      enddo
      ! if (DEBUG) write(*,*) 'DEBUG : lta_max = ', lta_max
      lta = lta_org
    endif

  ! calculate sta/lta for the envelope
    do i = 1, npts
      sta=(Cs*sta+extended_syn(n_extend+i))
      lta=(Cl*lta+extended_syn(n_extend+i))
      if (lta.gt.TOL) then
        if( TWEAK ) then
          ! additional envelope criteria
          if( lta .gt. TOL*lta_max ) then
            if( i >= thres_noise) then
              STA_LTA(i)=sta/lta
            endif
          endif
        else
          STA_LTA(i)=sta/lta
        endif
      else
         STA_LTA(i) = noise ! Match Pyflex: sta[lta < TOL] = noise
      endif

    enddo

    deallocate(extended_syn)

  end function STA_LTA

end module sta_lta_mod