module fwat_utils
  use constants
  use ma_constants, only: NCHI, TRBDNDW, APARM, IORD, PASSES
  use my_mpi
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: RPRE = SIZE_REAL;
  integer, parameter :: DPRE = SIZE_DOUBLE;
  integer, parameter :: IPRE = SIZE_INTEGER;

  interface write_timestamp_log
    module procedure log_str, log_real, log_int, log_double
  end interface

  type, public :: chi_table
    integer :: n_rows
    integer :: n_cols = 12+NCHI
    character(len=20), dimension(:), allocatable :: evtid,sta,net,chan
    integer, dimension(:), allocatable :: idx_rec, meas
    double precision, dimension(:,:), allocatable :: window_chi
    double precision, dimension(:), allocatable :: tr_chi, am_chi, T_pmax_dat, T_pmax_syn, tstart,tend

    contains
    procedure, public :: read_misfits
    procedure, public :: mean
    procedure, public :: std
    procedure, public :: init_table
    procedure, public :: distory
    procedure, public :: get_column
    procedure, public :: apply_weight
  end type

contains
  subroutine nowtime(timestamp)
    character(len=MAX_STRING_LEN), intent(out) :: timestamp
    integer, dimension(8) :: time_values

    call date_and_time(VALUES=time_values)
    write(timestamp,'(i4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2)') time_values(1),time_values(2), &
          time_values(3),time_values(5),time_values(6),time_values(7)
    
  end subroutine

  subroutine log_str(fid, value)
    character(len=*) :: value
    character(len=MAX_STRING_LEN) :: timestamp
    integer :: fid

    call nowtime(timestamp)
    write(fid, '("[",a,"] ",a)') trim(timestamp), trim(value) 
  end subroutine

  subroutine log_real(fid, value)
    real :: value
    character(len=MAX_STRING_LEN) :: timestamp
    integer :: fid
     
    call nowtime(timestamp)
    write(fid, '("[",a,"] ",F6.4)') trim(timestamp), value
  end subroutine

  subroutine log_int(fid, value)
    integer :: value
    character(len=MAX_STRING_LEN) :: timestamp
    integer :: fid
     
    call nowtime(timestamp)
    write(fid, '("[",a,"] ",I6)') trim(timestamp), value
  end subroutine

  subroutine log_double(value)
    double precision :: value
    character(len=MAX_STRING_LEN) :: timestamp
    integer :: fid

    call nowtime(timestamp)
    write(fid, '("[",a,"] ",F12.8)') trim(timestamp), value
  end subroutine

  subroutine read_misfits(this, fname)
    character(len=MAX_STRING_LEN) :: fname
    class(chi_table) :: this
    integer :: ierr, i

    open(unit=300,file=Trim(fname),status='old')
    this%n_rows = 0
    do 
      read(300, *, iostat=ierr)
      if (ierr /= 0) exit
      this%n_rows = this%n_rows+1
    enddo
    call this%init_table()
    rewind(300)
    do i = 1, this%n_rows
      read(300,'(a20,a8,a4,a5,i4,i4,2e14.6,20e14.6,2e14.6,2f14.6)') &
        this%evtid(i),this%sta(i),this%net(i),this%chan(i),this%idx_rec(i), &
        this%meas(i),this%tstart(i),this%tend(i),this%window_chi(i, :),&
        this%tr_chi(i),this%am_chi(i),this%T_pmax_dat(i),this%T_pmax_syn(i)
    enddo
    close(300)

  end subroutine

  subroutine init_table(this)
    class(chi_table) :: this

    allocate( &
      this%evtid(this%n_rows),this%sta(this%n_rows),this%net(this%n_rows),this%chan(this%n_rows), &
      this%idx_rec(this%n_rows), this%meas(this%n_rows), this%tstart(this%n_rows),this%tend(this%n_rows),&
      this%window_chi(this%n_rows,NCHI), &
      this%tr_chi(this%n_rows), this%am_chi(this%n_rows), &
      this%T_pmax_dat(this%n_rows), this%T_pmax_syn(this%n_rows) &
    )
  end subroutine

  subroutine distory(this)
   class(chi_table) :: this

   deallocate(this%evtid,this%sta,this%net,this%chan, &
      this%idx_rec, this%meas, this%tstart,this%tend, this%window_chi, &
      this%tr_chi, this%am_chi, &
      this%T_pmax_dat, this%T_pmax_syn &
    )
  end subroutine

  function get_column(this, col) result(array)
    class(chi_table) :: this
    integer :: col, i
    double precision, dimension(this%n_rows) :: array

    if (col <= 8 .or. col > this%n_cols) then
      stop 'mean value will be calculated when column > 8'
    elseif (col <= 28) then
      i = col - 8
      array(:) = this%window_chi(:, i)
    elseif (col == 29) then
      array(:) = this%tr_chi(:)
    elseif (col == 30) then
      array(:) = this%am_chi(:)
    elseif (col == 31) then
      array(:) = this%T_pmax_dat(:)
    elseif (col == 32) then
      array(:) = this%T_pmax_syn(:)
    endif

  end function

  function mean(this, col) result(mean_value)
    class(chi_table) :: this
    integer :: col
    double precision :: mean_value
    double precision, dimension(this%n_rows) :: array

    array = this%get_column(col)
    mean_value = sum(array)/this%n_rows
  end function

  function std(this, col) result(std_value)
  class(chi_table) :: this
  integer :: i, col
  double precision :: std_value, mean_value, variance
  double precision, dimension(this%n_rows) :: array

  array = this%get_column(col)
  mean_value = sum(array)/this%n_rows
  variance = 0.0
  do i = 1, this%n_rows
    variance = variance + (array(i)-mean_value)**2
  enddo
  std_value = sqrt(variance)
  end function

  function apply_weight(this, col, nevt, evtnames, weight) result(array)
    class(chi_table) :: this
    integer :: col, i, j, nevt
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: evtnames
    real(kind=CUSTOM_REAL), dimension(:), allocatable  :: weight
    double precision, dimension(this%n_rows) :: array

    array = this%get_column(col)

    do i = 1, this%n_rows
      do j = 1, nevt
        if (trim(this%evtid(i)) == trim(evtnames(j))) then
          array(i) = array(i) * weight(j) 
        endif
      enddo
    enddo

  end function

  subroutine read_evt_weight(evtset, nevt, evtnames, weight)

    character(len=MAX_STRING_LEN) :: evtset, evtset_file,line_junk
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: evtnames
    real(kind=CUSTOM_REAL), dimension(:), allocatable  :: weight
    integer :: nevt, ier,ievt
    real :: junk_la, junk_lo, junk_dp, junk_br

    evtset_file='src_rec/sources_'//trim(evtset)//'.dat'
    open(unit=400,file=trim(evtset_file),status='old',iostat=ier)
    if (ier /=0) then
      print *,'Error could not open source subset file: ',trim(evtset_file) 
      stop 'Error opening source subset file'
    endif
    nevt=0
    do 
      read(400,*,iostat=ier) line_junk  ! read event name
      !print*,trim(evtnm)
      if (ier /=0) exit
      nevt=nevt+1
    enddo
    rewind(400)
    if(allocated(evtnames)) deallocate(evtnames)
    if(allocated(weight)) deallocate(weight)
    allocate(evtnames(nevt), weight(nevt))
    do ievt=1,nevt
      read(400, '(a)', iostat=ier) line_junk
      if (ier /=0) exit
      read(line_junk,*,iostat=ier) evtnames(ievt),junk_la, junk_lo,&
                                   junk_dp, junk_br, weight(ievt)
      if (ier /= 0) weight(ievt) = 1.0
    enddo
    close(400)
  end subroutine

  subroutine get_band_name(SHORT_P, LONG_P, bandname)
    real :: SHORT_P, LONG_P, fl, fh
    character(len=MAX_STRING_LEN) :: bandname, bandstr1, bandstr2

    fh = 1./SHORT_P
    fl = 1./LONG_P
    if (SHORT_P < 1) then
      write(bandstr1, '(a1,i3.3)') 'F',int(fh)
    else
      write(bandstr1, '(a1,i3.3)') 'T',int(SHORT_P)
    endif
    if (LONG_P < 1) then
      write(bandstr2, '(a1,i3.3)') 'F',int(fl)
    else
      write(bandstr2, '(a1,i3.3)') 'T',int(LONG_P)
    endif
    bandname = trim(bandstr1)//'_'//trim(bandstr2)

  end subroutine

  ! Rotation of components
  subroutine rotate_ZNE_to_ZRT(vz,vn,ve,vz2,vr,vt,nt,bazi)
    use specfem_par, only: CUSTOM_REAL
  
    integer,                  intent(in) :: nt
    real(kind=CUSTOM_REAL),   intent(in) :: bazi

    real(kind=CUSTOM_REAL), dimension(nt),  intent(in) :: vz,  vn, ve
    real(kind=CUSTOM_REAL), dimension(nt), intent(out) :: vz2, vr, vt

    real(kind=CUSTOM_REAL) :: baz

    integer :: it

    baz = deg2rad * bazi

    do it = 1, nt
        vr(it) = -ve(it) * sin(baz) - vn(it) * cos(baz)
        vt(it) = -ve(it) * cos(baz) + vn(it) * sin(baz)
        vz2(it) = vz(it)
    enddo
  
  end subroutine rotate_ZNE_to_ZRT
  
  subroutine rotate_ZRT_to_ZNE(vz2,vr,vt,vz,vn,ve,nt,bazi)
    use specfem_par, only: CUSTOM_REAL
  
    integer,                  intent(in) :: nt
    real(kind=CUSTOM_REAL),   intent(in) :: bazi

    real(kind=CUSTOM_REAL), dimension(nt),  intent(in) :: vz2, vr, vt
    real(kind=CUSTOM_REAL), dimension(nt), intent(out) :: vz,  vn, ve

    real(kind=CUSTOM_REAL) :: baz
    integer :: it

    baz = deg2rad * bazi

    do it = 1, nt
        ve(it) = -vr(it) * sin(baz) - vt(it) * cos(baz)
        vn(it) = -vr(it) * cos(baz) + vt(it) * sin(baz)
        vz(it) = vz2(it)
    enddo
  
  end subroutine rotate_ZRT_to_ZNE

  function find_maxima(x) result(idx)
    implicit none

    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x
    integer, dimension(:), allocatable :: idx
    integer :: i, n, count
    
    logical, allocatable :: rising(:), falling(:)
    integer, allocatable :: max_idx(:)

    n = size(x)
    if (n < 3) then
      allocate(idx(0))
      return
    end if

    allocate(rising(n-1), falling(n-1))

    rising = (x(2:) - x(1:n-1)) > 0.0
    falling = (x(2:) - x(1:n-1)) < 0.0

    allocate(max_idx(n-2))
    count = 0
    do i = 1, n-2
      if (rising(i) .and. falling(i+1)) then
        count = count + 1
        max_idx(count) = i + 1
      end if
    end do

    allocate(idx(count))
    idx = max_idx(1:count)

  end function

  subroutine bandpass_iir(x,n,delta_t,f1,f2)
    ! a double-precision wrapper around sac xapiir()
    ! modified from FLEXWIN subroutines on 26-July-2009

    implicit none
    integer, intent(in) :: n
    double precision, intent(inout),  dimension(*) :: x
    double precision, intent(in) :: delta_t,f1,f2
    real, dimension(:), allocatable :: x_sngl

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
    call xapiir(x_sngl,n,'BU',sngl(TRBDNDW),sngl(APARM),IORD,'BP',sngl(f1),sngl(f2),sngl(delta_t),PASSES)

    x(1:n) = dble(x_sngl(1:n))

    deallocate(x_sngl)


  end subroutine bandpass_iir
end module