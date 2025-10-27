module travel_times_mod
  use config
  implicit none
  integer, parameter :: MAX_PHASES=60, phase_name_len=8
  
contains
  subroutine ttimes(dist_deg,evdp,nphases,names,times)

  real(kind=dp), intent(in) :: dist_deg, evdp
  integer, intent(out) :: nphases 
  character(len=phase_name_len), dimension(:), allocatable, intent(out) :: names
  real(kind=dp), dimension(:), allocatable, intent(out) :: times

  ! legacy variables needed to use the libtau routines
  logical :: prnt(3)
  character(len=phase_name_len) :: phlst(10)
  real :: usrc(2)
  real, dimension(MAX_PHASES) :: dtdd,dtdh,dddp
  real, dimension(MAX_PHASES) :: times_sngl
  character(len=phase_name_len), dimension(MAX_PHASES) :: phnames
  character(len=MAX_STRING_LEN) :: modnam

  ! ask for all phases
  phlst(1)="P"
  prnt(1)=.false.
  prnt(2)=.false.
  prnt(3)=.false.
  modnam=trim(TAUP_DATA_DIR)//'/'//'iasp91'
  ! modnam='iasp91'
  call tabin(1,modnam)
  call brnset(1,phlst,prnt)
  call depset(sngl(evdp),usrc)

  call trtm(sngl(dist_deg),MAX_PHASES,nphases,times_sngl,dtdd,dtdh,dddp,phnames)
  allocate(names(nphases))
  allocate(times(nphases))
  times(1:nphases) = dble(times_sngl(1:nphases))
  names(1:nphases) = phnames(1:nphases)
  end subroutine ttimes

end module travel_times_mod