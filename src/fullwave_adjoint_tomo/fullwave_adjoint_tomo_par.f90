!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
!

module fullwave_adjoint_tomo_par

  !! IMPORT VARIABLES ------------------------------------------------------------------------------------------------
  use specfem_par, only: CUSTOM_REAL, MAX_STRING_LEN, MAX_LENGTH_STATION_NAME, MAX_LENGTH_NETWORK_NAME
  use my_mpi

  !-------------------------------------------------------------------------------------------------------------------

  implicit none

  !! ------------------------------ compialtion config parameters -----------------------------------------------------------------

  !! log file for inversion
  integer,                       public, parameter  :: OUT_FWAT_LOG=888, FIIN=889, FIOUT=890
  logical,                       public, parameter  :: DEBUG_MODE=.false.
  integer,                        public, parameter :: NUM_INV_TYPE = 3
  character(len= MAX_STRING_LEN), public, parameter :: FWAT_PAR_FILE = 'DATA/FWAT.PAR'
  logical,                       public             :: VERBOSE_MODE=.false.
  real(kind=CUSTOM_REAL), public, dimension(:), allocatable  :: src_weight
  integer, public                                                   :: NRCOMP, NUM_FILTER
  character(len= MAX_STRING_LEN), public, dimension(:), allocatable :: RCOMPS
  character(len= MAX_STRING_LEN), public                            :: CH_CODE, dat_coord
  real(kind=CUSTOM_REAL), public, dimension(:),       allocatable     :: SHORT_P
  real(kind=CUSTOM_REAL), public, dimension(:),       allocatable     :: LONG_P
  character(len= MAX_STRING_LEN)                                   :: source_fname
  character(len= MAX_STRING_LEN)                                   :: station_fname
  character(len=MAX_STRING_LEN), public                            :: PRECOND_TYPE
  logical,                       public                             :: is_read_database = .true.

!################################################# WORKFLOW ######################################################################
  type fwat_acqui
    integer               :: nevents ! # of events
    character(len= MAX_STRING_LEN), dimension(:),       allocatable  :: station_file
    character(len= MAX_STRING_LEN), dimension(:),       allocatable  :: src_solution_file
    character(len= MAX_STRING_LEN), dimension(:),       allocatable  :: evtid_names
    character(len= MAX_STRING_LEN), dimension(:),       allocatable  :: out_fwd_path
    character(len= MAX_STRING_LEN), dimension(:),       allocatable  :: in_dat_path
    real(kind=CUSTOM_REAL), dimension(:), allocatable                :: evla,evlo,evdp

    contains
    procedure :: free => fwat_acqui_free
    procedure :: malloc => fwat_acqui_malloc
    procedure :: read_source_set => fwat_acqui_read_source_set
  end type fwat_acqui


!###################### FWAT PAR for Noise tomo #####################
  type fwat_noise_parameters
    integer, public                                                   :: NSCOMP
    character(len= MAX_STRING_LEN), public, dimension(:), allocatable :: SCOMPS
    integer, public                                                   :: NRCOMP
    character(len= MAX_STRING_LEN), public, dimension(:), allocatable :: RCOMPS
    integer, public                                                   :: NUM_FILTER, NSTEP
    character(len= MAX_STRING_LEN), public                            :: CH_CODE, dat_coord
    character(len= MAX_STRING_LEN), public                            :: PRECOND_TYPE = 'default'
    real(kind=CUSTOM_REAL), public                                    :: DT
    real(kind=CUSTOM_REAL), public, dimension(:),       allocatable   :: SHORT_P
    real(kind=CUSTOM_REAL), public, dimension(:),       allocatable   :: LONG_P
    real(kind=CUSTOM_REAL), public, dimension(:), allocatable         :: GROUPVEL_MIN !BinHe added advised by Kai
    real(kind=CUSTOM_REAL), public, dimension(:), allocatable         :: GROUPVEL_MAX  !BinHe added advised by Kai
    logical, public                                                   :: USE_NEAR_OFFSET !BinHe added
    ! logical,                       public :: SHOW_DETAILS !BinHe added
    ! logical,                       public :: UPDATE_RHO !BinHe added
    logical, public                                                   :: ADJ_SRC_NORM
    logical, public                                                   :: SUPPRESS_EGF
    integer, public                                                   :: IMEAS, ITAPER

    contains
    procedure :: read_noise_par_file => fwat_noise_read_par_file
  end type fwat_noise_parameters

!###################### FWAT PAR for Tele FWI #######################
  type fwat_tele_parameters
    integer, public                                                   :: NRCOMP
    character(len= MAX_STRING_LEN), public, dimension(:), allocatable :: RCOMPS
    integer, public                                                   :: NUM_FILTER, NSTEP
    real(kind=CUSTOM_REAL), public                                    :: DT
    character(len= MAX_STRING_LEN), public                            :: CH_CODE, dat_coord
    character(len= MAX_STRING_LEN), public                            :: PRECOND_TYPE = 'z_precond'
    real(kind=CUSTOM_REAL), public, dimension(:),       allocatable   :: SHORT_P
    real(kind=CUSTOM_REAL), public, dimension(:),       allocatable   :: LONG_P
    real(kind=CUSTOM_REAL), public                                    :: TW_BEFORE, TW_AFTER ! time window before and after ttp
    integer, public                                                   :: IMEAS, ITAPER
    logical, public                                                   :: USE_LOCAL_STF = .false.


    contains
    procedure :: read_tele_par_file => fwat_tele_read_par_file
  end type fwat_tele_parameters

!###################### FWAT PAR for RFAT #######################
  type fwat_rf_parameters
    real(kind=CUSTOM_REAL), public                                    :: MINDERR, RF_TSHIFT
    integer, public                                                   :: NGAUSS, MAXIT
    real(kind=CUSTOM_REAL), public, dimension(:), allocatable         :: F0

    contains
    procedure :: read_rf_par_file => fwat_rf_read_par_file
  end type fwat_rf_parameters

!###################### FWAT PAR for Local EQ FWI #######################
  type fwat_leq_parameters
    ! Other parameters need to be add.
    logical,                       public                              :: CMT3D_INV ! for cmt3d_fwd
    character(len=3),              public                              :: CMT_PAR(9)
    integer,                       public                              :: FLEX_NWIN ! number of windows made on each measurement
    integer, public                                                    :: NRCOMP
    character(len= MAX_STRING_LEN), public, dimension(:), allocatable  :: RCOMPS
    character(len= MAX_STRING_LEN), public                             :: CH_CODE, dat_coord
    character(len= MAX_STRING_LEN), public                            :: PRECOND_TYPE = 'default'
    integer, public                                                    :: NUM_FILTER,NSTEP
    real(kind=CUSTOM_REAL), public                                     :: DT
    real(kind=CUSTOM_REAL), public, dimension(:),       allocatable    :: SHORT_P
    real(kind=CUSTOM_REAL), public, dimension(:),       allocatable    :: LONG_P
    integer, public                                                    :: IMEAS, ITAPER

    contains
    procedure :: read_leq_par_file => fwat_leq_read_par_file
  end type fwat_leq_parameters

  type fwat_tomo_parameters
    logical,                       public                              :: SAVE_OUTPUT_EACH_EVENT=.false.
    logical,                       public                              :: USE_SPH_SMOOTH ! If smoothing spherical mesh
    logical,                       public                              :: USE_PDE_SMOOTH ! If use pde-based smoothing mesh 
    real(kind=CUSTOM_REAL),        public                              :: NOISE_SIGMA_H, NOISE_SIGMA_V ! for smoothing
    real(kind=CUSTOM_REAL),        public                              :: TELE_SIGMA_H, TELE_SIGMA_V ! for smoothing
    real(kind=CUSTOM_REAL),        public                              :: LEQ_SIGMA_H, LEQ_SIGMA_V ! for smoothing
    character(len= MAX_STRING_LEN), public                             :: OPT_METHOD
    integer,                       public                              :: LBFGS_M_STORE = 5
    integer,                       public                              :: ITER_START = 0
    integer,                       public                              :: NUM_STEP
    real(kind=CUSTOM_REAL),        public                              :: MAX_SLEN ! maximum step length for model update
    logical,                       public                              :: DO_LS
    real(kind=CUSTOM_REAL), public, dimension(:),   allocatable        :: STEP_LENS ! for model update and line search
    logical,                       public                              :: KERNEL_TAPER
    real(kind=CUSTOM_REAL),        public                              :: TAPER_H_SUPPRESS, TAPER_H_BUFFER
    real(kind=CUSTOM_REAL),        public                              :: TAPER_V_SUPPRESS, TAPER_V_BUFFER
    logical,                       public                              :: USE_RHO_SCALING_FWAT
    logical,                       public                              :: USE_RHO_SCALING_NOISE = .true.
    logical,                       public                              :: is_smooth = .true.
    ! Joint inversion
    ! integer,                       public                              :: INV_TYPE ! 1 for NOISE; 2 for TELE/RF; 3 for LEQ; 4 for NOISE/TELE TODO: 5 for LEQ/NOISE ...
    ! integer,                       public                              :: NUM_NOISE_SETS=1,NUM_TELE_SETS=1
    character(len=MAX_STRING_LEN), public, dimension(2)                :: NOISE_SET_RANGE
    character(len=MAX_STRING_LEN), public, dimension(2)                :: TELE_SET_RANGE
    character(len=MAX_STRING_LEN), public, dimension(2)                :: LEQ_SET_RANGE
    character(len=MAX_STRING_LEN), public, dimension(NUM_INV_TYPE)     :: INV_TYPE_NAME = (/'noise','tele ','leq  '/)
    character(len=MAX_STRING_LEN), public                              :: NORM_TYPE = 'max'
    real(kind=CUSTOM_REAL), public, dimension(NUM_INV_TYPE)            :: JOINT_WEIGHT
    logical,                public, dimension(NUM_INV_TYPE)            :: INV_TYPE ! 1 for NOISE; 2 for TELE/RF; 3 for LEQ;
    logical,                public                                     :: USE_RF
    real(kind=CUSTOM_REAL), public, dimension(2)                       :: VPVS_RATIO_RANGE = (1.3, 2.4)
    
    contains
    procedure :: read_tomo_par_file => fwat_tomo_read_par_file

  end type fwat_tomo_parameters

  contains 

subroutine fwat_acqui_read_source_set(this,model,evtset,simu_type)
  implicit none
  class(fwat_acqui),intent(inout) :: this 
  character(len=MAX_STRING_LEN),intent(in) :: model,evtset,simu_type
  integer :: myrank,ievt,ier 
  character(len=MAX_STRING_LEN) :: evtset_file,evtnm,line_junk
  real(kind=CUSTOM_REAL) :: junk_dp

  call world_rank(myrank)

  if(myrank == 0) then
    evtset_file='src_rec/sources_'//trim(evtset)//'.dat'
    open(unit=400,file=trim(evtset_file),status='old',iostat=ier)
    if (ier /=0) then
      print *,'Error could not open source subset file: ',trim(evtset_file) 
      call exit_mpi(myrank,'Error opening source subset file')
    endif
    ievt=0
    do 
      read(400,*,iostat=ier) evtnm  ! read event name
      !print*,trim(evtnm)
      if (ier /=0) exit
      ievt=ievt+1
    enddo
    close(400)
    this%nevents=ievt
    call this%malloc()
    
    open(unit=400,file=trim(evtset_file),status='old',iostat=ier)
    do ievt=1,this%nevents
      read(400, '(a)', iostat=ier) line_junk
      if (ier /=0) exit
      read(line_junk,*,iostat=ier) evtnm,this%evla(ievt),this%evlo(ievt), &
                                   this%evdp(ievt),junk_dp,&
                                   src_weight(ievt)
      if (ier /= 0) src_weight(ievt) = 1.0
      ! print *, evtnm,this%evla(ievt),this%evlo(ievt), &
      !       this%evdp(ievt),junk_dp
      if (simu_type=='noise') then 
        this%src_solution_file(ievt)='src_rec/FORCESOLUTION_'//trim(evtnm)
      elseif (simu_type=='leq') then
        this%src_solution_file(ievt)='src_rec/CMTSOLUTION_'//trim(evtnm)
      elseif (simu_type=='tele' .or. simu_type=='rf') then
        this%src_solution_file(ievt)='DATA/CMTSOLUTION'
      endif
      this%evtid_names(ievt)=trim(evtnm)
      this%station_file(ievt)='src_rec/STATIONS_'//trim(evtnm)
      this%out_fwd_path(ievt)='solver/'//trim(model)//'.'//trim(evtset)//'/'//trim(evtnm)
      this%in_dat_path(ievt)='fwat_data/'//trim(evtnm)
    enddo
    close(400)
  endif
    ! Not sure if the synchronize_all() is necessary 
  call synchronize_all()
  call bcast_all_singlei(this%nevents)
  if (myrank /= 0) then
    allocate(this%station_file(this%nevents))
    allocate(this%src_solution_file(this%nevents))
    allocate(this%evtid_names(this%nevents))
    allocate(this%out_fwd_path(this%nevents))
    allocate(this%in_dat_path(this%nevents))
    allocate(src_weight(this%nevents))
    allocate(this%evla(this%nevents))
    allocate(this%evlo(this%nevents))
    allocate(this%evdp(this%nevents))
  endif
  call bcast_all_ch_array(this%src_solution_file,this%nevents,MAX_STRING_LEN)
  call bcast_all_ch_array(this%evtid_names,this%nevents,MAX_STRING_LEN)
  call bcast_all_ch_array(this%station_file,this%nevents,MAX_STRING_LEN)
  call bcast_all_ch_array(this%out_fwd_path,this%nevents,MAX_STRING_LEN)
  call bcast_all_ch_array(this%in_dat_path,this%nevents,MAX_STRING_LEN)
  call bcast_all_cr(this%evla,this%nevents)
  call bcast_all_cr(this%evlo,this%nevents)
  call bcast_all_cr(this%evdp,this%nevents)
  ! call bcast_all_ch_array(model,1,MAX_STRING_LEN)
  call bcast_all_cr(src_weight,this%nevents)
end subroutine

subroutine fwat_acqui_free(this)
  !! free fwat_acqui class
  implicit none
  class(fwat_acqui),intent(inout) :: this 

  deallocate(this%in_dat_path,this%out_fwd_path,this%station_file)
  deallocate(this%src_solution_file,this%evtid_names)
  deallocate(this%evla,this%evlo,this%evdp,src_weight)
  
end subroutine fwat_acqui_free

subroutine fwat_acqui_malloc(this)
  !! allocate space for fwat_acqui
  implicit none
  class(fwat_acqui),intent(inout) :: this 
  integer :: nevents 

  ! local
  integer :: ier 
  ! allocate space 
  ! this%nevents = nevents
  nevents = this%nevents
  allocate(this%in_dat_path(nevents),this%station_file(nevents),&
          this%src_solution_file(nevents),this%out_fwd_path(nevents),&
          this%evtid_names(nevents),this%evla(nevents),this%evlo(nevents),&
          this%evdp(nevents),src_weight(nevents), stat=ier)
  if(ier /= 0) stop 'error in allocating fwat_acqui'

end subroutine fwat_acqui_malloc

subroutine fwat_noise_read_par_file(this)
  implicit none

  ! input variables
  class(fwat_noise_parameters),intent(inout) :: this 
  ! local
  integer            :: myrank,ipos0,ipos1,NUM_FILTER,icomp
  character(len=MAX_STRING_LEN)  :: line,keyw  
  ! get current proc
  call world_rank(myrank)

  ! begin reading params 
  if(myrank == 0) then 
    this%USE_NEAR_OFFSET=.true.
    this%IMEAS = 7
    this%ITAPER = 1
    open(666,file=FWAT_PAR_FILE) 
    do 
      read(666,'(a)',end=99) line
      if (is_blank_line(line) .or. line(1:1) == '#') cycle
      !! INDICES TO READ line -----------------------------------------------
      ipos0=index(line,':')+1
      ipos1=index(line,'#')-1
      if (ipos1 < 0 ) ipos1=len_trim(line)

      !! STORE KEYWORD ITEM -------------------------------------------------
      keyw=trim(adjustl(line(1:ipos0-2)))
               !! DIFFERENT ITEM TO READ ---------------------------------------------
      select case (trim(keyw))
        case('NOISE_NSCOMP')
          read(line(ipos0:ipos1),*) this%NSCOMP
          allocate(this%SCOMPS(this%NSCOMP))
        case('NOISE_SCOMPS')
          read(line(ipos0:ipos1),*) this%SCOMPS(1:this%NSCOMP)
        case('NOISE_NRCOMP')
          read(line(ipos0:ipos1),*) this%NRCOMP
          allocate(this%RCOMPS(this%NRCOMP))
        case('NOISE_RCOMPS')
          read(line(ipos0:ipos1),*) this%RCOMPS(1:this%NRCOMP)
          ! Mijian move dat_coord here as a public variable
          do icomp=1,this%NRCOMP
            if (trim(this%RCOMPS(icomp))=='R' .or. trim(this%RCOMPS(icomp))=='T') then
              this%dat_coord='ZRT'
            else 
              this%dat_coord='ZNE'
            endif
          enddo
        case('NOISE_CH_CODE')
          read(line(ipos0:ipos1),*) this%CH_CODE
        case('NOISE_NSTEP')
          read(line(ipos0:ipos1),*) this%NSTEP
        case('NOISE_DT')
          read(line(ipos0:ipos1),*) this%DT
        case('NOISE_NUM_FILTER')
          read(line(ipos0:ipos1),*) NUM_FILTER
          this%NUM_FILTER = NUM_FILTER
          allocate(this%SHORT_P(NUM_FILTER))
          allocate(this%LONG_P(NUM_FILTER))
          allocate(this%GROUPVEL_MIN(NUM_FILTER))!BinHe added
          allocate(this%GROUPVEL_MAX(NUM_FILTER))!BinHe added
        case('NOISE_SHORT_P')
          read(line(ipos0:ipos1),*) this%SHORT_P(:)
        case('NOISE_LONG_P')
          read(line(ipos0:ipos1),*) this%LONG_P(:)
        case('NOISE_GROUPVEL_MIN')!BinHe added
          read(line(ipos0:ipos1),*) this%GROUPVEL_MIN(:)
        case('NOISE_GROUPVEL_MAX')!BinHe added
          read(line(ipos0:ipos1),*) this%GROUPVEL_MAX(:)
        case('NOISE_ADJ_SRC_NORM')
          read(line(ipos0:ipos1),*) this%ADJ_SRC_NORM
        case('NOISE_USE_NEAR_OFFSET')!BinHe added
          read(line(ipos0:ipos1),*) this%USE_NEAR_OFFSET
        case('NOISE_SUPPRESS_EGF')
          read(line(ipos0:ipos1),*) this%SUPPRESS_EGF
        case('NOISE_PRECOND_TYPE')
          read(line(ipos0:ipos1),*) this%PRECOND_TYPE
      end select
    enddo
  endif
99 close(666) ! close par file
  call bcast_all_i(this%NSCOMP, 1,MPI_INTEGER)
  call bcast_all_i(this%NRCOMP, 1,MPI_INTEGER)
  call bcast_all_i(this%NUM_FILTER, 1,MPI_INTEGER)
  call bcast_all_i(this%IMEAS, 1,MPI_INTEGER)
  call bcast_all_i(this%ITAPER, 1,MPI_INTEGER)
  call bcast_all_singlei(this%NSTEP)
  call bcast_all_singlecr(this%DT)
  if (myrank>0) then
    allocate(this%SCOMPS(this%NSCOMP)) 
    allocate(this%RCOMPS(this%NRCOMP)) 
    allocate(this%SHORT_P(this%NUM_FILTER))
    allocate(this%LONG_P(this%NUM_FILTER))
    allocate(this%GROUPVEL_MIN(this%NUM_FILTER))!BinHe added
    allocate(this%GROUPVEL_MAX(this%NUM_FILTER))!BinHe added
  endif  
  call bcast_all_ch_array(this%SCOMPS, this%NSCOMP, MAX_STRING_LEN)
  call bcast_all_ch_array(this%RCOMPS, this%NRCOMP, MAX_STRING_LEN)
  call bcast_all_ch_array(this%dat_coord, 1, MAX_STRING_LEN)
  call bcast_all_ch_array(this%CH_CODE, 1, MAX_STRING_LEN)
  call bcast_all_cr(this%SHORT_P, this%NUM_FILTER)
  call bcast_all_cr(this%LONG_P, this%NUM_FILTER)
  call bcast_all_cr(this%GROUPVEL_MIN, this%NUM_FILTER)!BinHe added Kai changed int to real
  call bcast_all_cr(this%GROUPVEL_MAX, this%NUM_FILTER)!BinHe added Kai changed int to real
  call bcast_all_singlel(this%ADJ_SRC_NORM)
  call bcast_all_singlel(this%USE_NEAR_OFFSET)
  call bcast_all_singlel(this%SUPPRESS_EGF)
  call bcast_all_ch_array(this%PRECOND_TYPE, 1, MAX_STRING_LEN)
end subroutine fwat_noise_read_par_file

subroutine fwat_tele_read_par_file(this)
  implicit none

  ! input variables
  class(fwat_tele_parameters),intent(inout) :: this 
  ! local
  integer            :: myrank,ipos0,ipos1,NUM_FILTER,icomp
  character(len=MAX_STRING_LEN)  :: line,keyw  
  ! get current proc
  call world_rank(myrank)

  ! begin reading params 
  if(myrank == 0) then 
    this%TW_BEFORE=0.
    this%TW_AFTER=0.
    this%IMEAS = 2
    this%ITAPER = 2
    open(666,file=FWAT_PAR_FILE) 
    do 
      read(666,'(a)',end=99) line
      if (is_blank_line(line) .or. line(1:1) == '#') cycle
      !! INDICES TO READ line -----------------------------------------------
      ipos0=index(line,':')+1
      ipos1=index(line,'#')-1
      if (ipos1 < 0 ) ipos1=len_trim(line)

      !! STORE KEYWORD ITEM -------------------------------------------------
      keyw=trim(adjustl(line(1:ipos0-2)))
      !! DIFFERENT ITEM TO READ ---------------------------------------------
      select case (trim(keyw))
        case('TELE_NRCOMP')
          read(line(ipos0:ipos1),*) this%NRCOMP
          allocate(this%RCOMPS(this%NRCOMP))
        case('TELE_RCOMPS')
          read(line(ipos0:ipos1),*) this%RCOMPS(1:this%NRCOMP)
          ! Mijian move dat_coord here as a public variable
          do icomp=1,this%NRCOMP
            if (trim(this%RCOMPS(icomp))=='R' .or. trim(this%RCOMPS(icomp))=='T') then
              this%dat_coord='ZRT'
            else 
              this%dat_coord='ZNE'
            endif
          enddo
        case('TELE_CH_CODE')
          read(line(ipos0:ipos1),*) this%CH_CODE
        case('TELE_NSTEP')
          read(line(ipos0:ipos1),*) this%NSTEP
        case('TELE_DT')
          read(line(ipos0:ipos1),*) this%DT
        case('TELE_NUM_FILTER')
          read(line(ipos0:ipos1),*) NUM_FILTER
          this%NUM_FILTER = NUM_FILTER
          allocate(this%SHORT_P(NUM_FILTER))
          allocate(this%LONG_P(NUM_FILTER))
        case('TELE_SHORT_P')
          read(line(ipos0:ipos1),*) this%SHORT_P(:)
        case('TELE_LONG_P')
          read(line(ipos0:ipos1),*) this%LONG_P(:)
        case('TELE_TW_BEFORE')
          read(line(ipos0:ipos1),*) this%TW_BEFORE
        case('TELE_TW_AFTER')
          read(line(ipos0:ipos1),*) this%TW_AFTER
        case('TELE_PRECOND_TYPE')
          read(line(ipos0:ipos1),*) this%PRECOND_TYPE
        case('TELE_USE_LOCAL_STF')
          read(line(ipos0:ipos1),*) this%USE_LOCAL_STF
      end select
    enddo
  endif
99 close(666) ! close par file
  call bcast_all_singlei(this%NRCOMP)
  call bcast_all_singlei(this%NUM_FILTER)
  call bcast_all_singlei(this%IMEAS)
  call bcast_all_singlei(this%ITAPER)
  call bcast_all_singlei(this%NSTEP)
  call bcast_all_singlecr(this%DT)
  if (myrank>0) then
    allocate(this%RCOMPS(this%NRCOMP)) 
    allocate(this%SHORT_P(this%NUM_FILTER))
    allocate(this%LONG_P(this%NUM_FILTER))
  endif
  call bcast_all_ch_array(this%RCOMPS,this%NRCOMP, MAX_STRING_LEN)
  call bcast_all_ch_array(this%dat_coord, 1, MAX_STRING_LEN)
  call bcast_all_ch_array(this%CH_CODE, 1, MAX_STRING_LEN)
  call bcast_all_cr(this%SHORT_P, this%NUM_FILTER)
  call bcast_all_cr(this%LONG_P, this%NUM_FILTER)
  call bcast_all_singlecr(this%TW_BEFORE)
  call bcast_all_singlecr(this%TW_AFTER)
  call bcast_all_ch_array(this%PRECOND_TYPE, 1, MAX_STRING_LEN)
  call bcast_all_singlel(this%USE_LOCAL_STF)
end subroutine fwat_tele_read_par_file

subroutine fwat_rf_read_par_file(this)
  implicit none

  class(fwat_rf_parameters),intent(inout) :: this 
  ! local
  integer            :: myrank,ipos0,ipos1,NGAUSS
  character(len=MAX_STRING_LEN)  :: line,keyw  
  ! get current proc
  call world_rank(myrank)

  ! begin reading params 
  if(myrank == 0) then 
    open(666,file=FWAT_PAR_FILE) 
    do 
      read(666,'(a)',end=99) line
      if (is_blank_line(line) .or. line(1:1) == '#') cycle
      !! INDICES TO READ line -----------------------------------------------
      ipos0=index(line,':')+1
      ipos1=index(line,'#')-1
      if (ipos1 < 0 ) ipos1=len_trim(line)

      !! STORE KEYWORD ITEM -------------------------------------------------
      keyw=trim(adjustl(line(1:ipos0-2)))
      !! DIFFERENT ITEM TO READ ---------------------------------------------
      select case (trim(keyw))
        case('RF_NGAUSS')
          read(line(ipos0:ipos1),*) NGAUSS
          this%NGAUSS = NGAUSS
          allocate(this%F0(NGAUSS))
        case('RF_F0')
          read(line(ipos0:ipos1),*) this%F0(:)
        case('RF_MAXIT')
          read(line(ipos0:ipos1),*) this%MAXIT
        case('RF_MINDERR')
          read(line(ipos0:ipos1),*) this%MINDERR
        case('RF_TSHIFT')
          read(line(ipos0:ipos1),*) this%RF_TSHIFT
      end select
    enddo
  endif
99 close(666) ! close par file
  ! Mijian add for RF
  call bcast_all_singlei(this%NGAUSS)
  if (myrank>0) allocate(this%F0(this%NGAUSS))
  call bcast_all_cr(this%F0, this%NGAUSS)
  call bcast_all_singlecr(this%RF_TSHIFT)
  call bcast_all_singlecr(this%MINDERR)
  call bcast_all_singlei(this%MAXIT)
end subroutine fwat_rf_read_par_file

subroutine fwat_leq_read_par_file(this)
  implicit none

  class(fwat_leq_parameters),intent(inout) :: this 
  ! local
  integer            :: myrank,ipos0,ipos1,icomp
  character(len=MAX_STRING_LEN)  :: line,keyw  
  ! get current proc
  call world_rank(myrank)

  ! begin reading params 
  if(myrank == 0) then 
    this%CMT3D_INV=.false.
    this%IMEAS = 7
    this%ITAPER = 1
    open(666,file=FWAT_PAR_FILE) 
    do 
      read(666,'(a)',end=99) line
      if (is_blank_line(line) .or. line(1:1) == '#') cycle
      !! INDICES TO READ line -----------------------------------------------
      ipos0=index(line,':')+1
      ipos1=index(line,'#')-1
      if (ipos1 < 0 ) ipos1=len_trim(line)

      !! STORE KEYWORD ITEM -------------------------------------------------
      keyw=trim(adjustl(line(1:ipos0-2)))
      !! DIFFERENT ITEM TO READ ---------------------------------------------
      select case (trim(keyw))
        case('LEQ_CMT3D_INV')
          read(line(ipos0:ipos1),*) this%CMT3D_INV
        case('LEQ_NRCOMP')
          read(line(ipos0:ipos1),*) this%NRCOMP
          allocate(this%RCOMPS(this%NRCOMP))
        case('LEQ_RCOMPS')
          read(line(ipos0:ipos1),*) this%RCOMPS(1:this%NRCOMP)
          ! Mijian move dat_coord here as a public variable
          do icomp=1,this%NRCOMP
            if (trim(this%RCOMPS(icomp))=='R' .or. trim(this%RCOMPS(icomp))=='T') then
              this%dat_coord='ZRT'
            else 
              this%dat_coord='ZNE'
            endif
          enddo
        case('LEQ_CH_CODE')
          read(line(ipos0:ipos1),*) this%CH_CODE
        case('LEQ_NSTEP')
          read(line(ipos0:ipos1),*) this%NSTEP
        case('LEQ_DT')
          read(line(ipos0:ipos1),*) this%DT
        case('LEQ_NUM_FILTER')
          read(line(ipos0:ipos1),*) NUM_FILTER
          this%NUM_FILTER = NUM_FILTER
          allocate(this%SHORT_P(NUM_FILTER))
          allocate(this%LONG_P(NUM_FILTER))
        case('LEQ_SHORT_P')
          read(line(ipos0:ipos1),*) this%SHORT_P(:)
        case('LEQ_LONG_P')
          read(line(ipos0:ipos1),*) this%LONG_P(:)
        case('LEQ_IMEAS')
          read(line(ipos0:ipos1),*) this%IMEAS
        case('LEQ_ITAPER')
          read(line(ipos0:ipos1),*) this%ITAPER
        case('LEQ_PRECOND_TYPE')
          read(line(ipos0:ipos1),*) this%PRECOND_TYPE
      end select
    enddo
  endif
99 close(666) ! close par file
  call bcast_all_singlel(this%CMT3D_INV)
  call bcast_all_singlei(this%NUM_FILTER)
  call bcast_all_singlei(this%IMEAS)
  call bcast_all_singlei(this%ITAPER)
  call bcast_all_singlei(this%NSTEP)
  call bcast_all_singlecr(this%DT)
  if (myrank /= 0) then
    allocate(this%SHORT_P(this%NUM_FILTER))
    allocate(this%LONG_P(this%NUM_FILTER))
  endif
  call bcast_all_cr(this%SHORT_P, this%NUM_FILTER)
  call bcast_all_cr(this%LONG_P, this%NUM_FILTER)
  call bcast_all_ch_array(this%PRECOND_TYPE, 1, MAX_STRING_LEN)

end subroutine fwat_leq_read_par_file

subroutine fwat_tomo_read_par_file(this)
  implicit none

  class(fwat_tomo_parameters),intent(inout) :: this 
  ! local
  integer            :: myrank,ipos0,ipos1
  character(len=MAX_STRING_LEN)  :: line,keyw  
  ! get current proc
  call world_rank(myrank)

  ! begin reading params 
  if(myrank == 0) then 
    open(666,file=FWAT_PAR_FILE) 
    do 
      read(666,'(a)',end=99) line
      if (is_blank_line(line) .or. line(1:1) == '#') cycle
      !! INDICES TO READ line -----------------------------------------------
      ipos0=index(line,':')+1
      ipos1=index(line,'#')-1
      if (ipos1 < 0 ) ipos1=len_trim(line)

      !! STORE KEYWORD ITEM -------------------------------------------------
      keyw=trim(adjustl(line(1:ipos0-2)))
      !! DIFFERENT ITEM TO READ ---------------------------------------------
      select case (trim(keyw))
        case('SAVE_OUTPUT_EACH_EVENT')
          read(line(ipos0:ipos1),*) this%SAVE_OUTPUT_EACH_EVENT
        case('USE_SPH_SMOOTH')
          read(line(ipos0:ipos1),*) this%USE_SPH_SMOOTH
        case('USE_PDE_SMOOTH')
          read(line(ipos0:ipos1),*) this%USE_PDE_SMOOTH
        case('OPT_METHOD')
          read(line(ipos0:ipos1),*) this%OPT_METHOD
        case('LBFGS_M_STORE')
          read(line(ipos0:ipos1),*) this%LBFGS_M_STORE
        case('ITER_START')
          read(line(ipos0:ipos1),*) this%ITER_START
        case('DO_LS')
          read(line(ipos0:ipos1),*) this%DO_LS
        case('MAX_SLEN')
          read(line(ipos0:ipos1),*) this%MAX_SLEN
        case('NUM_STEP')
          read(line(ipos0:ipos1),*) this%NUM_STEP
          allocate(this%STEP_LENS(this%NUM_STEP))
        case('STEP_LENS')
          read(line(ipos0:ipos1),*) this%STEP_LENS(:)
        case('KERNEL_TAPER')
          read(line(ipos0:ipos1),*) this%KERNEL_TAPER
        case('TAPER_H_SUPPRESS')
          read(line(ipos0:ipos1),*) this%TAPER_H_SUPPRESS
        case('TAPER_H_BUFFER')
          read(line(ipos0:ipos1),*) this%TAPER_H_BUFFER
        case('TAPER_V_SUPPRESS')
          read(line(ipos0:ipos1),*) this%TAPER_V_SUPPRESS
        case('TAPER_V_BUFFER')
          read(line(ipos0:ipos1),*) this%TAPER_V_BUFFER
        case('USE_RHO_SCALING_FWAT')
          read(line(ipos0:ipos1),*) this%USE_RHO_SCALING_FWAT
        case('USE_RHO_SCALING_NOISE')
          read(line(ipos0:ipos1),*) this%USE_RHO_SCALING_NOISE
        case('IS_SMOOTH')
          read(line(ipos0:ipos1),*) this%IS_SMOOTH
        case('NOISE_SIGMA_H')
          read(line(ipos0:ipos1),*) this%NOISE_SIGMA_H
        case('NOISE_SIGMA_V')
          read(line(ipos0:ipos1),*) this%NOISE_SIGMA_V
        case('TELE_SIGMA_H')
          read(line(ipos0:ipos1),*) this%TELE_SIGMA_H
        case('TELE_SIGMA_V')
          read(line(ipos0:ipos1),*) this%TELE_SIGMA_V
        case('LEQ_SIGMA_H')
          read(line(ipos0:ipos1),*) this%LEQ_SIGMA_H
        case('LEQ_SIGMA_V')
          read(line(ipos0:ipos1),*) this%LEQ_SIGMA_V  
        case('INV_TYPE')
          read(line(ipos0:ipos1),*) this%INV_TYPE(:)
        case('NORM_TYPE')
          read(line(ipos0:ipos1),*) this%NORM_TYPE
        case('NOISE_SET_RANGE')
          read(line(ipos0:ipos1),*) this%NOISE_SET_RANGE(:)
        ! case('NUM_TELE_SETS')
        !   read(line(ipos0:ipos1),*) this%NUM_TELE_SETS
        !   allocate(this%TELE_SET_NAMES(this%NUM_TELE_SETS))
        case('TELE_SET_RANGE')
          read(line(ipos0:ipos1),*) this%TELE_SET_RANGE(:)
        case('LEQ_SET_RANGE')
          read(line(ipos0:ipos1),*) this%LEQ_SET_RANGE(:)
        case('USE_RF')
          read(line(ipos0:ipos1),*) this%USE_RF
          if (this%USE_RF) this%INV_TYPE_NAME(2)='rf'
        case('JOINT_WEIGHT')
          read(line(ipos0:ipos1),*) this%JOINT_WEIGHT(:)
        case('VPVS_RATIO_RANGE')
          read(line(ipos0:ipos1),*) this%VPVS_RATIO_RANGE(:)
      end select
    enddo
  endif
99 close(666) ! close par file
  call bcast_all_singlei(this%NUM_STEP)
  ! call bcast_all_singlei(this%INV_TYPE)
  call bcast_all_singlei(this%LBFGS_M_STORE)
  call bcast_all_singlei(this%ITER_START)
  call bcast_all_singlel(this%SAVE_OUTPUT_EACH_EVENT)
  call bcast_all_singlel(this%USE_SPH_SMOOTH)
  call bcast_all_singlel(this%USE_PDE_SMOOTH)
  call bcast_all_singlel(this%USE_RHO_SCALING_FWAT)
  call bcast_all_singlel(this%USE_RHO_SCALING_NOISE)
  call bcast_all_singlel(this%IS_SMOOTH)
  call bcast_all_singlel(this%USE_RF)
  call bcast_all_singlecr(this%NOISE_SIGMA_H)
  call bcast_all_singlecr(this%NOISE_SIGMA_V)
  call bcast_all_singlecr(this%TELE_SIGMA_H)
  call bcast_all_singlecr(this%TELE_SIGMA_V)
  call bcast_all_singlecr(this%LEQ_SIGMA_H)
  call bcast_all_singlecr(this%LEQ_SIGMA_V)
  call bcast_all_ch_array(this%OPT_METHOD, 1, MAX_STRING_LEN)
  call bcast_all_singlel(this%DO_LS)
  call bcast_all_singlecr(this%MAX_SLEN)
  if (myrank /= 0) then
    allocate(this%STEP_LENS(this%NUM_STEP))
  endif
  call bcast_all_cr(this%STEP_LENS,this%NUM_STEP)
  call bcast_all_cr(this%JOINT_WEIGHT,NUM_INV_TYPE)
  call bcast_all_l_array(this%INV_TYPE,NUM_INV_TYPE)
  call bcast_all_ch_array(this%INV_TYPE_NAME, NUM_INV_TYPE, MAX_STRING_LEN)
  call bcast_all_ch_array(this%TELE_SET_RANGE, 2, MAX_STRING_LEN)
  call bcast_all_ch_array(this%NOISE_SET_RANGE, 2, MAX_STRING_LEN)
  call bcast_all_ch_array(this%LEQ_SET_RANGE, 2, MAX_STRING_LEN)
  call bcast_all_cr(this%VPVS_RATIO_RANGE,2)
  call bcast_all_singlel(this%KERNEL_TAPER)
  call bcast_all_singlecr(this%TAPER_H_SUPPRESS)
  call bcast_all_singlecr(this%TAPER_H_BUFFER)
  call bcast_all_singlecr(this%TAPER_V_SUPPRESS)
  call bcast_all_singlecr(this%TAPER_V_BUFFER)

end subroutine fwat_tomo_read_par_file

subroutine fwat_common_read_par_file()
  implicit none

  integer            :: myrank,ipos0,ipos1
  character(len=MAX_STRING_LEN)  :: line,keyw  
  ! get current proc
  call world_rank(myrank)

  ! begin reading params 
  if(myrank == 0) then 
    open(666,file=FWAT_PAR_FILE) 
    do 
      read(666,'(a)',end=99) line
      if (is_blank_line(line) .or. line(1:1) == '#') cycle
      !! INDICES TO READ line -----------------------------------------------
      ipos0=index(line,':')+1
      ipos1=index(line,'#')-1
      if (ipos1 < 0 ) ipos1=len_trim(line)

      !! STORE KEYWORD ITEM -------------------------------------------------
      keyw=trim(adjustl(line(1:ipos0-2)))
      !! DIFFERENT ITEM TO READ ---------------------------------------------
      select case (trim(keyw))
        case('VERBOSE_MODE')
          read(line(ipos0:ipos1),*) VERBOSE_MODE
      end select
    enddo
  endif
99 close(666) ! close par file
  call bcast_all_singlel(VERBOSE_MODE)
end subroutine fwat_common_read_par_file

logical function is_blank_line(line)
  character(len=MAX_STRING_LEN), intent(in) :: line
  is_blank_line=.false.
  if (len(trim(adjustl(line))) == 0) is_blank_line = .true.
  if (INDEX(trim(adjustl(line)),'#') == 1) is_blank_line = .true.
end function is_blank_line

end module fullwave_adjoint_tomo_par

!################################################# measure_adj ####################################################################

module ma_constants

  ! number of entries in window_chi output file
  integer, parameter :: N_MEASUREMENT = 5
  integer, parameter :: NCHI = 3*(N_MEASUREMENT-1) + 8

  ! constants
  double precision, parameter :: PI = 3.141592653589793d+00
  double precision, parameter :: TWOPI = 2.0 * PI
  complex(kind=8), parameter :: CCI = cmplx(0.,1.)
  double precision, parameter :: LARGE_VAL = 1.0d8

  ! FFT parameters
  integer, parameter :: LNPT = 16, NPT = 2**LNPT, NDIM = 80000
  double precision, parameter :: FORWARD_FFT = 1.0
  double precision, parameter :: REVERSE_FFT = -1.0

  ! phase correction control parameters, set this between (PI, 2PI),
  ! use a higher value for conservative phase wrapping
  double precision, parameter :: PHASE_STEP = 1.5 * PI

  ! filter parameters for xapiir bandpass subroutine (filter type is BP)
  ! (These should match the filter used in pre-processing.)
  double precision, parameter :: TRBDNDW = 0.3
  double precision, parameter :: APARM = 30.
  integer, parameter :: IORD = 4
  integer, parameter :: PASSES = 2

  ! takes waveform of first trace dat_dtw, without taking the difference waveform to the second trace syn_dtw
  ! this is useful to cissor out later reflections which appear in data (no synthetics needed)
  logical, parameter :: NO_WAVEFORM_DIFFERENCE = .false.

  ! constructs adjoint sources for a "ray density" kernel, where all misfits are equal to one
  logical, parameter :: DO_RAY_DENSITY_SOURCE = .false.

end module ma_constants

module ma_variables

  use ma_constants
!
! multi-taper measurements
!
! Ying Zhou: The fit between the recovered data and the data can be improved
! by either increasing the window width (HWIN above) or by decreasing NPI.
! In her experience, NPI = 2.5 is good for noisy data.
! For synthetic data, we can use a lower NPI.
! number of tapers should be fixed as twice NPI -- see Latex notes
!
! See write_par_file.pl and measure_adj.f90

  character(len=512) :: OUT_DIR

  double precision, public :: TLONG, TSHORT
  double precision, public :: WTR, NPI, DT_FAC, ERR_FAC, DT_MAX_SCALE, NCYCLE_IN_WINDOW
  !double precision :: BEFORE_QUALITY, AFTER_QUALITY, BEFORE_TSHIFT, AFTER_TSHIFT
  double precision, public :: TSHIFT_MIN, TSHIFT_MAX, DLNA_MIN, DLNA_MAX, CC_MIN
  double precision, public :: DT_SIGMA_MIN, DLNA_SIGMA_MIN

  integer, public :: ntaper, ipwr_t, ipwr_w, ERROR_TYPE
  integer, public :: imeas0, imeas, itaper, is_mtm0, is_mtm

  logical, public :: DISPLAY_DETAILS,OUTPUT_MEASUREMENT_FILES,RUN_BANDPASS,COMPUTE_ADJOINT_SOURCE,USE_PHYSICAL_DISPERSION
  
  contains

end module ma_variables


module ma_weighting

! module for weighting/normalizing measurements

  logical,parameter :: DO_WEIGHTING = .false.

  ! transverse, radial and vertical weights
  double precision :: weight_T, weight_R, weight_Z
  ! body waves: number of picks on vertical, radial and transverse component
  double precision :: num_P_SV_V,num_P_SV_R,num_SH_T
  ! surface waves: number of pick on vertical, radial and transverse
  double precision :: num_Rayleigh_V,num_Rayleigh_R,num_Love_T

  ! typical surface wave speed in km/s, to calculate surface wave arrival times
  ! Love waves faster than Rayleigh
  double precision, parameter :: surface_vel = 4.0

  ! wave type pick
  integer, parameter :: P_SV_V = 1
  integer, parameter :: P_SV_R = 2
  integer, parameter :: SH_T = 3
  integer, parameter :: Rayleigh_V = 4
  integer, parameter :: Rayleigh_R = 5
  integer, parameter :: Love_T = 6

end module ma_weighting

!################################################# flexwin ####################################################################

module user_parameters
 
  use ma_constants, only: TRBDNDW,APARM,IORD,PASSES,NDIM
  !===================================================================
  ! filter parameters for xapiir subroutine (filter type is BP)
  !double precision, parameter :: TRBDNDW = 0.3
  !double precision, parameter :: APARM = 30.0
  !integer, parameter :: IORD = 5
  !integer, parameter :: PASSES = 2

  ! -------------------------------------------------------------
  ! array dimensions
  ! note that some integer arrays (iM,iL,iR) are NWINDOWS * NWINDOWS
  ! THESE SHOULD PROBABLY BE USER PARAMETERS, SINCE THEY WILL AFFECT
  ! THE SPEED OF THE PROGRAM (ALTHOUGH NOT THE OUTPUT).
  integer, parameter :: NWINDOWS = 2500

  ! -------------------------------------------------------------
  ! miscellaneous - do not modify!
  ! -------------------------------------------------------------

  ! mathematical constants
  double precision, parameter :: PI = 3.1415926535897
  double precision, parameter :: E  = 2.7182818284590

  ! filter types
  integer, parameter :: HANNING = 1
  integer, parameter :: HAMMING = 2
  integer, parameter :: COSINE  = 3

  ! -------------------------------------------------------------

end module user_parameters

