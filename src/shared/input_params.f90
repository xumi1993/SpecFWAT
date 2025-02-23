module imput_params
  use config

  implicit none

  type acqui_params
    integer               :: nevents ! # of events
    character(len= MAX_STRING_LEN), dimension(:),       allocatable  :: station_file
    character(len= MAX_STRING_LEN), dimension(:),       allocatable  :: src_solution_file
    character(len= MAX_STRING_LEN), dimension(:),       allocatable  :: evtid_names
    character(len= MAX_STRING_LEN), dimension(:),       allocatable  :: out_fwd_path
    character(len= MAX_STRING_LEN), dimension(:),       allocatable  :: in_dat_path
    real(kind=CUSTOM_REAL), dimension(:), allocatable                :: evla,evlo,evdp

    ! contains
    ! procedure :: free => fwat_acqui_free
    ! procedure :: malloc => fwat_acqui_malloc
    ! procedure :: read_source_set => fwat_acqui_read_source_set
  end type acqui_params

  type rf_params
    real(kind=CUSTOM_REAL), public                                    :: MINDERR, RF_TSHIFT
    integer, public                                                   :: NGAUSS, MAXIT
    real(kind=CUSTOM_REAL), public, dimension(:), allocatable         :: F0
  end type rf_params

  type sim_params
    integer :: NRCOMP, NSCOMP, NUM_FILTER, NSTEP, IMEAS, ITAPER, PRECOND_TYPE
    character(len= MAX_STRING_LEN), dimension(:), allocatable :: RCOMPS, SCOMPS
    character(len= MAX_STRING_LEN) :: CH_CODE, dat_coord
    real(kind=cr) :: DT, TW_BEFORE, TW_AFTER
    real(kind=cr), dimension(:), allocatable :: SHORT_P, LONG_P, GROUPVEL_MIN, GROUPVEL_MAX
    logical :: USE_NEAR_OFFSET, ADJ_SRC_NORM, SUPPRESS_EGF, USE_LOCAL_STF
    type(rf_params) :: rf
  end type sim_params


end module imput_params