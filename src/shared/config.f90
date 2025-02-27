module fwat_constants
  include "constants.h"

  integer, parameter :: cr = SIZE_REAL
  integer, parameter :: dp = SIZE_DOUBLE

  character(len=MAX_STRING_LEN), parameter :: FWAT_PAR_FILE = trim(IN_DATA_FILES)//"fwat_params.yml"

  ! Path to the data directory
  character(len=MAX_STRING_LEN), parameter :: DATA_DIR = "fwat_data"
  character(len=MAX_STRING_LEN), parameter :: SRC_REC_DIR = "src_rec"
  character(len=MAX_STRING_LEN), parameter :: SOLVER_DIR = "solver"
  character(len=MAX_STRING_LEN), parameter :: MISFITS_DIR = "misfits"
  character(len=MAX_STRING_LEN), parameter :: OPT_DIR = "optimize"

  ! Prefix for the src_rec files
  character(len=MAX_STRING_LEN), parameter :: FKMODEL_PREFIX = "FKmodel"
  character(len=MAX_STRING_LEN), parameter :: SRC_PREFIX = "sources"
  character(len=MAX_STRING_LEN), parameter :: STATIONS_PREFIX = "STATIONS"
  character(len=MAX_STRING_LEN), parameter :: CMTSOLUTION_PREFIX = "CMTSOLUTION"
  character(len=MAX_STRING_LEN), parameter :: FORCESOLUTION_PREFIX = "FORCESOLUTION"
  character(len=MAX_STRING_LEN), parameter :: OUTPUT_PATH = "OUTPUT_FILES"
  character(len=MAX_STRING_LEN), parameter :: ADJOINT_PATH = "SEM"
  character(len=MAX_STRING_LEN), parameter :: EKERNEL_PATH = "EKERNEL"

  ! Preconditioner parameters
  integer, parameter :: DEFAULT_PRECOND = 1
  integer, parameter :: Z_PRECOND = 2
  character(len=MAX_STRING_LEN), parameter :: SIMU_TYPE_NOISE = 'noise'
  character(len=MAX_STRING_LEN), parameter :: SIMU_TYPE_TELE = 'tele'
  integer, parameter :: FORWARD_ONLY = 1
  integer, parameter :: FORWARD_MEASADJ = 2
  integer, parameter :: FORWARD_ADJOINT = 3

end module fwat_constants

module config
  use fwat_constants, only: cr, MAX_STRING_LEN

  integer :: worldrank, worldsize
  integer :: noderank, nodesize
  integer, dimension(:,:), allocatable :: rank_map

  logical :: single_run = .false.
  integer :: event_index = 0
  
  character(len=MAX_STRING_LEN) :: dat_coord, local_path_backup
  character(len=MAX_STRING_LEN) :: simu_type, dat_type, model_name
  real(kind=cr), dimension(:), allocatable :: rmassx_copy,rmassy_copy,rmassz_copy,rmass_copy,rmass_acoustic_copy
  real(kind=cr), dimension(:,:,:,:), allocatable :: sum_rho_kl 
  real(kind=cr), dimension(:,:,:,:,:), allocatable :: sum_cijkl_kl 
  real(kind=cr), dimension(:,:,:,:), allocatable :: sum_mu_kl 
  real(kind=cr), dimension(:,:,:,:), allocatable :: sum_kappa_kl 
  real(kind=cr), dimension(:,:,:,:), allocatable :: sum_hess_kl 

end module config

module ma_constants
  use fwat_constants, only: dp, PI

  ! number of entries in window_chi output file
  integer, parameter :: N_MEASUREMENT = 5
  integer, parameter :: NCHI = 3*(N_MEASUREMENT-1) + 8

  ! constants
  real(kind=dp), parameter :: TWOPI = 2.0 * PI
  complex(kind=dp), parameter :: CCI = cmplx(0.,1.)
  real(kind=dp), parameter :: LARGE_VAL = 1.0d8

  ! FFT parameters
  integer, parameter :: LNPT = 16, NPT = 2**LNPT, NDIM = 80000
  real(kind=dp), parameter :: FORWARD_FFT = 1.0
  real(kind=dp), parameter :: REVERSE_FFT = -1.0

  ! phase correction control parameters, set this between (PI, 2PI),
  ! use a higher value for conservative phase wrapping
  real(kind=dp), parameter :: PHASE_STEP = 1.5 * PI

  ! filter parameters for xapiir bandpass subroutine (filter type is BP)
  ! (These should match the filter used in pre-processing.)
  real(kind=dp), parameter :: TRBDNDW = 0.3
  real(kind=dp), parameter :: APARM = 30.
  integer, parameter :: IORD = 2
  integer, parameter :: PASSES = 2

  ! takes waveform of first trace dat_dtw, without taking the difference waveform to the second trace syn_dtw
  ! this is useful to cissor out later reflections which appear in data (no synthetics needed)
  logical, parameter :: NO_WAVEFORM_DIFFERENCE = .false.

  ! constructs adjoint sources for a "ray density" kernel, where all misfits are equal to one
  logical, parameter :: DO_RAY_DENSITY_SOURCE = .false.

end module ma_constants

module ma_variables
  use fwat_constants, only: MAX_STRING_LEN
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

  character(len=MAX_STRING_LEN) :: OUT_DIR

  real(kind=dp), public :: TLONG, TSHORT
  real(kind=dp), public :: WTR, NPI, DT_FAC, ERR_FAC, DT_MAX_SCALE, NCYCLE_IN_WINDOW
  !double precision :: BEFORE_QUALITY, AFTER_QUALITY, BEFORE_TSHIFT, AFTER_TSHIFT
  real(kind=dp), public :: TSHIFT_MIN, TSHIFT_MAX, DLNA_MIN, DLNA_MAX, CC_MIN
  real(kind=dp), public :: DT_SIGMA_MIN, DLNA_SIGMA_MIN

  integer, public :: ntaper, ipwr_t, ipwr_w, ERROR_TYPE
  integer, public :: imeas0, imeas, itaper, is_mtm0, is_mtm

  logical, public :: DISPLAY_DETAILS,OUTPUT_MEASUREMENT_FILES,RUN_BANDPASS,COMPUTE_ADJOINT_SOURCE,USE_PHYSICAL_DISPERSION
  
  contains

end module ma_variables