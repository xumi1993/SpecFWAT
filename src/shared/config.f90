module fwat_constants
  ! include "constants.h"
  use constants

  integer, parameter :: cr = SIZE_REAL
  integer, parameter :: dp = SIZE_DOUBLE
  integer, parameter :: targ = 1000

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
  character(len=MAX_STRING_LEN), parameter :: HESS_PREFIX = "hess_inv"
  character(len=MAX_STRING_LEN), parameter :: OUTPUT_PATH = "OUTPUT_FILES"
  character(len=MAX_STRING_LEN), parameter :: ADJOINT_PATH = "SEM"
  character(len=MAX_STRING_LEN), parameter :: EKERNEL_PATH = "EKERNEL"

  ! Kernel types
  character(len=MAX_STRING_LEN), dimension(3), parameter :: KERNEL_ISO = ['alpha', 'beta ', 'rhop ']
  character(len=MAX_STRING_LEN), dimension(3), parameter :: KERNEL_AZI_ANI = ['L ', 'Gc', 'Gs']

  ! model types
  character(len=MAX_STRING_LEN), dimension(3), parameter :: MODEL_ISO = ['vp ', 'vs ', 'rho']
  character(len=MAX_STRING_LEN), dimension(3), parameter :: MODEL_AZI_ANI = ['L ', 'Gc', 'Gs']

  real(kind=cr),parameter :: THRESHOLD_HESS = 1.e-3, RHO_SCALING_FAC = 0.33
  integer, parameter :: NUM_INV_TYPE = 2

  ! Preconditioner parameters
  integer, parameter :: DEFAULT_PRECOND = 1
  integer, parameter :: Z_PRECOND = 2
  integer, parameter :: Z_SQRT_PRECOND = 3
  character(len=MAX_STRING_LEN), parameter :: SIMU_TYPE_NOISE = 'noise'
  character(len=MAX_STRING_LEN), parameter :: SIMU_TYPE_TELE = 'tele'
  character(len=MAX_STRING_LEN), parameter :: SIMU_TYPE_LEQ = 'leq'
  integer, parameter :: FORWARD_ONLY = 1
  integer, parameter :: FORWARD_MEASADJ = 2
  integer, parameter :: FORWARD_ADJOINT = 3
  character(len=MAX_STRING_LEN), dimension(3), parameter :: INV_TYPE_NAMES = [SIMU_TYPE_NOISE, SIMU_TYPE_TELE, SIMU_TYPE_LEQ]


  ! math
  real(kind=cr), parameter :: deg2rad = PI / 180.0
  real(kind=cr), parameter :: rad2deg = 180.0 / PI
  real(kind=cr), parameter :: km2deg = 1.0 /(6371.0d0*pi/180.0d0)

  ! meas_adj
  integer, parameter :: NITER = 200

end module fwat_constants

module config
  use fwat_constants, only: cr, MAX_STRING_LEN

  integer :: worldrank, worldsize
  integer :: noderank, nodesize
  integer, dimension(:,:), allocatable :: rank_map

  logical :: single_run = .false.
  integer :: event_index = 0, run_mode, compress_level
  
  character(len=MAX_STRING_LEN) :: mesh_par_file
  character(len=MAX_STRING_LEN) :: dat_coord, local_path_backup, local_path_fwat, output_files_backup
  character(len=MAX_STRING_LEN) :: simu_type, dat_type, model_name, model_name_ls
  real(kind=cr), dimension(:), allocatable :: rmassx_copy,rmassy_copy,rmassz_copy,rmass_copy,rmass_acoustic_copy
  real(kind=cr), dimension(:,:,:,:), allocatable :: sum_rho_kl 
  real(kind=cr), dimension(:,:,:,:,:), allocatable :: sum_cijkl_kl 
  real(kind=cr), dimension(:,:,:,:), allocatable :: sum_mu_kl 
  real(kind=cr), dimension(:,:,:,:), allocatable :: sum_kappa_kl 
  real(kind=cr), dimension(:,:,:,:), allocatable :: sum_hess_kl

  ! output
  logical :: is_output_preproc, is_output_adj_src, IS_OUTPUT_EVENT_KERNEL,&
             is_output_sum_kernel, is_output_direction, is_output_inv_grid

  ! mesh
  real(kind=cr) :: x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                   elemsize_min_glob,elemsize_max_glob, &
                   distance_min_glob,distance_max_glob
  integer :: NSPEC_FWAT, NGLOB_FWAT
  real(kind=cr), dimension(:), allocatable :: xstore_fwat, ystore_fwat, zstore_fwat
  integer, dimension(:,:,:,:), allocatable :: ibool_fwat

  ! post
  logical :: is_joint = .false.
  integer :: nkernel, parameter_type
  character(len=MAX_STRING_LEN) :: kernel_type
  character(len=MAX_STRING_LEN), dimension(:), allocatable :: kernel_names, parameter_names
  character(len=MAX_STRING_LEN) :: model_prev, model_next, model_start, model_current
  real(kind=cr) :: step_len

contains

  function select_iproc_for_rec(iproc)
    integer, intent(in) :: iproc
    integer :: select_iproc_for_rec
    
    select_iproc_for_rec = mod(iproc-1, worldsize)
  end function select_iproc_for_rec

  function get_num_recs_per_proc(total_recs, iproc) result(nrecs)
    integer, intent(in) :: total_recs, iproc
    integer :: nrecs, base, remainder

    base = total_recs / worldsize
    remainder = mod(total_recs, worldsize)
    if (iproc < remainder) then
      nrecs = base + 1
    else
      nrecs = base
    endif
  end function get_num_recs_per_proc

  function select_global_id_for_rec(irec_local)
    integer :: irec_local
    integer :: select_global_id_for_rec

    select_global_id_for_rec = worldrank + 1 + (irec_local - 1) * worldsize
  end function select_global_id_for_rec

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
  integer, parameter :: LNPT = 16, NPT = 2**LNPT, NDIM_MA = 80000
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
  integer, parameter :: MAX_STR_CHI = 30

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

  logical, public :: DISPLAY_DETAILS=.false.,OUTPUT_MEASUREMENT_FILES=.false.,&
                     RUN_BANDPASS=.false.,COMPUTE_ADJOINT_SOURCE=.true.,USE_PHYSICAL_DISPERSION=.false.
  
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