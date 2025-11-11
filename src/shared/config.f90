module fwat_constants
  use constants
  implicit none
  include "fwat_config.h"

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
  character(len=MAX_STRING_LEN), dimension(5), parameter :: KERNEL_AZI_ANI = ['alpha', 'beta ', 'rhop ', 'gcp  ', 'gsp  ']

  ! model types
  character(len=MAX_STRING_LEN), dimension(3), parameter :: MODEL_NAME_ISO = ['vp ', 'vs ', 'rho']
  character(len=MAX_STRING_LEN), dimension(5), parameter :: MODEL_NAME_AZI_ANI = ['vp ', 'vs ', 'rho', 'gcp', 'gsp']

  real(kind=cr),parameter :: THRESHOLD_HESS = 1.e-3, RHO_SCALING_FAC = 0.33
  integer, parameter :: NUM_INV_TYPE = 3

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
  integer, parameter :: FORWARD_SAVE = 4
  integer, parameter :: ADJOINT_ONLY = 5
  character(len=MAX_STRING_LEN), dimension(3), parameter :: INV_TYPE_NAMES = [SIMU_TYPE_NOISE, SIMU_TYPE_TELE, SIMU_TYPE_LEQ]
  integer, parameter :: WIN_SELECTOR_TYPE = 1
  integer, parameter :: WIN_GROUPVEL_TYPE = 2
  integer, parameter :: WIN_ARRIVAL_TYPE = 3
  integer, parameter :: PHASE_NAME_LEN = 8

  ! math
  real(kind=cr), parameter :: EARTH_RADIUS = 6371.0d0
  real(kind=cr), parameter :: deg2rad = PI / 180.0d0
  real(kind=cr), parameter :: rad2deg = 180.0d0 / PI
  real(kind=cr), parameter :: km2deg = 1.0d0 /(EARTH_RADIUS*pi/180.0d0)

  ! meas_adj
  integer, parameter :: NITER = 200
  logical :: supp_stf = .false.
  real(kind=dp), parameter :: TRBDNDW = 0.3
  real(kind=dp), parameter :: APARM = 30.
  integer, parameter :: IORD = 2
  integer, parameter :: PASSES = 2
  real, parameter :: MAX_TIME_SHIFT = 10.0 ! in seconds

  ! Measurement types
  integer, parameter :: IMEAS_WAVEFORM = 1
  integer, parameter :: IMEAS_WAVEFORM_CONV = 2
  integer, parameter :: IMEAS_RF = 3
  integer, parameter :: IMEAS_EXP_PHASE = 4
  integer, parameter :: IMEAS_CCC = 5
  integer, parameter :: IMEAS_CC_TT = 11
  integer, parameter :: IMEAS_CC_DLNA = 12
  integer, parameter :: IMEAS_CC_TT_MT = 13
  integer, parameter :: IMEAS_CC_DLNA_MT = 14

end module fwat_constants

module config
  use fwat_constants

  integer :: worldrank, worldsize
  integer :: noderank, nodesize
  integer, dimension(:,:), allocatable :: rank_map

  logical :: single_run = .false.
  integer :: event_index = 0, run_mode = 0, compress_level
  
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
  ! real(kind=cr), dimension(:), allocatable :: xstore_fwat, ystore_fwat, zstore_fwat
  ! integer, dimension(:,:,:,:), allocatable :: ibool_fwat

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