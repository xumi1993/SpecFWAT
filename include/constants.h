!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! setup/constants.h.  Generated from constants.h.in by configure.

!
! solver in single or double precision depending on the machine (4 or 8 bytes)
!
! ALSO CHANGE FILE precision.h ACCORDINGLY
!
  integer, parameter :: SIZE_REAL = 4
  integer, parameter :: SIZE_DOUBLE = 8

! usually the size of integer and logical variables is the same as regular single-precision real variable
  integer, parameter :: SIZE_INTEGER = SIZE_REAL
  integer, parameter :: SIZE_LOGICAL = SIZE_REAL

! set to SIZE_REAL to run in single precision
! set to SIZE_DOUBLE to run in double precision (increases memory size by 2)
  integer, parameter :: CUSTOM_REAL = SIZE_REAL
  integer, parameter :: CUSTOM_DOUBLE = SIZE_DOUBLE

!----------- parameters that can be changed by the user -----------

! for validation of undoing of attenuation versus an exact solution saved to disk
! should never be needed any more, it was developed and used for Figure 3 of
! D. Komatitsch, Z. Xie, E. Bozdag, E. Sales de Andrade, D. Peter, Q. Liu and J. Tromp,
! Anelastic sensitivity kernels with parsimonious storage for adjoint tomography and full waveform inversion,
! Geophysical Journal International, vol. 206(3), p. 1467-1478, doi: 10.1093/gji/ggw224 (2016).
  logical, parameter :: EXACT_UNDOING_TO_DISK = .false.
! ID of the huge file in which we dump all the time steps of the simulation
  integer, parameter :: IFILE_FOR_EXACT_UNDOING = 244

! set to .false.  if running on a Beowulf-type machine with local disks
! set to .true. if running on a shared-memory machine with common file system
! if running on a Beowulf, also modify name of nodes in filter_machine_file.f90
  logical, parameter :: LOCAL_PATH_IS_ALSO_GLOBAL = .true.

! maximum length of strings used for paths, reading from files, etc.
  integer, parameter :: MAX_STRING_LEN = 512
  integer, parameter :: MAX_NAME_LEN = 30

!!-----------------------------------------------------------
!!
!! Gauss-Lobatto-Legendre resolution
!!
!!-----------------------------------------------------------

! number of GLL points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX
  integer, parameter :: NGLLZ = NGLLX

!!-----------------------------------------------------------
!!
!! I/O output
!!
!!-----------------------------------------------------------

! input, output and main MPI I/O files
  integer, parameter :: ISTANDARD_OUTPUT = 6
  integer, parameter :: IIN = 40,IOUT = 41
! uncomment this to write messages to a text file
  integer, parameter :: IMAIN = 42
! uncomment this to write messages to the screen
! integer, parameter :: IMAIN = ISTANDARD_OUTPUT
! I/O unit for source and receiver vtk file
  integer, parameter :: IOVTK = 98
! I/O number for absorbing boundary snapshots (opened in C for speed)
! Not a Fortran I/O unit! See write_c_binary.c for more information.
  integer, parameter :: IOABS = 0
  integer, parameter :: IOABS_AC = 1
! I/O unit for plotting source time function
  integer, Parameter :: IOSTF = 71
! I/O unit for noise simulations
  integer, parameter :: IIN_NOISE = 44, IOUT_NOISE = 45
! I/O unit for SU formatted seismograms
  integer, parameter :: IOUT_SU = 46
! I/O unit for SU formatted adjoint sources
  integer, parameter :: IIN_SU1 = 47, IIN_SU2 = 48, IIN_SU3 = 49

! ignore variable name field (junk) at the beginning of each input line
  logical, parameter :: IGNORE_JUNK = .true.,DONT_IGNORE_JUNK = .false.

  ! size of the ADIOS buffer to use
  integer, parameter :: ADIOS_BUFFER_SIZE_IN_MB = 200
  character(len=*), parameter :: ADIOS_TRANSPORT_METHOD = 'MPI'

  ! ASDF parameter
  logical, parameter :: OUTPUT_PROVENANCE = .false.

!!-----------------------------------------------------------
!!
!! smoothing of the sensitivity kernels
!!
!!-----------------------------------------------------------

! using this is more precise but more costly (but will be replaced with a cheap Bessel smoothing in the near future)
! see https://github.com/geodynamics/specfem3d/issues/621
  logical, parameter :: USE_QUADRATURE_RULE_FOR_SMOOTHING = .true.

!!-----------------------------------------------------------
!!
!! gravity integral calculations
!!
!!-----------------------------------------------------------

!! DK DK for gravity integrals
  logical, parameter :: GRAVITY_INTEGRALS = .false. !!! .true.

! check for negative Jacobians in the calculation of integrals or not
! (can safely be done once to check that the mesh is OK at a given resolution, and then permanently
!  turned off in future runs because the mesh does not change)
  logical, parameter :: CHECK_FOR_NEGATIVE_JACOBIANS = .true.

! file that contains the observation grid (which is a simple ASCII file, with x y z of each point on separate lines)
  character(len=*), parameter :: OBSERVATION_GRID_FILE = './DATA/observation_grid_to_use_for_gravity.txt'

! number of points in the observation grid
  integer, parameter :: NTOTAL_OBSERVATION = 4

! a receiver at which we check and output the gravity field (for display purposes only; you can thus safely leave it unchanged)
  integer, parameter :: iobs_receiver = 1

! how often (every how many spectral elements computed) we print a timestamp to monitor the behavior of the code
  integer, parameter :: NSPEC_DISPLAY_INTERVAL = 1000

! definition of an Eotvos compared to S.I. units.
! The unit of gravity gradient is the Eotvos, which is equivalent to 1e-9 s-2 (or 1e-4 mGal/m).
! A person walking at a distance of 2 meters provides a gravity gradient signal of approximately one Eotvos.
! Mountains can create signals of several hundred Eotvos.
  double precision, parameter :: SI_UNITS_TO_EOTVOS = 1.d+9

!!-----------------------------------------------------------
!!
!! source/receiver setup
!!
!!-----------------------------------------------------------

! sources

! flag to print the details of source location
  logical, parameter :: SHOW_DETAILS_LOCATE_SOURCE = .false.

! maximum length of station and network name for receivers
  integer, parameter :: MAX_LENGTH_STATION_NAME = 32
  integer, parameter :: MAX_LENGTH_NETWORK_NAME = 8

! number of sources to be gathered by MPI_Gather
  integer, parameter :: NGATHER_SOURCES = 100

! we mimic a triangle of half duration equal to half_duration_triangle
! using a Gaussian having a very close shape, as explained in Figure 4.2
! of the manual. This source decay rate to mimic an equivalent triangle
! was found by trial and error
  double precision, parameter :: SOURCE_DECAY_MIMIC_TRIANGLE = 1.628d0

! receivers

! use this t0 as earliest starting time rather than the automatically calculated one
! (must be positive and bigger than the automatically one to be effective;
!  simulation will start at t = - t0)
  double precision, parameter :: USER_T0 = 0.0d0

! the receivers can be located inside the model
  logical, parameter :: RECEIVERS_CAN_BE_BURIED = .true.
  logical, parameter :: SOURCES_CAN_BE_BURIED = .true.

! the seismograms are normal to the surface.
! the Z record corresponds to the normal, while E and N are two tangent vectors
! that completes an orthonormal.
  logical, parameter :: EXTERNAL_MESH_RECEIVERS_NORMAL = .false.

! maximum number of sources and receivers to locate simultaneously
  integer, parameter :: NSOURCES_SUBSET_MAX = 100
  integer, parameter :: NREC_SUBSET_MAX = 200

! use brute-force search (loops over all elements) or search in kd-tree for locating target point in mesh
! (brute-force search will be slower but even more accurate)
  logical, parameter :: DO_BRUTE_FORCE_POINT_SEARCH = .true.

!!-----------------------------------------------------------
!!
!! CPML perfectly matched absorbing layers
!!
!!-----------------------------------------------------------

! power (exponent) used for the PML damping profile
  double precision, parameter :: NPOWER = 1.d0

! C-PML theoretical reflection coefficient
! (INRIA research report section 6.1:  http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf)
  double precision, parameter :: CPML_Rcoef = 0.001d0

! flags for the seven CPML regions
  integer, parameter :: CPML_X_ONLY = 1
  integer, parameter :: CPML_Y_ONLY = 2
  integer, parameter :: CPML_Z_ONLY = 3
  integer, parameter :: CPML_XY_ONLY = 4
  integer, parameter :: CPML_XZ_ONLY = 5
  integer, parameter :: CPML_YZ_ONLY = 6
  integer, parameter :: CPML_XYZ = 7

! tolerance to automatically detect where the PMLs are in the case of the internal mesher
  double precision, parameter :: SMALL_PERCENTAGE_TOLERANCE = 1.005d0

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

!! DK DK this paragraph seems to be from Zhinan or from ChangHua:
! However, enforcing explicitly potential_dot_dot_acoustic, potential_dot_acoustic, potential_acoustic
! to be zero on outer boundary of PML help to improve the accuracy of absorbing low-frequency wave components
! in case of long-time simulation.

! impose Dirichlet conditions for the potential (i.e. Neumann for displacement) on the outer edges of the C-PML layers
  logical, parameter :: SET_NEUMANN_RATHER_THAN_DIRICHLET_FOR_FLUID_PMLs = .true.

!!-----------------------------------------------------------
!!
!! For coupling with EXTERNAL CODE
!!
!!-----------------------------------------------------------

! add support for AXISEM and FK as external codes, not only DSM
  integer, parameter :: INJECTION_TECHNIQUE_IS_DSM    = 1
  integer, parameter :: INJECTION_TECHNIQUE_IS_AXISEM = 2
  integer, parameter :: INJECTION_TECHNIQUE_IS_FK     = 3

! Big storage version of coupling with DSM (will always set to .false. after the light storage version will be validated)
  logical, parameter :: old_DSM_coupling_from_Vadim = .true.

! First run to store on boundaries for using reciprocty calculate Kirchoff-Helmholtz integral
  logical, parameter :: SAVE_RUN_BOUN_FOR_KH_INTEGRAL = .false.

! for subroutine compute_vol_or_surf_integral_on_whole_domain
  integer, parameter :: Surf_or_vol_integral    = 1  !!! 1 = Surface integral, 2 = Volume, 3 = Both

! some old tests (currently unstable; do not remove them though, we might fix this one day)
  integer, parameter :: Ntime_step_dsm = 100
  integer, parameter :: IIN_veloc_dsm = 51, IIN_tract_dsm = 52
  integer, parameter :: IIN_displ_axisem = 54

! stores mesh files as cubit-exported files into directory MESH/ (for single process run)
! todo: we could put this parameter into the Mesh_Par_file
  logical, parameter :: SAVE_MESH_AS_CUBIT = .false.


!!-----------------------------------------------------------
!!
!! directory structure
!!
!!-----------------------------------------------------------

! paths for inputs and outputs files
  character(len=*), parameter :: IN_DATA_FILES = './DATA/'
  character(len=*), parameter :: MF_IN_DATA_FILES = './DATA/meshfem3D_files/'
  character(len=*), parameter :: OUTPUT_FILES_BASE = './OUTPUT_FILES/'


!!-----------------------------------------------------------
!!
!! image output
!!
!!-----------------------------------------------------------

! plots VTK cross-section planes instead of model surface
! (EXPERIMENTAL feature)
! (requires MOVIE_SURFACE set to .true. in Par_file)
  logical, parameter :: PLOT_CROSS_SECTIONS = .false.
  real(kind=CUSTOM_REAL),parameter :: CROSS_SECTION_X = 67000.0_CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: CROSS_SECTION_Y = 65500.0_CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: CROSS_SECTION_Z = -30000.0_CUSTOM_REAL

! plots PNM cross-section image
! (EXPERIMENTAL feature)
! (cross-section plane parameters can be specified in create_color_image.f90)
!! DK DK Jan 2013: here for performance and to reduce the size of the files, one day
!! DK DK Jan 2013: we should switch to using the JPEG library directly, as already implemented in SPECFEM2D
  logical, parameter :: PNM_IMAGE = .false.

! geometry tolerance parameter to calculate number of independent grid points
! sensitive to actual size of model, assumes reference sphere of radius 1
! this is an absolute value for normalized coordinates in the Earth
  double precision, parameter :: SMALLVAL_TOL = 1.d-10

!!-----------------------------------------------------------
!!
!! mesh coloring for GPUs (not needed any more, obsolete)
!!
!!-----------------------------------------------------------

! add mesh coloring for the GPU + MPI implementation
! this is needed on NVIDIA hardware up to FERMI boards, included.
! Starting on KEPLER boards you can leave it off because on KEPLER hardware
! or higher atomic reduction operations have become as fast as resorting to mesh coloring.
  logical, parameter :: USE_MESH_COLORING_GPU = .false.
  integer, parameter :: MAX_NUMBER_OF_COLORS = 1000

! enhanced coloring:
!
! using Droux algorithm
! try several times with one more color before giving up
  logical, parameter :: USE_DROUX_OPTIMIZATION = .false.
  integer, parameter :: MAX_NB_TRIES_OF_DROUX_1993 = 15
!
! postprocess the colors to balance them if Droux (1993) algorithm is not used
  logical, parameter :: BALANCE_COLORS_SIMPLE_ALGO = .false.

!!-----------------------------------------------------------
!!
!! in-situ visualization (using VTK)
!!
!!-----------------------------------------------------------

! use low-res/hi-res mesh point locations (corners only or all GLL points)
  logical, parameter :: VTK_USE_HIRES         = .false.

! show free surface points
  logical, parameter :: VTK_SHOW_FREESURFACE  = .true.

! show volumetric field
  logical, parameter :: VTK_SHOW_VOLUME       = .true.

!!-----------------------------------------------------------
!! Option to run several events within the same MPI slice when using GPU acoustic simulations
!!-----------------------------------------------------------

  integer, parameter :: NB_RUNS_ACOUSTIC_GPU = 1

!------------------------------------------------------
!----------- do not modify anything below -------------
!------------------------------------------------------

!! DK DK August 2018: on modern processors this does not happen any more,
!! DK DK August 2018: and thus no need to purposely lose accuracy to avoid underflows; thus turning it off by default
  logical, parameter :: FIX_UNDERFLOW_PROBLEM = .false. ! .true.

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: TWO_PI = 2.d0 * PI

  double precision, parameter :: DEG2RAD = PI / 180.d0
  double precision, parameter :: RAD2DEG = 180.d0 / PI
  double precision, parameter :: DEG2KM = 111.19492664455873d0
  double precision, parameter :: KM2DEG = 1.d0/(6371.0d0*pi/180.0d0)

! 3-D simulation
  integer, parameter :: NDIM = 3

! dimension of the boundaries of the slices
  integer, parameter :: NDIM2D = 2

! this to only take the corners (i.e. extract a QUAD4 or HEX8 in some routines even if we use higher-order elements)
  integer, parameter :: NGNOD_EIGHT_CORNERS = 8
  integer, parameter :: NGNOD2D_FOUR_CORNERS = 4

! number of points in each AVS or OpenDX quadrangular cell for movies
  integer, parameter :: NGNOD2D_FOUR_CORNERS_AVS_DX = NGNOD2D_FOUR_CORNERS

! number of points per surface element
  integer, parameter :: NGLLSQUARE = NGLLX * NGLLY

! for optimized routines by Deville et al. (2002)
  integer, parameter :: m1 = NGLLX, m2 = NGLLX * NGLLY
  integer, parameter :: NGLLCUBE = NGLLX * NGLLY * NGLLZ

! mid-points inside a GLL element
  integer, parameter :: MIDX = (NGLLX+1)/2
  integer, parameter :: MIDY = (NGLLY+1)/2
  integer, parameter :: MIDZ = (NGLLZ+1)/2

! a few useful constants
  double precision, parameter :: ZERO = 0.d0, ONE = 1.d0, TWO = 2.d0, HALF = 0.5d0

  real(kind=CUSTOM_REAL), parameter :: &
    ONE_THIRD   = 1._CUSTOM_REAL/3._CUSTOM_REAL, &
    FOUR_THIRDS = 4._CUSTOM_REAL/3._CUSTOM_REAL

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9

! tiny real value declared independently of the machine
  real(kind=CUSTOM_REAL), parameter :: TINYVAL_SNGL = 1.e-25_CUSTOM_REAL

! number of standard linear solids in parallel for attenuation
  integer, parameter :: N_SLS = 3

! computation of standard linear solids
! ATTENUATION_COMP_RESOLUTION: Number of Digits after decimal
! ATTENUATION_COMP_MAXIMUM:    Maximum Q Value
  integer, parameter :: ATTENUATION_COMP_RESOLUTION = 1
  integer, parameter :: ATTENUATION_COMP_MAXIMUM    = 9000

! define flag for regions for anisotropy
  integer, parameter :: IANISOTROPY_MODEL1 = 1
  integer, parameter :: IANISOTROPY_MODEL2 = 2

! smallest real number on the Pentium and the SGI =  1.1754944E-38
! largest real number on the Pentium and the SGI  =  3.4028235E+38
! small negligible initial value to avoid very slow underflow trapping
! but not too small to avoid trapping on velocity and acceleration in Newmark
  real(kind=CUSTOM_REAL), parameter :: VERYSMALLVAL = 1.E-24_CUSTOM_REAL

! displacement threshold above which we consider the code became unstable
  real(kind=CUSTOM_REAL), parameter :: STABILITY_THRESHOLD = 1.E+25_CUSTOM_REAL

! geometrical tolerance for boundary detection
  double precision, parameter :: SMALLVAL = 0.00001d0

! do not use tags for MPI messages, use dummy tag instead
  integer, parameter :: itag = 0,itag2 = 0

! for the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0

! number of lines per source in CMTSOLUTION file
  integer, parameter :: NLINES_PER_CMTSOLUTION_SOURCE = 13
  integer, parameter :: NLINES_PER_FORCESOLUTION_SOURCE = 10

! number of iterations to solve the system for xi and eta
! setting it to 5 instead of 4 ensures that the result obtained is not compiler dependent
! (when using 4 only some small discrepancies were observed)
  integer, parameter :: NUM_ITER = 5

! flag to exclude elements that are too far from target in topography detection
  logical, parameter :: USE_DISTANCE_CRITERION_TOPO = .true.

! flag for projection from latitude/longitude to UTM, and back
  integer, parameter :: ILONGLAT2UTM = 0, IUTM2LONGLAT = 1

! for total energy calculation
  integer, parameter :: IOUT_ENERGY = 937  ! file number for the energy file

!!-----------------------------------------------------------
!!
!! Topography file
!!
!!-----------------------------------------------------------
! specifies topography file to determine the (exact) elevation of each GLL point at the free surface
! in case TOPOGRAPHY = .true. (used for ocean load approximation)

!! Socal1D example
!! size of topography and bathymetry file for EXAMPLE/meshfem3D_examples/socal1D/DATA/meshfem3D_files/interface4.dat
  integer, parameter :: NX_TOPO_FILE = 2,NY_TOPO_FILE = 2
  double precision, parameter :: ORIG_LAT_TOPO = 32.d0
  double precision, parameter :: ORIG_LONG_TOPO = -121.d0
  double precision, parameter :: DEGREES_PER_CELL_TOPO = 6.d0
  character(len=*), parameter :: TOPO_FILE = 'DATA/meshfem3D_files/interface4.dat'

!! Many interfaces example
!! size of topography and bathymetry file for EXAMPLE/meshfem3D_examples/socal1D/DATA/meshfem3D_files/interface4.dat
!  integer, parameter :: NX_TOPO_FILE = 221,NY_TOPO_FILE = 200
!  double precision, parameter :: ORIG_LAT_TOPO = 0.d0
!  double precision, parameter :: ORIG_LONG_TOPO = 0.d0
!  double precision, parameter :: DEGREES_PER_CELL_TOPO = 150.d0
!  character(len=*), parameter :: TOPO_FILE = 'DATA/meshfem3D_files/example_topo.dat'

!! LA Southern California
!! size of topography and bathymetry file for Southern California
!  integer, parameter :: NX_TOPO_FILE = 1401,NY_TOPO_FILE = 1001
!  double precision, parameter :: ORIG_LAT_TOPO = 32.d0
!  double precision, parameter :: ORIG_LONG_TOPO = -121.d0
!  double precision, parameter :: DEGREES_PER_CELL_TOPO = 5.d0 / 1000.d0
!  character(len=*), parameter :: TOPO_FILE = 'DATA/la_topography/topo_bathy_final.dat'

!! Piero's model
!! size of topography and bathymetry file for Piero Basini's model
!  integer, parameter :: NX_TOPO_FILE = 787, NY_TOPO_FILE = 793
!  double precision, parameter :: ORIG_LAT_TOPO = -102352.d0
!  double precision, parameter :: ORIG_LONG_TOPO = 729806.d0
!! for Piero Basini's model this is the resolution in meters of the topo file
!  double precision, parameter :: DEGREES_PER_CELL_TOPO = 250.d0
!  character(len=*), parameter :: TOPO_FILE = 'DATA/piero_model/dem_EV_UTM_regular_250_reordered.dat'

!!-----------------------------------------------------------
!!
!! APPROXIMATE_OCEAN_LOAD load approximation
!!
!!-----------------------------------------------------------
! minimum thickness in meters to include the effect of the oceans
! to avoid taking into account spurious oscillations in topography model
  double precision, parameter :: MINIMUM_THICKNESS_3D_OCEANS = 10.d0
! density of sea water
  real(kind=CUSTOM_REAL), parameter :: RHO_APPROXIMATE_OCEAN_LOAD = 1020.0_CUSTOM_REAL

!!-----------------------------------------------------------
!!
!! GRAVITY
!!
!!-----------------------------------------------------------
! gravitational constant
  double precision, parameter :: GRAV = 6.6723d-11
! number of layers in PREM
  integer, parameter :: NR = 640
! R_EARTH is the radius of the bottom of the oceans (radius of Earth in m)
  double precision, parameter :: R_EARTH = 6371000.d0
! same radius in km
  double precision, parameter :: R_EARTH_KM = R_EARTH / 1000.d0
! radius of the Earth for gravity calculation
  double precision, parameter :: R_EARTH_GRAVITY = 6371000.d0
! radius of the ocean floor for gravity calculation
  double precision, parameter :: ROCEAN_GRAVITY = 6368000.d0
! average density in the full Earth to normalize equation
  double precision, parameter :: RHOAV = 5514.3d0

!!-----------------------------------------------------------
!!
!! DOMAINS
!!
!!-----------------------------------------------------------
! material domain ids
  integer, parameter :: IDOMAIN_ACOUSTIC    = 1
  integer, parameter :: IDOMAIN_ELASTIC     = 2
  integer, parameter :: IDOMAIN_POROELASTIC = 3

! model ids
  integer, parameter :: IMODEL_DEFAULT          = 1
  integer, parameter :: IMODEL_1D_PREM          = 2
  integer, parameter :: IMODEL_1D_SOCAL         = 3
  integer, parameter :: IMODEL_1D_CASCADIA      = 4
  integer, parameter :: IMODEL_TOMO             = 5
  integer, parameter :: IMODEL_USER_EXTERNAL    = 6
  integer, parameter :: IMODEL_GLL              = 7
  integer, parameter :: IMODEL_SALTON_TROUGH    = 8
  integer, parameter :: IMODEL_1D_PREM_PB       = 9
  integer, parameter :: IMODEL_IPATI            = 10
  integer, parameter :: IMODEL_IPATI_WATER      = 11
  integer, parameter :: IMODEL_SEP              = 12
  integer, parameter :: IMODEL_COUPLED          = 13

!!-----------------------------------------------------------
!!
!! time stamp information
!!
!!-----------------------------------------------------------
! print date and time estimate of end of run in another country,
! in addition to local time.
! For instance: the code runs at Caltech in California but the person
! running the code is connected remotely from France, which has 9 hours more.
! The time difference with that remote location can be positive or negative
  logical, parameter :: ADD_TIME_ESTIMATE_ELSEWHERE = .false.
  integer, parameter :: HOURS_TIME_DIFFERENCE = +9
  integer, parameter :: MINUTES_TIME_DIFFERENCE = +0

!!-----------------------------------------------------------
!!
!! for LDDRK high-order time scheme
!!
!!-----------------------------------------------------------
  integer, parameter :: NSTAGE = 6

  real(kind=CUSTOM_REAL), dimension(NSTAGE), parameter :: ALPHA_LDDRK = &
    (/0.0_CUSTOM_REAL,-0.737101392796_CUSTOM_REAL, -1.634740794341_CUSTOM_REAL, &
      -0.744739003780_CUSTOM_REAL,-1.469897351522_CUSTOM_REAL,-2.813971388035_CUSTOM_REAL/)

  real(kind=CUSTOM_REAL), dimension(NSTAGE), parameter :: BETA_LDDRK = &
    (/0.032918605146_CUSTOM_REAL,0.823256998200_CUSTOM_REAL,0.381530948900_CUSTOM_REAL, &
      0.200092213184_CUSTOM_REAL,1.718581042715_CUSTOM_REAL,0.27_CUSTOM_REAL/)

  real(kind=CUSTOM_REAL), dimension(NSTAGE), parameter :: C_LDDRK = &
    (/0.0_CUSTOM_REAL,0.032918605146_CUSTOM_REAL,0.249351723343_CUSTOM_REAL, &
      0.466911705055_CUSTOM_REAL,0.582030414044_CUSTOM_REAL,0.847252983783_CUSTOM_REAL/)

!!----------------------
!!
!!  VM VM test dipole source in fluid regionn
!!
!!  if .true. this will change the way the force solution is set in fluid region, trying to put a dipole in fluids
!!  if .false. the behaviour of code remains unchanged
!!
!!-----------
  logical, parameter  ::    DIPOLE_SOURCE_IN_FLUID = .false.
