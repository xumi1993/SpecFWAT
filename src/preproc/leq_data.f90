module leq_data
  use config
  use misfit_mod
  use adjoint_source, only: calculate_adjoint_source
  use adj_config, only: AdjointMeasurement
  use noise_data, only: NoiseData
  use input_params, fpar => fwat_par_global
  use shared_parameters, only: SUPPRESS_UTM_PROJECTION
  use specfem_par, only: T0, DT, NSTEP,OUTPUT_FILES
  use fwat_mpi
  use utils, only: zeros_dp, zeros
  use distaz_lib
  use common_lib, only: get_band_name, get_icomp_syn, rotate_ZRT_to_ZNE, mkdir
  use signal, only: interpolate_func_dp, dif1, detrend, demean, bandpass_dp
  use sacio
  use logger, only: log
  implicit none

  character(len=MAX_STRING_LEN), private :: msg
  integer, private :: ier
 
  type, extends(NoiseData) :: LEQData
  contains
    ! procedure :: semd2sac, preprocess, finalize
    ! procedure, private :: calc_distaz, measure_adj, write_in_preocess
  end type LEQData

contains


end module leq_data