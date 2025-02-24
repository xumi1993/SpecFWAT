program fwat_fwd_measure_adj

use config
use fwat_mpi
use specfem_par, only: myrank
use imput_params, fpar => fwat_par_global

implicit none

call init_mpi()
call init_mpi_tomo()

call fpar%read(FWAT_PAR_FILE)

if (worldrank == 0) then
  print *, tele_par%RCOMPS
endif


end program fwat_fwd_measure_adj