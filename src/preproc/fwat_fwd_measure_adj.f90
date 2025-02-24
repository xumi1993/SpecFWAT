program fwat_fwd_measure_adj
use config
use fwat_mpi
use imput_params, fpar => fwat_par_global

implicit none

call init_mpi()
call init_mpi_fwat()

call fpar%read(FWAT_PAR_FILE)

! if (worldrank == 0) then
!   print *, tele_par%RCOMPS
! endif

call finalize_mpi()

end program fwat_fwd_measure_adj