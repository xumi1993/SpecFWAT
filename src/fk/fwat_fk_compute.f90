! Compute one selected event or all events using the shared FK preparation flow.
program fwat_fk_compute
  use fwat_mpi, only: init_mpi_single_group, init_mpi_fwat
  use argparse, only: parse_args_fk
  use prepare_fk, only: PrepareFK

  implicit none

  type(PrepareFK) :: fk
  integer :: ievt

  call init_mpi_single_group()
  call init_mpi_fwat()
  call parse_args_fk()

  ! Read FWAT parameters and the mesh once; -e limits the event range to one entry.
  call fk%init()
  do ievt = fk%first_event, fk%last_event
    call fk%prepare_for_event(ievt)
  enddo

  call fk%finalize()
  call finalize_mpi()
end program fwat_fk_compute
