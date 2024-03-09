!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.0
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
!

program xflexwin
use seismo_variables
implicit none

character(len=256) :: syn_name,obs_name,basename
integer :: n_seis, i
integer :: ier

! read the parameter file
call read_parameter_file_flexwin()

write(*,*) 'Number of seismograms to measure : '
read(*,*) n_seis
write(*,*) n_seis

write(*,*) 'For each seismogram: '
do i = 1, n_seis
  write(*,*) 'Observed seismogram'
  read(*,'(a)') obs_name
  write(*,'(a)') trim(obs_name)
  write(*,*) 'Synthetic seismogram'
  read(*,'(a)') syn_name
  write(*,'(a)') trim(syn_name)
  write(*,*) 'Output basename'
  read(*,'(a)') basename
  write(*,'(a)') trim(basename)

  if (DEBUG) write(*,*) 'DEBUG : reading sac files'
  call read_sac_files(trim(syn_name),trim(obs_name),ier)
  if( ier /= 0 ) then
    write(*,*) 'error files read in:',trim(syn_name)
    write(*,*) 'skipping files'
    cycle
  endif

  if (DEBUG) write(*,*) 'DEBUG : selecting windows'
  call select_windows_stalta2()

  if(MAKE_SEISMO_PLOTS) then
     if (DEBUG) write(*,*) 'DEBUG : writing output seismos'
     call write_seismos_gmt(basename)
  endif

  if(MAKE_WINDOW_FILES) then
     if (DEBUG) write(*,*) 'DEBUG : writing mt input'
     call write_mt_input_2(basename,obs_name,syn_name)
  endif

enddo



  
end program xflexwin


  subroutine usage()

  print *, 'Usage: xflexwin datafile synfile basename'
  stop

  end subroutine usage
