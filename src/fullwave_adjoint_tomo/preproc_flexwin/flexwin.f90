!----------------------------------------------------
subroutine flexwin_fwat()

use seismo_variables


! read the parameter file
!call read_parameter_file_flexwin()

!write(*,*) 'Number of seismograms to measure : '
!read(*,*) n_seis
!write(*,*) n_seis

!write(*,*) 'For each seismogram: '
!do i = 1, 1!n_seis
!  write(*,*) 'Observed seismogram'
!  !read(*,'(a)') obs_name(i)
!  write(*,'(a)') trim(obs_name)
!  write(*,*) 'Synthetic seismogram'
!  !read(*,'(a)') syn_name(i)
!  write(*,'(a)') trim(syn_name)
!  write(*,*) 'Output basename'
!  !read(*,'(a)') basename(i)
!  write(*,'(a)') trim(basename)
!
!  if (DEBUG) write(*,*) 'DEBUG : reading sac files'
!  call read_sac_files(trim(syn_name),trim(obs_name),ier)
!  if( ier /= 0 ) then
!    write(*,*) 'error files read in:',trim(syn_name)
!    write(*,*) 'skipping files'
!    cycle
!  endif

  if (DEBUG) write(*,*) 'DEBUG : selecting windows'
  call select_windows_stalta2()


  !if(MAKE_SEISMO_PLOTS) then
  !   if (DEBUG) write(*,*) 'DEBUG : writing output seismos'
  !   call write_seismos_gmt(basename)
  !endif

  !if(MAKE_WINDOW_FILES) then
  !   if (DEBUG) write(*,*) 'DEBUG : writing mt input'
  !   call write_mt_input_2(basename,obs_name,syn_name)
  !endif

!enddo

end subroutine flexwin_fwat
