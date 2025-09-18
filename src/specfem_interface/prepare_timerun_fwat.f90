!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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


  subroutine prepare_timerun_fwat(ievt)

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use fk_coupling, only: read_fk_coupling_file, check_fk_files, &
                        couple_with_injection_prepare_boundary_fwat
  use input_params, only: fpar => fwat_par_global
  use logger, only: log

  !! solving wavefield discontinuity problem with non-split-node scheme
  use wavefield_discontinuity_solver, only: &
               prepare_timerun_wavefield_discontinuity

  implicit none

  ! local parameters
  double precision :: tCPU,tstart
  double precision, external :: wtime
  integer :: ievt

  ! checks if anything to do
  if (.not. IO_compute_task) return

  ! synchonizes
  call synchronize_all()

  ! get MPI starting time
  tstart = wtime()

  ! user output infos
  call prepare_timerun_user_output()

  ! sets up mass matrices
  call prepare_timerun_mass_matrices()

  ! sets up time increments
  call prepare_timerun_constants()

  ! initializes arrays
  call prepare_wavefields()

  ! initializes fault rupture arrays
  call prepare_timerun_faults()

  ! prepares attenuation arrays
  call prepare_attenuation()

  ! prepares gravity arrays
  call prepare_gravity()

  ! ZN I do not use if (USE_LDDRK) call prepare_timerun_lddrk()
  ! ZN in order to avoid the error of using unallocated arrays later on in the code,
  ! ZN since R_**_lddrk are arguments in subroutine compute_forces_viscoelastic
  call prepare_timerun_lddrk()

  ! prepares C-PML arrays
  if (PML_CONDITIONS) call prepare_timerun_pml()

  ! Stacey boundaries
  call prepare_timerun_stacey()

  ! wavefield discontinuity
  call prepare_timerun_wavefield_discontinuity()

  ! prepares ADJOINT simulations
  call prepare_timerun_adjoint()

  ! prepares noise simulations
  call prepare_noise()

  ! prepares coupling with injection boundary
  ! call couple_with_injection_prepare_boundary()
  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE==1) then
    if (check_fk_files(fpar%acqui%evtid_names(ievt))) then
      call log%write('Read FK wavefield for event '//trim(fpar%acqui%evtid_names(ievt)), .true.)
      call read_fk_coupling_file(fpar%acqui%evtid_names(ievt))
    else
      call log%write('Calculating FK wavefield for event '//trim(fpar%acqui%evtid_names(ievt)), .true.)
      if (.not. fpar%sim%SAVE_FK) then
        call couple_with_injection_prepare_boundary()
      else
        call couple_with_injection_prepare_boundary_fwat(fpar%acqui%evtid_names(ievt))
      endif
    endif
    call log%write('Finished FK wavefield preparation', .true.)
  endif

  ! prepares GPU arrays
  if (GPU_MODE) call prepare_GPU()

  ! optimizes array memory layout for better performance
  call prepare_optimized_arrays()

  ! synchronize all the processes
  call synchronize_all()

  ! elapsed time since beginning of preparation
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for preparing timerun in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) '************'
    write(IMAIN,*) ' time loop'
    write(IMAIN,*) '************'
    if (USE_LDDRK) then
      write(IMAIN,*) '              scheme:         LDDRK with',NSTAGE_TIME_SCHEME,'stages'
    else
      write(IMAIN,*) '              scheme:         Newmark'
    endif
    write(IMAIN,*)
    write(IMAIN,*) '           time step: ',sngl(DT),' s'
    write(IMAIN,*) 'number of time steps: ',NSTEP
    write(IMAIN,*) 'total simulated time: ',sngl(NSTEP*DT),' seconds'
    write(IMAIN,*) 'start time:',sngl(-t0),' seconds'
    write(IMAIN,*)

    ! flushes file buffer for main output file (IMAIN)
    call flush_IMAIN()

    !#TODO: total time estimation based on number of elements and number of time steps before time loop starts?
    !  average time per element per time step:
    !     elastic elements    ~ dt = 1.17e-05 s (intel xeon 2.6GHz, stand 2013)
    !                              = 3.18e-07 s (Kepler K20x, stand 2013)
    !
    !  total time per time step:
    !     T_total = dt * nspec_total
    !
    !  total time using nproc processes (slices) for NSTEP time steps:
    !     T_simulation = T_total * NSTEP / nproc

  endif

  ! synchronize all the processes
  call synchronize_all()

  end subroutine prepare_timerun_fwat

!