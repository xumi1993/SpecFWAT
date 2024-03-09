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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine prepare_timerun_fwat(old_local_path)

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use fault_solver_dynamic, only: BC_DYNFLT_init
  use fault_solver_kinematic, only: BC_KINFLT_init

  implicit none

  ! local parameters
  double precision :: tCPU,tstart
  double precision, external :: wtime
  character(len=MAX_STRING_LEN) :: old_local_path, local_path_copy

  ! get MPI starting time
  tstart = wtime()

  ! user output infos
  call prepare_timerun_user_output()

  ! sets up mass matrices
  call prepare_timerun_mass_matrices()

  ! initializes arrays
  call prepare_wavefields()

  ! Loading kinematic and dynamic fault solvers.
  call BC_DYNFLT_init(prname,DT,myrank)

  call BC_KINFLT_init(prname,DT,myrank)

  ! sets up time increments
  call prepare_timerun_constants()

  ! prepares attenuation arrays
  local_path_copy = LOCAL_PATH
  LOCAL_PATH = old_local_path
  call prepare_attenuation()
  LOCAL_PATH = local_path_copy
  call create_name_database(prname,myrank,LOCAL_PATH) 

  ! prepares gravity arrays
  call prepare_gravity()

  ! ZN I do not use if (USE_LDDRK) call prepare_timerun_lddrk()
  ! ZN in order to avoid the error of using unallocated arrays later on in the code,
  ! ZN since R_**_lddrk are arguments in subroutine compute_forces_viscoelastic
  call prepare_timerun_lddrk()

  ! prepares C-PML arrays
  if (PML_CONDITIONS) call prepare_timerun_pml()

  ! prepares ADJOINT simulations
  call prepare_timerun_adjoint()

  ! prepares noise simulations
  call prepare_noise()

  ! prepares GPU arrays
  if (GPU_MODE) call prepare_GPU()

#ifdef USE_OPENMP
  ! prepares arrays for OpenMP
  call prepare_timerun_OpenMP()
#endif

  ! prepars coupling with injection boundary
  call couple_with_injection_prepare_boundary()

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

    !daniel debug: total time estimation
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


