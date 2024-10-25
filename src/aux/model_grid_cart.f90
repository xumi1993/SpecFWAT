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

  program model_grid_cart

! combines the database files on several slices, project it on a regular grid
! and saves it in a binary file that can be read with Paraview in RAW format
!
! works for external, unregular meshes
! The region of regular grid should totally contain the SEM domain

  use constants
  use shared_parameters
  use specfem_par, only: xigll, yigll, zigll, wxgll, wygll, wzgll
  use hdf5_interface
  use utils
  use projection_on_FD_grid_fwat
  use shared_input_parameters, only: SUPPRESS_UTM_PROJECTION
  use MeshGrid

  implicit none

  ! integer, parameter :: NARGS = 6

  integer :: i, j, k, ier
  character(len=MAX_STRING_LEN) :: arg(9), indir, outdir, replace_zero
  character(len=MAX_STRING_LEN) :: data_filename, fname, model_name, hstr
  character(len=MAX_STRING_LEN*2) :: local_data_file
  logical :: BROADCAST_AFTER_READ
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            x_min_glob, x_max_glob, &
                            y_min_glob, y_max_glob, &
                            z_min_glob, z_max_glob, &
                            x_min, y_min, z_min, &
                            x_max, y_max, z_max
  type(ReglGrid) :: rg
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: model_on_FD_grid
  real(kind=CUSTOM_REAL) :: dx, dy, dz
  type(hdf5_file) :: h5file

  ! MPI initialization
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then
    print *
    print *,'Projecting volumetric data on a regular grid'
    print *
  endif

  ! needs local_path for mesh files
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! parse command line arguments
  if (command_argument_count() ==0) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xcombine_vol_data_on_regular_grid dx dy dz &
               data_filename input_dir output_dir model_name is_replace_zero'
      stop 'Please check command line arguments'
    endif
  endif

  ! reads in arguments
  ! do i = 1, NARGS
  !   call get_command_argument(i,arg(i))
  ! enddo

  call get_command_argument(1,arg(1))
  call get_command_argument(2,arg(2))
  call get_command_argument(3,arg(3))
  call get_command_argument(4,data_filename)
  call get_command_argument(5,indir)
  call get_command_argument(6,outdir)
  call get_command_argument(7,model_name)
  if (command_argument_count() == 8) then
    call get_command_argument(8,replace_zero)
  elseif (command_argument_count() == 7) then
    replace_zero = 'true'
  endif

  read(arg(1), *) dx
  read(arg(2), *) dy
  read(arg(3), *) dz

  call rg%read_mesh()

  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
         x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
         elemsize_min_glob,elemsize_max_glob, &
         distance_min_glob,distance_max_glob)
  call min_all_all_cr(x_min_glob, x_min)
  call min_all_all_cr(y_min_glob, y_min)
  call min_all_all_cr(z_min_glob, z_min)
  call max_all_all_cr(x_max_glob, x_max)
  call max_all_all_cr(y_max_glob, y_max)
  call max_all_all_cr(z_max_glob, z_max)
  ! z_max = 0._CUSTOM_REAL
  if (myrank == 0) then
    print *, 'x_min_glob = ', x_min
    print *, 'x_max_glob = ', x_max
    print *, 'y_min_glob = ', y_min
    print *, 'y_max_glob = ', y_max
    print *, 'z_min_glob = ', z_min
    print *, 'z_max_glob = ', z_max
  endif

  call rg%init(x_min, x_max, y_min, y_max, z_min, z_max, dx, dy, dz)

  call rg%griddata(indir, data_filename, model_on_FD_grid)

  if (myrank==0) print *, 'Finish grid data'
  ! Write output on a Fortran binary file
  if (myrank == 0) then

    if (replace_zero=='true') then
      do i=1, rg%nx
        do j=1, rg%ny
          do k=1, rg%nz
            if (model_on_FD_grid(i,j,k) == 0) then
              if (rg%zfd(k)>0) then
                call find_nearestZ_nonzero(model_on_FD_grid,i,j,k,rg%nx,rg%ny,rg%nz) 
              else
                call find_nearestXY_nonzero(model_on_FD_grid,i,j,k,rg%nx,rg%ny,rg%nz)
              endif
            endif
          enddo
        enddo
      enddo
    endif

    fname =  trim(outdir)//'/'//trim(data_filename)//'_'//trim(model_name)//'.h5'
    call h5file%open(fname, status='new', action='write')
    call h5file%add('/x', rg%xfd)
    call h5file%add('/y', rg%yfd)
    call h5file%add('/z', rg%zfd)
    call h5file%add('/'//trim(data_filename), transpose_3(model_on_FD_grid))
    print *, 'Done writing ', trim(fname)
  endif

  call finalize_mpi()

end program model_grid_cart
