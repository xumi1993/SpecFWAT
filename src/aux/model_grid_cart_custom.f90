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

  program model_grid_cart_custom

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
  use input_params, only: fpar => fwat_par_global

  implicit none

  ! integer, parameter :: NARGS = 6

  ! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  double precision,dimension(:,:,:,:),allocatable :: data_dp
  ! real array for data
  real,dimension(:,:,:,:),allocatable :: data_sp
  integer :: i, j, k, ier
  character(len=MAX_STRING_LEN) :: arg(9), indir, outdir, replace_zero
  character(len=MAX_STRING_LEN) :: prname, prname_lp, data_filename, fname, grid_file, model_name, hstr
  character(len=MAX_STRING_LEN*2) :: local_data_file
  logical :: BROADCAST_AFTER_READ
  integer :: ios
  integer :: sizeprocs, NSPEC_IRREGULAR
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            x_min_glob, x_max_glob, &
                            y_min_glob, y_max_glob, &
                            z_min_glob, z_max_glob, &
                            x_min, y_min, z_min, &
                            x_max, y_max, z_max

  type(profd)  :: projection_fd
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: model_on_FD_grid
  double precision, dimension(:), allocatable :: xfd, yfd, zfd
  double precision :: ox, oy
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
  call read_parameter_file(BROADCAST_AFTER_READ)
  call fpar%read('DATA/fwat_params.yml')

  ! parse command line arguments
  if (command_argument_count() ==0) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xcombine_vol_data_on_regular_grid &
               data_filename input_dir output_dir model_name [is_replace_zero]'
      stop 'Please check command line arguments'
    endif
  endif

  ! reads in arguments
  ! do i = 1, NARGS
  !   call get_command_argument(i,arg(i))
  ! enddo

  call get_command_argument(1,data_filename)
  call get_command_argument(2,indir)
  call get_command_argument(3,outdir)
  call get_command_argument(4,model_name)
  if (command_argument_count() == 5) then
    call get_command_argument(5,replace_zero)
  elseif (command_argument_count() == 4) then
    replace_zero = 'true'
  endif


    hx_fd_proj = fpar%grid%regular_grid_interval(1)
    hy_fd_proj = fpar%grid%regular_grid_interval(2)
    hz_fd_proj = fpar%grid%regular_grid_interval(3)
    ox_fd_proj = fpar%grid%regular_grid_min_coord(1)
    oy_fd_proj = fpar%grid%regular_grid_min_coord(2)
    oz_fd_proj = fpar%grid%regular_grid_min_coord(3)
    nx_fd_proj = fpar%grid%regular_grid_size(1)
    ny_fd_proj = fpar%grid%regular_grid_size(2)
    nz_fd_proj = fpar%grid%regular_grid_size(3)

  if (myrank==0) print *, 'Reading GLL mesh...'
  ! Get dimensions of current model, stored in proc******_external_mesh.bin
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'
  open(unit=27,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  read(27) NSPEC_AB
  read(27) NGLOB_AB
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1102')
  if (ier /= 0) stop 'error allocating array ibool'
  allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1103')
  if (ier /= 0) stop 'error allocating array xstore etc.'

  read(27) NSPEC_IRREGULAR
  read(27) ibool
  read(27) xstore
  read(27) ystore
  read(27) zstore
  close(27)

  ! Get regular grid properties, and precomputes all interpolation coefficients
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
 
  if (myrank==0) print *, 'Initialize regular grid...'
  call compute_interpolation_coeff_FD_SEM(projection_fd, myrank)
  ! allocate(model_on_FD_grid(projection_fd%nx, projection_fd%ny, projection_fd%nz),stat=ier)
  allocate(xfd(projection_fd%nx))
  allocate(yfd(projection_fd%ny))
  allocate(zfd(projection_fd%nz))
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1104')

  if (myrank == 0) then
    print *, 'Grid size are : ',projection_fd%nx, projection_fd%ny, projection_fd%nz
  endif

  ! Get data to project
  allocate(data_sp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1105')
  if (ier /= 0) stop 'error allocating single precision data array'
  if (CUSTOM_REAL == SIZE_DOUBLE) then
    allocate(data_dp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1106')
    if (ier /= 0) stop 'error allocating double precision data array'
  endif

  ! data file
  write(prname,'(a,i6.6,a)') trim(indir)//'/proc',myrank,'_'
  local_data_file = trim(prname) // trim(data_filename) // '.bin'
  open(unit = 28,file = trim(local_data_file),status='old', &
        action='read',form ='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening ',trim(local_data_file)
    stop
  endif

  ! Read either SP or DP floating point numbers.
  if (CUSTOM_REAL == SIZE_DOUBLE) then
    read(28) data_dp
  else
    read(28) data_sp
  endif
  close(28)

  ! uses conversion to real values
  if (CUSTOM_REAL == SIZE_DOUBLE) then
    data_sp(:,:,:,:) = sngl(data_dp(:,:,:,:))
    deallocate(data_dp)
  endif
  call synchronize_all()
  if (myrank==0) print *, 'Finish reading data'

  ! put data from SEM mesh to a regular grid
  call Project_model_SEM2FD_grid(data_sp, model_on_FD_grid, projection_fd, myrank)
  call synchronize_all()
  if (myrank==0) print *, 'Finish grid data'
  ! Write output on a Fortran binary file
  if (myrank == 0) then
    ox = dble(projection_fd%ox)
    oy = dble(projection_fd%oy)

    do i=1,projection_fd%nx
      xfd(i) = ox + (i-1) * projection_fd%hx
    enddo
    do i=1,projection_fd%ny
      yfd(i) = oy + (i-1) * projection_fd%hy
    enddo
    do i=1,projection_fd%nz
      zfd(i) = (projection_fd%oz + (i-1) * projection_fd%hz)
    enddo

    if (replace_zero=='true') then
      do i=1,projection_fd%nx
        do j=1,projection_fd%ny
          do k=1,projection_fd%nz
            if (model_on_FD_grid(i,j,k) == 0) then
              if (zfd(k)>=0) then
                call find_nearestZ_nonzero(model_on_FD_grid,i,j,k,&
                      projection_fd%nx,projection_fd%ny,projection_fd%nz)
              else
                call find_nearestXY_nonzero(model_on_FD_grid,i,j,k,&
                      projection_fd%nx,projection_fd%ny,projection_fd%nz)
              endif
            endif
          enddo
        enddo
      enddo
    endif

    fname =  trim(outdir)//'/'//trim(data_filename)//'_'//trim(model_name)//'.h5'
    call h5file%open(fname, status='new', action='write')
    call h5file%add('/x', xfd)
    call h5file%add('/y', yfd)
    call h5file%add('/z', zfd)
    call h5file%add('/'//trim(data_filename), transpose_3(model_on_FD_grid))
    print *, 'Done writing ', trim(fname)
  endif

  call finalize_mpi()

end program model_grid_cart_custom
