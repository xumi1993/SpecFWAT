!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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
! Revised by Kai Wang, 2018 Dec 19
! Add two arguments: xsum_kernels ekernel_dir_list output_dir 
!==========================================================
! sum_kernels
!
! this program can be used for event kernel summation,
! where it sums up transverse isotropic kernel files:
!
!   - proc***_reg1_bulk_c_kernel.bin
!   - proc***_reg1_bulk_betav_kernel.bin
!   - proc***_reg1_bulk_betah_kernel.bin
!   - proc***_reg1_eta_kernel.bin
!
! input file: kernels_list.txt
!   lists all event kernel directories which should be summed together
!
! input directory:  INPUT_KERNELS/
!    contains links to all event kernel directories (listed in "kernels_list.txt")
!
! output directory: OUTPUT_SUM/
!    the resulting kernel files will be stored in this directory
!
!
! DEPRECATION WARNING: Eventually, all of the following routines, or at lesast
! some the subroutines, will be merged with src/tomography/xcombine_sem
!

program sum_kernels

  use tomography_par, only: MAX_STRING_LEN,MAX_KERNEL_PATHS,IIN, &
    myrank,sizeprocs,NGLOB,NSPEC,USE_ALPHA_BETA_RHO,USE_ALPHA_BETA_RHO_TISO,USE_ISO_KERNELS

  use shared_parameters

  implicit none

  character(len=MAX_STRING_LEN) :: kernel_list(MAX_KERNEL_PATHS)
  character(len=MAX_STRING_LEN) :: sline, kernel_name,prname_lp
  character(len=MAX_STRING_LEN) :: output_dir,ekernel_dir_list
  integer :: nker
  integer :: ier

  logical :: BROADCAST_AFTER_READ

  ! ============ program starts here =====================

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  call get_command_argument(1,ekernel_dir_list)
  call get_command_argument(2,output_dir)

  if (trim(ekernel_dir_list) == '' &
     .or. trim(output_dir) == '' ) then
     call exit_mpi(myrank,'USAGE: xsum_kernels ekernel_dir_list output_dir')
  endif


  if (myrank == 0) then
     write(*,*) 'SUM EVENT KERNELS TO GET MISFIT KENRELS'
     write(*,*) 'INPUT EVENT DIRECTORY:',ekernel_dir_list
     write(*,*) 'OUTPUT DIRECTORY:',output_dir
  endif
  call synchronize_all()

  ! reads in event list
  nker=0
  open(unit = IIN, file = trim(ekernel_dir_list), status = 'old',iostat = ier)
  if (ier /= 0) then
     print *,'Error opening ',trim(ekernel_dir_list),myrank
     stop 1
  endif
  do while (1 == 1)
     read(IIN,'(a)',iostat=ier) sline
     if (ier /= 0) exit
     nker = nker+1
     if (nker > MAX_KERNEL_PATHS) stop 'Error number of kernels exceeds MAX_KERNEL_PATHS'
     kernel_list(nker) = sline
  enddo
  close(IIN)
  if (myrank == 0) then
    write(*,*) '  ',nker,' events'
    write(*,*)
  endif

  ! needs local_path for mesh files
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! checks if number of MPI process as specified
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *
      print *,'Error: run xsum_kernels with the same number of MPI processes '
      print *,'       as specified in Par_file by NPROC when slices were created'
      print *
      print *,'for example: mpirun -np ',NPROC,' ./xsum_kernels ...'
      print *
    endif
    call synchronize_all()
    stop 'Error total number of slices'
  endif
  call synchronize_all()

  ! reads mesh file
  !
  ! needs to get array dimensions

  ! opens external mesh file
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//'external_mesh.bin'
  open(unit=27,file=trim(prname_lp), &
          status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error: could not open database '
    print *,'path: ',trim(prname_lp)
    stop 'Error reading external mesh file'
  endif

  ! gets number of elements and global points for this partition
  read(27) NSPEC
  read(27) NGLOB

  close(27)

  ! user output
  if (myrank == 0) then
    print *,'summing kernels in event kernel directories:'
    print *,kernel_list(1:nker)
    print *
  endif

  ! synchronizes
  call synchronize_all()

  ! sums up kernels
  if (USE_ISO_KERNELS) then

    !  isotropic kernels
    if (myrank == 0) write(*,*) 'isotropic kernels: bulk_c, bulk_beta, rho'

    kernel_name = 'bulk_c_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'bulk_beta_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'rho_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

  else if (USE_ALPHA_BETA_RHO) then

    ! isotropic kernels
    if (myrank == 0) write(*,*) 'isotropic kernels: alpha, beta, rho'

    kernel_name = 'alpha_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'beta_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'rho_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)
  !!! Kai added !!!
  else if (USE_ALPHA_BETA_RHO_TISO) then

    ! transverse isotropic kernels 
    if (myrank == 0) write(*,*) 'transverse isotropic kernels: alphav, alphah, bulk_betav, bulk_betah,eta'

    !kernel_name = 'alphav_kernel'
    !call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    !kernel_name = 'alphah_kernel'
    !call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'betav_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'betah_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    !kernel_name = 'eta_kernel'
    !call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    !kernel_name = 'alpha_kernel'
    !call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    !kernel_name = 'beta_kernel'
    !call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    !kernel_name = 'rho_kernel'
    !call sum_kernel(output_dir,kernel_name,kernel_list,nker)

  !!! Kai  !!!

  else

    ! transverse isotropic kernels
    if (myrank == 0) write(*,*) 'transverse isotropic kernels: bulk_c, bulk_betav, bulk_betah,eta'

    kernel_name = 'bulk_c_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'bulk_betav_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'bulk_betah_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'eta_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

  endif

  if (myrank == 0) write(*,*) 'done writing all kernels, see directory'//trim(output_dir)

  ! stop all the processes, and exit
  call finalize_mpi()

end program sum_kernels

!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel(output_dir,kernel_name,kernel_list,nker)

  use tomography_par

  implicit none

  character(len=MAX_STRING_LEN) :: output_dir,kernel_name,kernel_list(MAX_KERNEL_PATHS)
  integer :: nker

  ! local parameters
  character(len=MAX_STRING_LEN*2) :: k_file
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel,total_kernel
  double precision :: norm,norm_sum
  integer :: ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: mask_source

  ! initializes arrays
  allocate(kernel(NGLLX,NGLLY,NGLLZ,NSPEC), &
           total_kernel(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel arrays'

  if (USE_SOURCE_MASK) then
    allocate( mask_source(NGLLX,NGLLY,NGLLZ,NSPEC) )
    mask_source(:,:,:,:) = 1.0_CUSTOM_REAL
  endif

  ! loops over all event kernels
  total_kernel = 0._CUSTOM_REAL
  do iker = 1, nker
    ! user output
    if (myrank == 0) then
      write(*,*) 'reading in event kernel for: ',trim(kernel_name)
      write(*,*) '    ',iker, ' out of ', nker
    endif

    ! sensitivity kernel / frechet derivative
    kernel = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') trim(kernel_list(iker)) &
                          //'/proc',myrank,trim(REG)//trim(kernel_name)//'.bin'

    open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ier)
    if (ier /= 0) then
      write(*,*) '  kernel not found: ',trim(k_file)
      stop 'Error kernel file not found'
    endif
    read(IIN) kernel
    close(IIN)

    ! outputs norm of kernel
    norm = sum( kernel * kernel )
    call sum_all_dp(norm, norm_sum)
    if (myrank == 0) then
      print *,'  norm kernel: ',sqrt(norm_sum)
      print *
    endif

    ! source mask
    if (USE_SOURCE_MASK) then
      ! reads in mask
      write(k_file,'(a,i6.6,a)') trim(kernel_list(iker)) &
                            //'/proc',myrank,trim(REG)//'mask_source.bin'
      open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0) then
        write(*,*) '  file not found: ',trim(k_file)
        stop 'Error source mask file not found'
      endif
      read(IIN) mask_source
      close(IIN)

      ! masks source elements
      kernel = kernel * mask_source
    endif

    ! sums all kernels from each event
    total_kernel = total_kernel + kernel
  enddo

  ! stores summed kernels
  if (myrank == 0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name)

  write(k_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,trim(REG)//trim(kernel_name)//'.bin'

  open(IOUT,file=trim(k_file),form='unformatted',status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'Error kernel not written:',trim(k_file)
    stop 'Error kernel write'
  endif
  write(IOUT) total_kernel
  close(IOUT)

  if (myrank == 0) write(*,*)

  ! frees memory
  deallocate(kernel,total_kernel)
  if (USE_SOURCE_MASK) deallocate(mask_source)

end subroutine sum_kernel


