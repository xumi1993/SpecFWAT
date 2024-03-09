!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
subroutine compute_ctrl_grp(ekernel_dir_list,output_dir,model,nctrl)

  use tomography_par, only: MAX_STRING_LEN,MAX_KERNEL_PATHS,IIN, &
    NGLLX,NGLLY,NGLLZ,CUSTOM_REAL, &
    myrank,sizeprocs,NGLOB,NSPEC,USE_ALPHA_BETA_RHO,USE_ALPHA_BETA_RHO_TISO,USE_ISO_KERNELS
  !!! for reading Database
  use specfem_par
  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  use specfem_par_poroelastic, only: ispec_is_poroelastic

  
  use tomography_kernels_iso, only: total_kernel
  use shared_parameters
  use constants, only: PI

  implicit none

  real(kind=4), parameter       :: angle_thold=22.5
  character(len=MAX_STRING_LEN) :: model 
  character(len=MAX_STRING_LEN) :: kernel_list(MAX_KERNEL_PATHS)
  character(len=MAX_STRING_LEN) :: sline, kernel_name
  character(len=MAX_STRING_LEN) :: output_dir,ekernel_dir_list
  integer :: nker
  integer :: ier

  logical :: BROADCAST_AFTER_READ
  ! for selecting control group from a mini-batch
  character(len=MAX_STRING_LEN*2) :: k_file
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel,ctrl_kernel
  character(len=MAX_STRING_LEN) :: kernel_list_left(MAX_KERNEL_PATHS)
  integer :: iker,ictrl,nctrl,iker_reject_loc
  integer :: glob_num_left(MAX_KERNEL_PATHS),iker_reject_glob
  real(kind=CUSTOM_REAL) :: norm_tot, norm_sum_tot,norm_ctrl,norm_sum_ctrl, &
                            norm_tot_ctrl, norm_sum_tot_ctrl
  real(kind=CUSTOM_REAL) :: angle_diff,angle_diff_max
  ! write station list of the control group
   character(len=50) :: kname
   integer :: i,nsrc
   real(kind=4) :: lat,lon,ele,bur  
   character(len=50), dimension(:), allocatable  :: evnm
   real(kind=4), dimension(:), allocatable       :: evla,evlo,evdp,evbur

  ! ============ program starts here =====================

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  !call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  !call get_command_argument(1,ekernel_dir_list)
  !call get_command_argument(2,output_dir)

  !if (trim(ekernel_dir_list) == '' &
  !   .or. trim(output_dir) == '' ) then
  !   call exit_mpi(myrank,'USAGE: xsum_kernels ekernel_dir_list output_dir')
  !endif


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
  !if (sizeprocs /= NPROC) then
  !  if (myrank == 0) then
  !    print *
  !    print *,'Error: run xsum_kernels with the same number of MPI processes '
  !    print *,'       as specified in Par_file by NPROC when slices were created'
  !    print *
  !    print *,'for example: mpirun -np ',NPROC,' ./xsum_kernels ...'
  !    print *
  !  endif
  !  call synchronize_all()
  !  stop 'Error total number of slices'
  !endif
  call synchronize_all()

  ! reads mesh file
  !
  ! needs to get array dimensions

  !-------------------------------------------------------------------------
  ! opens external mesh file
  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some
  ! array sizes below
  call read_mesh_for_init()

  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),irregular_element_number(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 980')

  if (NSPEC_IRREGULAR > 0) then
    ! allocate arrays for storing the databases
    allocate(xix(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 981')
    allocate(xiy(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 982')
    allocate(xiz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 983')
    allocate(etax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 984')
    allocate(etay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 985')
    allocate(etaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 986')
    allocate(gammax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 987')
    allocate(gammay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 988')
    allocate(gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 989')
    allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 990')
   else
       ! allocate arrays for storing the databases
    allocate(xix(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 991')
    allocate(xiy(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 992')
    allocate(xiz(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 993')
    allocate(etax(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 994')
    allocate(etay(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 995')
    allocate(etaz(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 996')
    allocate(gammax(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 997')
    allocate(gammay(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 998')
    allocate(gammaz(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 999')
    allocate(jacobian(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1000')
  endif
  ! mesh node locations
  allocate(xstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1001')
  allocate(ystore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1002')
  allocate(zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1003')
  if (ier /= 0) stop 'Error allocating arrays for mesh nodes'

  ! material properties
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1004')
  allocate(mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1005')
  if (ier /= 0) stop 'Error allocating arrays for material properties'

  ! material flags
  allocate(ispec_is_acoustic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1006')
  allocate(ispec_is_elastic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1007')
  allocate(ispec_is_poroelastic(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1008')
  if (ier /= 0) stop 'Error allocating arrays for material flags'
  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.


  ! reads in external mesh
  call read_mesh_databases()

  NSPEC=NSPEC_AB
  NGLOB=NGLOB_AB
! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
 
  !-------------------------------------------------------------------------
  ! initializes arrays
  allocate(kernel(NGLLX,NGLLY,NGLLZ,NSPEC), &
           total_kernel(NGLLX,NGLLY,NGLLZ,NSPEC), &
           ctrl_kernel(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel arrays'

  ! user output
  ! if (myrank == 0) then
  !   print *,'summing kernels in event kernel directories:'
  !   print *,trim(kernel_list(1:nker))
  !   print *
  ! endif

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
    if (myrank == 0) write(*,*) 'isotropic kernels: beta'

    kernel_name = 'beta_kernel'
    call sum_kernel_array(kernel_name,kernel_list,nker)

    ! loop over each kernel added to control group
    !nctrl=nker/2 ! set size of the control group to be half of the mini-batch size
    if (nker <2*nctrl) stop 'nker should be at least 2*nctrl'
    glob_num_left=0
    do iker=1,nker
       glob_num_left(iker)=iker
    enddo
    kernel_list_left=kernel_list
    ctrl_kernel=total_kernel ! Inital ctrol group gradient = mini-batch gradient
    call calc_inner_product(total_kernel,total_kernel,norm_tot)
    call sum_all_cr(norm_tot,norm_sum_tot)
    if (myrank == 0) then
       norm_sum_tot = sqrt(norm_sum_tot)
    endif
 
    do ictrl=1, nker-nctrl
       if (myrank==0) print *, 'ictrl= ',ictrl
       ! loop over each kernel remained in the mini-batch
       angle_diff_max=-1000.
       iker_reject_loc=0
       do iker=1, nker-ictrl+1
          kernel = 0._CUSTOM_REAL
          write(k_file,'(a,i6.6,a)') trim(kernel_list_left(iker)) &
                          //'/proc',myrank,'_'//trim(kernel_name)//'.bin'

          open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ier)
          if (ier /= 0) then
          write(*,*) '  kernel not found: ',trim(k_file)
          stop 'Error kernel file not found'
          endif
          
          read(IIN) kernel
          close(IIN)
          ! remove one kernel every time
          kernel = ctrl_kernel - kernel
          call calc_inner_product(kernel,kernel,norm_ctrl)
          call sum_all_cr(norm_ctrl,norm_sum_ctrl)
          call calc_inner_product(total_kernel,kernel,norm_tot_ctrl)
          call sum_all_cr(norm_tot_ctrl,norm_sum_tot_ctrl)
          if (myrank == 0) then
             norm_sum_ctrl = sqrt(norm_sum_ctrl)
             !print *, 'norm_sum_tot,norm_sum_ctrl=',norm_sum_tot,norm_sum_ctrl
             !print *, 'norm_sum_tot_ctrl',norm_sum_tot_ctrl
             angle_diff = norm_sum_tot_ctrl / (norm_sum_ctrl * norm_sum_tot) 
             print *, '  iker,iker_glob,angle_diff,kername= ',iker,glob_num_left(iker), &
                      dble(acos(angle_diff))/PI * 180.,trim(kernel_list_left(iker))
             if (angle_diff_max<angle_diff) then
                angle_diff_max=angle_diff
                iker_reject_loc= iker
             endif
          endif
       enddo 
       !update ctrl_kernel and kernel_list_left 
       call bcast_all_singlei(iker_reject_loc)
       iker_reject_glob=glob_num_left(iker_reject_loc)
       if (myrank==0) write(*,*) "max angle diff=",dble(acos(angle_diff_max))/PI * 180.
       call bcast_all_singlecr(angle_diff_max)
       if (dble(acos(angle_diff_max))/PI * 180. > angle_thold) then
          if (myrank==0) print *, "Terminate. Angle diff reach the thresold ", angle_thold 
          nctrl=nker-ictrl
          exit
       endif
       if (myrank==0) write(*,*) " Reject ",iker_reject_glob," event kernel ",trim(kernel_list_left(iker_reject_loc))
       kernel = 0._CUSTOM_REAL
       write(k_file,'(a,i6.6,a)') trim(kernel_list(iker_reject_glob)) &
                       //'/proc',myrank,'_'//trim(kernel_name)//'.bin'

       open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ier)
       if (ier /= 0) then
       write(*,*) '  kernel not found: ',trim(k_file)
       stop 'Error kernel file not found'
       endif
       
       read(IIN) kernel
       close(IIN)
       ctrl_kernel=ctrl_kernel - kernel

       do iker=iker_reject_loc, nker-ictrl+1 
          glob_num_left(iker)=glob_num_left(iker+1)
          kernel_list_left(iker)=kernel_list_left(iker+1)
       enddo 

    enddo
    ! write station list of the control group
    kernel_name='src_rec/sources_batch.'//trim(model)//'.dat'
    open(9,file=trim(kernel_name),status='old')
    nsrc=0
    do
      read(9,*,iostat=ier) kname,lat,lon,ele,bur
      !write(*,*) trim(kname),lat,lon,ele,bur
      if (ier < 0 ) exit
      nsrc=nsrc+1
    enddo
    close(9)

    if (nsrc.ne.nker) then
      stop 'nker is not equal to nsrc of minibatch'
    endif
    allocate(evnm(nsrc))
    allocate(evla(nsrc))
    allocate(evlo(nsrc))
    allocate(evdp(nsrc))
    allocate(evbur(nsrc))
    open(9,file=trim(kernel_name),status='old')
    i=1
    do
      read(9,*,iostat=ier) evnm(i),evla(i),evlo(i),evdp(i),evbur(i)
      if (ier < 0 ) exit
      !write(*,*) trim(evnm(i)),evla(i),evlo(i),evdp(i),evbur(i)
      i=i+1
    enddo
    close(9)
 
    kernel_name='src_rec/sources_ctrlgrp.'//trim(model)//'.dat'
    open(9,file=trim(kernel_name))
    do i=1,nctrl
      write(9,"(A16,4F22.4)") evnm(glob_num_left(i)),evla(glob_num_left(i)),evlo(glob_num_left(i)), &
                       evdp(glob_num_left(i)),evbur(glob_num_left(i))
    enddo
    close(9)
    !!! chooose 8 events from the control group for line search
    kernel_name='src_rec/sources_ls.dat'
    open(9,file=trim(kernel_name))
    do i=1,int(0.8*nctrl)
      write(9,"(A16,4F22.4)") evnm(glob_num_left(i)),evla(glob_num_left(i)),evlo(glob_num_left(i)), &
                       evdp(glob_num_left(i)),evbur(glob_num_left(i))
    enddo
    close(9)
    deallocate(evnm,evla,evlo,evdp,evbur)
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
  !call finalize_mpi()
  ! frees memory
  deallocate(kernel,total_kernel)


end subroutine compute_ctrl_grp

subroutine calc_inner_product(vect1, vect2, q)
  use specfem_par, only: jacobian, irregular_element_number, &
                         wxgll, wygll, wzgll, NGLLX, NGLLY, NGLLZ, &
                         CUSTOM_REAL
  use tomography_par, only: NSPEC
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC), &
                          intent(in):: vect1, vect2
  real(kind=CUSTOM_REAL), intent(out) :: q
  ! local variables
  integer :: i,j,k,ispec,ispec_irreg
  real(kind=CUSTOM_REAL) :: jacobianl, weight, coeff_n1, coeff_n2
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
                          vect1n, vect2n

  ! nornalize
  coeff_n1 = maxval(abs(vect1(:,:,:,:)))
  if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1=1._CUSTOM_REAL
  vect1n(:,:,:,:) = vect1(:,:,:,:) / coeff_n1

  coeff_n2 = maxval(abs(vect2(:,:,:,:)))
  if (coeff_n2 == 0._CUSTOM_REAL) coeff_n2=1._CUSTOM_REAL
  vect2n(:,:,:,:) = vect2(:,:,:,:) / coeff_n2

  q = 0._CUSTOM_REAL

  do ispec = 1, NSPEC
    ispec_irreg = irregular_element_number(ispec)
    do k=1, NGLLZ; do j=1,NGLLY; do i=1,NGLLX
      weight = wxgll(i)*wygll(j)*wzgll(k)
      if (ispec_irreg /= 0) jacobianl = jacobian(i,j,k,ispec_irreg)
      q = q + jacobianl * weight * vect1n(i,j,k,ispec) * &
              vect2n(i,j,k,ispec)
    enddo;enddo;enddo
  enddo
  q = q * coeff_n1 * coeff_n2

end subroutine calc_inner_product
!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel_array(kernel_name,kernel_list,nker)

  use tomography_par
  use tomography_kernels_iso, only: total_kernel

  implicit none

  character(len=MAX_STRING_LEN) :: kernel_name,kernel_list(MAX_KERNEL_PATHS)
  integer :: nker

  ! local parameters
  character(len=MAX_STRING_LEN*2) :: k_file
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel!,total_kernel
  integer :: ier

  ! initializes arrays
  allocate(kernel(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel arrays'

  ! loops over all event kernels
  total_kernel = 0._CUSTOM_REAL
  do iker = 1, nker
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
    ! sums all kernels from each event
    total_kernel = total_kernel + kernel
  enddo

  ! frees memory
  deallocate(kernel)

end subroutine sum_kernel_array


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


subroutine sum_kernels_fwat(ekernel_dir_list,output_dir)

  use tomography_par, only: MAX_STRING_LEN,MAX_KERNEL_PATHS,IIN, &
                            myrank,sizeprocs,NGLOB,NSPEC,USE_ALPHA_BETA_RHO,USE_ISO_KERNELS,&
                            USE_ALPHA_BETA_RHO_TISO

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
  !call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  !call get_command_argument(1,ekernel_dir_list)
  !call get_command_argument(2,output_dir)

  !if (trim(ekernel_dir_list) == '' &
  !   .or. trim(output_dir) == '' ) then
  !   call exit_mpi(myrank,'USAGE: xsum_preconditioned_kernels ekernel_dir_list output_dir')
  !endif


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
  ! if (myrank == 0) then
  !   print *,'summing kernels in event kernel directories:'
  !   print *,trim(kernel_list(1:nker))
  !   print *
  ! endif

  ! synchronizes
  call synchronize_all()

  if (USE_ALPHA_BETA_RHO) then

    ! isotropic kernels
    if (myrank == 0) write(*,*) 'isotropic kernels: alpha, beta, rho'

    kernel_name = 'alpha_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'beta_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)

    kernel_name = 'rhop_kernel'
    call sum_kernel(output_dir,kernel_name,kernel_list,nker)
  endif
  if (myrank == 0) write(*,*) 'done writing all kernels, see directory '//trim(output_dir)
end subroutine sum_kernels_fwat


subroutine sum_preconditioned_kernels_fwat(ekernel_dir_list,output_dir,norm_hess)

  use tomography_par, only: MAX_STRING_LEN,MAX_KERNEL_PATHS,IIN, &
                            myrank,sizeprocs,NGLOB,NSPEC,USE_ALPHA_BETA_RHO,USE_ISO_KERNELS,&
                            USE_ALPHA_BETA_RHO_TISO

  use shared_parameters

  implicit none

  character(len=MAX_STRING_LEN) :: kernel_list(MAX_KERNEL_PATHS)
  character(len=MAX_STRING_LEN) :: sline, kernel_name,prname_lp
  character(len=MAX_STRING_LEN) :: output_dir,ekernel_dir_list
  integer :: nker
  integer :: ier

  logical :: BROADCAST_AFTER_READ,norm_hess

  ! ============ program starts here =====================

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  !call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  !call get_command_argument(1,ekernel_dir_list)
  !call get_command_argument(2,output_dir)

  !if (trim(ekernel_dir_list) == '' &
  !   .or. trim(output_dir) == '' ) then
  !   call exit_mpi(myrank,'USAGE: xsum_preconditioned_kernels ekernel_dir_list output_dir')
  !endif


  if (myrank == 0) then
     write(*,*) 'SUM EVENT KERNELS TO GET MISFIT KENRELS'
     write(*,*) 'INPUT EVENT DIRECTORY:',trim(ekernel_dir_list)
     write(*,*) 'OUTPUT DIRECTORY:',trim(output_dir)
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
  ! if (myrank == 0) then
  !   print *,'summing kernels in event kernel directories:'
  !   print *,trim(kernel_list(1:nker))
  !   print *
  ! endif

  ! synchronizes
  call synchronize_all()

  ! sums up kernels
  if (USE_ISO_KERNELS) then

    !  isotropic kernels
    if (myrank == 0) write(*,*) 'isotropic kernels: bulk_c, bulk_beta, rho'

    kernel_name = 'bulk_c_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'bulk_beta_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'rho_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

  else if (USE_ALPHA_BETA_RHO) then

    ! isotropic kernels
    if (myrank == 0) write(*,*) 'isotropic kernels: alpha, beta, rho'

    kernel_name = 'alpha_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'beta_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'rhop_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)
  !!! Kai added !!!
  else if (USE_ALPHA_BETA_RHO_TISO) then

    ! transverse isotropic kernels 
    if (myrank == 0) write(*,*) 'transverse isotropic kernels: alphav, alphah, bulk_betav, bulk_betah,eta'

    kernel_name = 'alphav_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'alphah_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'betav_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'betah_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'eta_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'rho_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

  !!! Kai  !!!

  else

    ! transverse isotropic kernels
    if (myrank == 0) write(*,*) 'transverse isotropic kernels: bulk_c, bulk_betav, bulk_betah,eta'

    kernel_name = 'bulk_c_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'bulk_betav_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'bulk_betah_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

    kernel_name = 'eta_kernel'
    call sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

  endif

  if (myrank == 0) write(*,*) 'done writing all kernels, see directory '//trim(output_dir)

  ! stop all the processes, and exit
  !call finalize_mpi()

end subroutine sum_preconditioned_kernels_fwat

!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel_pre(output_dir,kernel_name,kernel_list,nker,norm_hess)

  use tomography_par

  implicit none

  character(len=MAX_STRING_LEN) :: output_dir,kernel_name,kernel_list(MAX_KERNEL_PATHS)
  integer :: nker

  ! local parameters
  character(len=MAX_STRING_LEN*2) :: k_file
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel,hess,total_kernel
  double precision :: norm,norm_sum
  integer :: ier
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: total_hess,mask_source
  logical :: norm_hess

  ! initializes arrays
  allocate(kernel(NGLLX,NGLLY,NGLLZ,NSPEC), &
           hess(NGLLX,NGLLY,NGLLZ,NSPEC), &
           total_kernel(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel arrays'

  if (USE_HESS_SUM) then
    allocate( total_hess(NGLLX,NGLLY,NGLLZ,NSPEC) )
    total_hess(:,:,:,:) = 0.0_CUSTOM_REAL
  endif

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
      write(*,*) 'and preconditioner         : ','hess_kernel'
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
      print *,'  norm kernel        : ',sqrt(norm_sum)
    endif

    ! approximate Hessian
    hess = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') trim(kernel_list(iker)) &
                          //'/proc',myrank,trim(REG)//'hess_kernel.bin'

    open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ier)
    if (ier /= 0) then
      write(*,*) '  Hessian kernel not found: ',trim(k_file)
      stop 'Error hess_kernel.bin files not found'
    endif

    read(IIN) hess
    close(IIN)

    ! outputs norm of preconditioner
    norm = sum( hess * hess )
    call sum_all_dp(norm, norm_sum)
    if (myrank == 0) then
      print *,'  norm preconditioner: ',sqrt(norm_sum)
    endif

    ! note: we take absolute values for Hessian (as proposed by Yang)
    hess = abs(hess)

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

    ! precondition
    if (USE_HESS_SUM) then

      ! sums up Hessians first
      total_hess = total_hess + hess

    else

      ! inverts Hessian
      call invert_hess( hess )

      ! preconditions each event kernel with its Hessian
      kernel = kernel * hess

    endif

    ! sums all kernels from each event
    total_kernel = total_kernel + kernel

    if (myrank == 0) print *
  enddo

  ! preconditions summed kernels with summed Hessians
  if (USE_HESS_SUM) then
    ! inverts Hessian matrix
    call invert_hess( total_hess )

    if (norm_hess) then
      call synchronize_all()
      norm = dble(maxval( abs(total_hess) ))
      call max_all_dp(norm, norm_sum)
      call bcast_all_singledp(norm_sum)
      total_hess = total_hess/norm_sum
    endif

    ! preconditions kernel
    total_kernel = total_kernel * total_hess

  endif
  norm = sum( total_kernel * total_kernel )
  call sum_all_dp(norm, norm_sum)
  if (myrank == 0) then
    print *,'  total_kernel      : ',sqrt(norm_sum)
  endif

  ! stores summed kernels
  if (myrank == 0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name)

  ! outputs summed kernel
  write(k_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,trim(REG) // trim(kernel_name) // '.bin'
  open(IOUT,file=trim(k_file),form='unformatted',status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'Error kernel not written: ',trim(k_file)
    stop 'Error kernel write'
  endif
  write(IOUT) total_kernel
  close(IOUT)

  ! outputs summed Hessian
  if (USE_HESS_SUM) then
    if (myrank == 0) write(*,*) 'writing out summed kernel for: ','hess_inv_kernel'
    write(k_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,trim(REG) // 'hess_inv_kernel' // '.bin'
    open(IOUT,file=trim(k_file),form='unformatted',status='unknown',action='write',iostat=ier)
    if (ier /= 0) then
      write(*,*) 'Error kernel not written: ',trim(k_file)
      stop 'Error kernel write'
    endif
    write(IOUT) total_hess
    close(IOUT)
  endif

  if (myrank == 0) write(*,*)

  ! frees memory
  deallocate(kernel,hess,total_kernel)
  if (USE_HESS_SUM) deallocate(total_hess)
  if (USE_SOURCE_MASK) deallocate(mask_source)

end subroutine sum_kernel_pre

!
!-------------------------------------------------------------------------------------------------
!

subroutine invert_hess( hess_matrix )

! inverts the Hessian matrix
! the approximate Hessian is only defined for diagonal elements: like
! H_nn = \frac{ \partial^2 \chi }{ \partial \rho_n \partial \rho_n }
! on all GLL points, which are indexed (i,j,k,ispec)

  use tomography_par
  use fullwave_adjoint_tomo_par, only: OUT_FWAT_LOG, PRECOND_TYPE

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: hess_matrix

  ! local parameters
  real(kind=CUSTOM_REAL) :: maxh,maxh_all

  ! maximum value of Hessian
  maxh = maxval( abs(hess_matrix) )

  ! determines maximum from all slices on master
  call max_all_all_cr(maxh, maxh_all)

  ! user output
  if (myrank == 0) then
    print *
    print *,'Hessian maximum: ',maxh_all
    print *
  endif

  ! normalizes Hessian
  if (maxh_all < 1.e-18) then
    ! Hessian is zero, re-initializes
    hess_matrix = 1.0_CUSTOM_REAL
    !stop 'Error Hessian too small'
  else
    ! since Hessian has absolute values, this scales between [0,1]
    hess_matrix = hess_matrix / maxh_all
  endif


  ! inverts Hessian values
  if (PRECOND_TYPE=='default') then
    ! write(*,*) 'USE THRESHOLD_HESS: ',THRESHOLD_HESS
    ! write(OUT_FWAT_LOG,*) 'USE THRESHOLD_HESS: ',THRESHOLD_HESS
    where( abs(hess_matrix(:,:,:,:)) > THRESHOLD_HESS )
      hess_matrix = 1.0_CUSTOM_REAL / hess_matrix
    elsewhere
      hess_matrix = 1.0_CUSTOM_REAL / THRESHOLD_HESS
    endwhere
  else
    ! write(OUT_FWAT_LOG,*) 'NO THRESHOLD_HESS'
    hess_matrix = 1.0_CUSTOM_REAL / hess_matrix
  endif


  ! rescales Hessian
  !hess_matrix = hess_matrix * maxh_all

end subroutine invert_hess

!----------------------------------------------------------------------------------

subroutine smooth_sem_fwat(sigma_h,sigma_v,kernel_names_comma_delimited,input_dir,&
                      output_dir,USE_GPU)

  use constants, only: USE_QUADRATURE_RULE_FOR_SMOOTHING
  use postprocess_par, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,NGLLSQUARE, &
    MAX_STRING_LEN,IIN,IOUT,GAUSSALPHA,GAUSSBETA,PI,TWO_PI,MAX_KERNEL_NAMES

  use specfem_par
  use specfem_par_movie
  use taper3d
  use fullwave_adjoint_tomo_par, only: OUT_FWAT_LOG
  use fwat_input, only: tomo_par

  implicit none

  integer, parameter :: NARGS = 6

  ! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat,dat_smooth
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dummy ! for jacobian read
  integer :: NSPEC_N, NGLOB_N ,NSPEC_IRREGULAR_N !BinHe

  integer :: i,j,k,iglob,ier,ispec2,ispec,inum,ispec_irreg
  integer :: icounter,num_slices
  integer :: iproc,ncuda_devices
  !integer(kind=8) :: Container

  integer,parameter :: MAX_NODE_LIST = 300
  integer :: node_list(MAX_NODE_LIST)
  logical :: do_include_slice

  !character(len=MAX_STRING_LEN) :: arg(6)
  character(len=MAX_STRING_LEN) :: kernel_name, input_dir, output_dir
  character(len=MAX_STRING_LEN) :: prname_lp
  character(len=MAX_STRING_LEN*2) :: local_data_file


  !character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  integer :: nker
  real t1,t2

  ! smoothing parameters
  character(len=MAX_STRING_LEN*2) :: ks_file

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2, sigma_h3, sigma_v, sigma_v2, sigma_v3
  real(kind=CUSTOM_REAL) :: sigma_h2_inv,sigma_v2_inv
  real(kind=CUSTOM_REAL) :: sigma_h3_sq,sigma_v3_sq

  real(kind=CUSTOM_REAL) :: x0, y0, z0, norm, norm_h, norm_v
  real(kind=CUSTOM_REAL) :: center_x0, center_y0, center_z0
  real(kind=CUSTOM_REAL) :: center_x, center_y, center_z

  real(kind=CUSTOM_REAL) :: max_old,max_new,max_old_all,max_new_all
  real(kind=CUSTOM_REAL) :: min_old,min_new,min_old_all,min_new_all
  real(kind=CUSTOM_REAL) :: max_abs_all
  real(kind=CUSTOM_REAL), dimension(2) :: max_abs

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: exp_val,factor
  real(kind=CUSTOM_REAL) :: jacobianl

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: tk, bk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: xl, yl, zl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: xx, yy, zz

  real(kind=CUSTOM_REAL), dimension(:),allocatable :: cx0, cy0, cz0
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: cx, cy, cz

  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: element_size

  ! reference slice
  real(kind=CUSTOM_REAL) :: x_min_ref,x_max_ref
  real(kind=CUSTOM_REAL) :: y_min_ref,y_max_ref
  real(kind=CUSTOM_REAL) :: z_min_ref,z_max_ref
  real(kind=CUSTOM_REAL) :: x_min,x_max
  real(kind=CUSTOM_REAL) :: y_min,y_max
  real(kind=CUSTOM_REAL) :: z_min,z_max
  real(kind=CUSTOM_REAL) :: dim_x,dim_y,dim_z

  logical :: BROADCAST_AFTER_READ, USE_GPU
  ! type(taper_cls) :: taper_3d
  ! real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: taper_val

  !call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) print *,"Running XSMOOTH_SEM"
  ! if (myrank == 0) write(OUT_FWAT_LOG, *) "Running XSMOOTH_SEM ", trim(kernel_names_comma_delimited)
  call synchronize_all()
  call cpu_time(t1)

  ! parse command line arguments
  !if (command_argument_count() /= NARGS) then
  !  if (myrank == 0) then
  !      print *,'USAGE:  mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR GPU_MODE'
  !    stop 'Please check command line arguments'
  !  endif
  !endif
  call synchronize_all()

  !do i = 1, NARGS
  !  call get_command_argument(i,arg(i), status=ier)
  !enddo

  !read(arg(1),*) sigma_h
  !read(arg(2),*) sigma_v
  !kernel_names_comma_delimited = arg(3)
  !input_dir= arg(4)
  !output_dir = arg(5)
  !read(arg(6),*) USE_GPU

  !call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  kernel_name = trim(kernel_names_comma_delimited)

  if (USE_GPU) call initialize_cuda_device(myrank,ncuda_devices)

  nker = 1
  if (nker > 1) then
    if (myrank == 0) then
      ! The machinery for reading multiple names from the command line is in place,
      ! but the smoothing routines themselves have not yet been modified to work
      !  on multiple arrays.
      if (myrank == 0) then
        print *,'Smoothing only first name in list: ',trim(kernel_name)
        print *
      endif
    endif
  endif
  call synchronize_all()

  ! check smoothing radii
  sigma_h2 = 2.0 * sigma_h ** 2  ! factor two for Gaussian distribution with standard variance sigma
  sigma_v2 = 2.0 * sigma_v ** 2

  if (sigma_h2 < 1.e-18) stop 'Error sigma_h2 zero, must non-zero'
  if (sigma_v2 < 1.e-18) stop 'Error sigma_v2 zero, must non-zero'

  ! adds margin to search radius
  element_size = max(sigma_h,sigma_v) * 0.5

  ! search radius
  sigma_h3 = 3.0  * sigma_h + element_size
  sigma_v3 = 3.0  * sigma_v + element_size

  ! helper variables
  sigma_h2_inv = 1.0_CUSTOM_REAL / sigma_h2
  sigma_v2_inv = 1.0_CUSTOM_REAL / sigma_v2

  sigma_h3_sq = sigma_h3 * sigma_h3
  sigma_v3_sq = sigma_v3 * sigma_v3

  ! theoretic normal value
  ! (see integral over -inf to +inf of exp[- x*x/(2*sigma) ] = sigma * sqrt(2*pi) )
  ! note: smoothing is using a Gaussian (ellipsoid for sigma_h /= sigma_v),
  norm_h = real(2.0*PI*sigma_h**2,kind=CUSTOM_REAL)
  norm_v = real(sqrt(2.0*PI) * sigma_v,kind=CUSTOM_REAL)
  norm   = norm_h * norm_v

  ! user output
  if (myrank == 0) then
    print *,'command line arguments:'
    print *,'  smoothing sigma_h , sigma_v                : ',sigma_h,sigma_v
    ! scale length: approximately S ~ sigma * sqrt(8.0) for a Gaussian smoothing
    print *,'  smoothing scalelengths horizontal, vertical: ',sigma_h*sqrt(8.0),sigma_v*sqrt(8.0)
    print *,'  input dir : ',trim(input_dir)
    print *,'  output dir: ',trim(output_dir)
    print *,"  GPU_MODE: ", USE_GPU
    print *
  endif

  ! reads the parameter file
  ! BROADCAST_AFTER_READ = .true.
  ! call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  if (ADIOS_ENABLED) stop 'Flag ADIOS_ENABLED not supported yet for smoothing, please rerun program...'

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *,'Error number of processors supposed to run on: ',NPROC
      print *,'Error number of MPI processors actually run on: ',sizeprocs
      print *
      print *,'Please rerun with: mpirun -np ',NPROC,' bin/xsmooth_sem .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif

  ! for smoothing, we use cell centers to find and locate nearby elements
  !
  ! sets the location of the center of the elements and local points
  allocate(xl(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           yl(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           zl(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
          !  taper_val(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           cx0(NSPEC_AB), &
           cy0(NSPEC_AB), &
           cz0(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array xl etc.'

  ! sets element center location
  do ispec = 1, nspec_AB

    do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX

      iglob = ibool(i,j,k,ispec)
      xl(i,j,k,ispec) = xstore(iglob)
      yl(i,j,k,ispec) = ystore(iglob)
      zl(i,j,k,ispec) = zstore(iglob)
    enddo; enddo; enddo ! NGLLZ,NGLLY,NGLLX

    cx0(ispec) = (xl(1,1,1,ispec) + xl(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
    cy0(ispec) = (yl(1,1,1,ispec) + yl(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
    cz0(ispec) = (zl(1,1,1,ispec) + zl(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
  enddo

  ! reference slice dimension
  x_min_ref = minval(xstore)
  x_max_ref = maxval(xstore)
  y_min_ref = minval(ystore)
  y_max_ref = maxval(ystore)
  z_min_ref = minval(zstore)
  z_max_ref = maxval(zstore)

  dim_x = abs(x_max_ref - x_min_ref)
  dim_y = abs(y_max_ref - y_min_ref)
  dim_z = abs(z_max_ref - z_min_ref)

  ! user output
  if (myrank == 0) then
    print *,'smoothing:'
    print *,'  single slice dimensions in x/y/z-direction: ',dim_x,' / ',dim_y,' / ',dim_z
    print *,'  Gaussian search radius horizontal =',sigma_h3,' vertical =',sigma_v3
    print *
  endif

  ! checks if Gaussian support exceeds slice dimensions
  if (sigma_h3 >= dim_x .or. sigma_h3 >= dim_y .or. sigma_v3 >= dim_z) then
    ! Gaussian support is likely larger than the direct neighbour and has support in a much wider area
    ! user output
    if (myrank == 0) then
      print *,'  using large Gaussian with respect to slice dimension for smoothing'
      print *
    endif
  endif

  ! frees memory
  deallocate(xstore,ystore,zstore)
  deallocate(ibool)
  !deallocate(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
  if(allocated(irregular_element_number)) deallocate(irregular_element_number)
  if(allocated(jacobian)) deallocate(jacobian)
  if(allocated(dummy)) deallocate(dummy)

  call synchronize_all()

  ! initializes node list
  node_list(:) = -1

  ! adds this partition itself first
  icounter = 1
  node_list(icounter) = myrank

  ! sets up slices to process
  do iproc = 0,NPROC-1
    ! skip own process slice, has already been added
    if (iproc == myrank) cycle

    ! checks if slice is a direct neighbor
    do_include_slice = .false.
    do i = 1,num_interfaces_ext_mesh
      if (iproc == my_neighbors_ext_mesh(i)) then
        ! found a neighbor slice
        do_include_slice = .true.
        exit
      endif
    enddo
    if (.not. do_include_slice) then
      ! note: Gaussian support might be larger than closest neighbor slices
      !       we add all slices close enough to still have an influence

      ! checks distances to this slice
      ! reads in slice mesh
      ! neighbor database file
      call create_name_database(prname,iproc,LOCAL_PATH)
      prname_lp = prname(1:len_trim(prname))//'external_mesh.bin'

      ! gets number of elements and global points for this partition
      open(unit=IIN,file=trim(prname_lp),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'Error could not open database file: ',trim(prname_lp)
        call exit_mpi(myrank, 'Error reading neighbors external mesh file')
      endif
      read(IIN) NSPEC_N
      read(IIN) NGLOB_N
      read(IIN) NSPEC_IRREGULAR_N!BinHe

      ! allocates mesh arrays
      allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
      if (ier /= 0) stop 'Error allocating array ibool'
      allocate(xstore(NGLOB_N),ystore(NGLOB_N),zstore(NGLOB_N),stat=ier)
      if (ier /= 0) stop 'Error allocating array xstore etc.'

      ! ibool file
      read(IIN) ibool

      ! global point arrays
      read(IIN) xstore
      read(IIN) ystore
      read(IIN) zstore
      close(IIN)

      ! determines min/max values of slice mesh
      x_min = minval(xstore)
      x_max = maxval(xstore)
      y_min = minval(ystore)
      y_max = maxval(ystore)
      z_min = minval(zstore)
      z_max = maxval(zstore)

      ! slice dimensions
      dim_x = abs(x_max - x_min)
      dim_y = abs(y_max - y_min)
      dim_z = abs(z_max - z_min)

      ! frees memory
      deallocate(ibool)
      deallocate(xstore,ystore,zstore)

      ! re-evalutes distances between slices
      ! checks if edges are within search radius from reference slice
      ! note: for slices with irregular shapes (e.g. from scotch decomposition), this assumes a box slice shape
      !       and might add slices even if further separated; a conservative approach
      if ((x_min < x_max_ref + sigma_h3 .and. x_max > x_min_ref - sigma_h3) &
          .and. (y_min < y_max_ref + sigma_h3 .and. y_max > y_min_ref - sigma_h3) &
          .and. (z_min < z_max_ref + sigma_v3 .and. z_max > z_min_ref - sigma_v3)) then
        do_include_slice = .true.
      endif

    endif

    ! adds to smoothing list neighbors
    if (do_include_slice) then
      icounter = icounter + 1

      ! checks bounds
      if (icounter > MAX_NODE_LIST) stop 'Error number of interfaces exceeds MAX_NODE_LIST'

      ! adds slice to smoothing list
      node_list(icounter) = iproc
    endif
  enddo
  ! set total number of slices to loop over
  num_slices = icounter

  ! synchronizes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    print *,'slices:',num_slices
    print *,'  rank ',myrank,'  has smoothing slices:'
    print *,node_list(1:num_slices)
    print *
  endif

  !do i=0,sizeprocs-1
  !  if (myrank == i) then
  !    print *,'rank:',myrank,'  smoothing slices'
  !    print *,node_list(1:num_slices)
  !    print *
  !  endif
  !enddo

  ! for jacobian and weights
  ! GLL points weights
  if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
    call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          wgll_cube(i,j,k) = wxgll(i)*wygll(j)*wzgll(k)
        enddo
      enddo
    enddo
  endif

  ! loops over slices
  ! each process reads in his own neighbor slices and Gaussian filters the values
  allocate(tk(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           bk(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array tk and bk'

  tk = 0.0_CUSTOM_REAL
  bk = 0.0_CUSTOM_REAL

  ! GPU setup
  !if (USE_GPU) then
  !  call prepare_GPU_smooth(Container,xl,yl,zl,sigma_h2_inv,sigma_v2_inv,sigma_h3_sq,sigma_v3_sq,NSPEC_AB,nker,wgll_cube)

  !  ! synchronizes all processes
  !  call synchronize_all()
  !endif


  do inum = 1,num_slices

    iproc = node_list(inum)

    if (myrank == 0) print *,'  reading slice:',iproc

    ! neighbor database file
    call create_name_database(prname,iproc,LOCAL_PATH)
    prname_lp = prname(1:len_trim(prname))//'external_mesh.bin'

    ! gets number of point locations
    open(unit=IIN,file=trim(prname_lp),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname_lp)
      call exit_mpi(myrank, 'Error reading neighbors external mesh file')
    endif
    !===========BinHe==========! 
    !copy from src/tomography/post/sem_smooth.F90
    read(IIN) NSPEC_N
    read(IIN) NGLOB_N
    read(IIN) NSPEC_IRREGULAR_N
    ! allocates arrays
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1023')
    if (ier /= 0) stop 'Error allocating array ibool'
    allocate(xstore(NGLOB_N),ystore(NGLOB_N),zstore(NGLOB_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1024')
    if (ier /= 0) stop 'Error allocating array xstore etc.'

    if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
      if (NSPEC_IRREGULAR_N > 0) then
        allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR_N),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1025')
        if (ier /= 0) stop 'Error allocating array jacobian'
        allocate(dummy(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1026')
        if (ier /= 0) stop 'Error allocating array dummy'
      else
        allocate(jacobian(1,1,1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1027')
        if (ier /= 0) stop 'Error allocating array jacobian'
        allocate(dummy(1,1,1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1028')
        if (ier /= 0) stop 'Error allocating array dummy'
      endif
      allocate(irregular_element_number(NSPEC_N),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1029')
      if (ier /= 0) stop 'Error allocating array irregular_element_number'
    endif

    ! ibool file
    read(IIN) ibool

    ! global point arrays
    read(IIN) xstore
    read(IIN) ystore
    read(IIN) zstore

    if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
    ! reads in jacobian
      read(IIN) irregular_element_number
      read(IIN) xix_regular
      read(IIN) jacobian_regular
      read(IIN) dummy ! xix
      read(IIN) dummy ! xiy
      read(IIN) dummy ! xiz
      read(IIN) dummy ! etax
      read(IIN) dummy ! etay
      read(IIN) dummy ! etaz
      read(IIN) dummy ! gammax
      read(IIN) dummy ! gammay
      read(IIN) dummy ! gammaz
      read(IIN) jacobian
    endif

    close(IIN)
    !===========BinHe==========! 

    ! get the location of the center of the elements and local points
    allocate(xx(NGLLX,NGLLY,NGLLZ,NSPEC_N), &
             yy(NGLLX,NGLLY,NGLLZ,NSPEC_N), &
             zz(NGLLX,NGLLY,NGLLZ,NSPEC_N), &
             cx(NSPEC_N), &
             cy(NSPEC_N), &
             cz(NSPEC_N),stat=ier)
    if (ier /= 0) stop 'Error allocating array xx etc.'

    ! sets element center location
    do ispec = 1, NSPEC_N

      do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX

        iglob = ibool(i,j,k,ispec)
        xx(i,j,k,ispec) = xstore(iglob)
        yy(i,j,k,ispec) = ystore(iglob)
        zz(i,j,k,ispec) = zstore(iglob)

      enddo; enddo; enddo ! NGLLZ,NGLLY,NGLLX

      ! calculate element center location
      cx(ispec) = (xx(1,1,1,ispec) + xx(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
      cy(ispec) = (yy(1,1,1,ispec) + yy(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
      cz(ispec) = (zz(1,1,1,ispec) + zz(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
    enddo

    deallocate(xstore,ystore,zstore)
    deallocate(ibool)

    ! data file
    write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_'
    local_data_file = trim(prname) // trim(kernel_name) // '.bin'

    open(unit = IIN,file = trim(local_data_file),status='old',action='read',form ='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening data file: ',trim(local_data_file)
      stop 'Error opening data file'
    endif

    allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) stop 'Error allocating dat array'

    read(IIN) dat
    close(IIN)

    ! statistics
    if (iproc == myrank) then
      min_old = minval(dat(:,:,:,:))
      max_old = maxval(dat(:,:,:,:))
    endif

    ! finds closest elements for smoothing
    !if (myrank==0) print *, '  start looping over elements and points for smoothing ...'
    if (USE_GPU) then
      stop
      !call compute_smooth(Container,jacobian,xx,yy,zz,dat,NSPEC_N)
    else
      ! loop over elements to be smoothed in the current slice
      do ispec = 1, NSPEC_AB

        ! element center position
        center_x0 = cx0(ispec)
        center_y0 = cy0(ispec)
        center_z0 = cz0(ispec)

        ! --- only double loop over the elements in the search radius ---
        do ispec2 = 1, NSPEC_N

          ! search element center position
          center_x = cx(ispec2)
          center_y = cy(ispec2)
          center_z = cz(ispec2)

          ! calculates horizontal and vertical distance between two element centers
          ! (squared distances)
          call get_distance_vec(dist_h,dist_v,center_x0,center_y0,center_z0,center_x,center_y,center_z)

          ! checks distance between centers of elements
          if (dist_h > sigma_h3_sq .or. dist_v > sigma_v3_sq) cycle

          ! integration factors
          if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
            ispec_irreg = irregular_element_number(ispec2)
            if (ispec_irreg == 0) jacobianl = jacobian_regular
            do k=1,NGLLZ
              do j=1,NGLLY
                do i=1,NGLLX
                  if (ispec_irreg /= 0) jacobianl = jacobian(i,j,k,ispec)
                  factor(i,j,k) = jacobianl * wgll_cube(i,j,k)
                enddo
              enddo
            enddo

          endif

          ! loop over GLL points of the elements in current slice (ispec)
          do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX

            ! reference location
            ! current point (i,j,k,ispec) location, Cartesian coordinates
            x0 = xl(i,j,k,ispec)
            y0 = yl(i,j,k,ispec)
            z0 = zl(i,j,k,ispec)

            ! calculate weights based on Gaussian smoothing
            call smoothing_weights_vec(x0,y0,z0,sigma_h2_inv,sigma_v2_inv,exp_val, &
                                       xx(:,:,:,ispec2),yy(:,:,:,ispec2),zz(:,:,:,ispec2))

            ! adds GLL integration weights
            if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
              exp_val(:,:,:) = exp_val(:,:,:) * factor(:,:,:)
            endif

            ! adds contribution of element ispec2 to smoothed kernel values
            tk(i,j,k,ispec) = tk(i,j,k,ispec) + sum(exp_val(:,:,:) * dat(:,:,:,ispec2))

            ! normalization, integrated values of Gaussian smoothing function
            bk(i,j,k,ispec) = bk(i,j,k,ispec) + sum(exp_val(:,:,:))

          enddo; enddo; enddo ! NGLLZ,NGLLY,NGLLX

        enddo ! ispec2
      enddo ! ispec

    endif ! GPU_MODE

    ! frees arrays
    deallocate(dat)
    deallocate(xx,yy,zz)
    deallocate(cx,cy,cz)

    if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
      deallocate(irregular_element_number)
      deallocate(jacobian)
      deallocate(dummy)
    endif

  enddo ! iproc

  ! normalizes/scaling factor
  if (myrank == 0) then
    print *
    print *,'Scaling values: min/max = ',minval(bk),maxval(bk)
  endif

  allocate(dat_smooth(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array dat_smooth'

  dat_smooth(:,:,:,:) = 0.0_CUSTOM_REAL

  if (USE_GPU) then
    !call get_smooth(Container,dat_smooth)
    stop
  else
    do ispec = 1, NSPEC_AB

      do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX

        ! checks the normalization criterion
        !if (abs(bk(i,j,k,ispec) - norm) > 1.e-4) then
        !  print *, 'Problem norm here --- ', ispec, i, j, k, bk(i,j,k,ispec), norm
        !endif
        if (abs(bk(i,j,k,ispec)) < 1.e-18) then
          print *, 'Problem norm here --- ', ispec, i, j, k, bk(i,j,k,ispec), norm
        endif

        ! normalizes smoothed kernel values by integral value of Gaussian weighting
        ! dat_smooth(i,j,k,ispec) = taper_val(i,j,k,ispec)* tk(i,j,k,ispec) / bk(i,j,k,ispec)
        dat_smooth(i,j,k,ispec) = tk(i,j,k,ispec) / bk(i,j,k,ispec)

      enddo; enddo; enddo ! NGLLZ,NGLLY,NGLLX

    enddo !  ispec

  endif ! GPU_MODE

  ! frees memory
  deallocate(tk,bk)
  deallocate(xl,yl,zl)
  deallocate(cx0,cy0,cz0)
  ! deallocate(taper_val)

  ! statistics
  min_new = minval(dat_smooth(:,:,:,:))
  max_new = maxval(dat_smooth(:,:,:,:))

  ! the min/maximum value for the smoothed kernel
  call min_all_cr(min_old, min_old_all)
  call min_all_cr(min_new, min_new_all)
  call max_all_cr(max_old, max_old_all)
  call max_all_cr(max_new,max_new_all)
  call bcast_all_singlecr(min_new_all)
  call bcast_all_singlecr(max_new_all)

  if (myrank == 0) then
    print *
    print *,'Minimum data value before smoothing = ', min_old_all
    print *,'Minimum data value after smoothing  = ', min_new_all
    print *
    print *,'Maximum data value before smoothing = ', max_old_all
    print *,'Maximum data value after smoothing  = ', max_new_all
    print *
    close(IMAIN)
  endif

  ! ! Mijian added for normalization
  ! max_abs(1) = abs(min_new_all)
  ! max_abs(2) = abs(max_new_all)
  ! max_abs_all = maxval(max_abs)
  ! do ispec = 1, NSPEC_AB
  !   do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX
  !     dat_smooth(i,j,k,ispec) = dat_smooth(i,j,k,ispec)/max_abs_all
  !   enddo;enddo;enddo
  ! enddo

  ! file output
  ! smoothed kernel file name
  write(ks_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_name)//'_smooth.bin'

  open(IOUT,file=trim(ks_file),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error opening smoothed kernel file'
  write(IOUT) dat_smooth(:,:,:,:)
  close(IOUT)
  if (myrank == 0) print *,'written: ',trim(ks_file)

  ! frees memory
  deallocate(dat_smooth)

  ! synchronizes
  call synchronize_all()

  call cpu_time(t2)

  !if (USE_GPU) then
  !  print *,'Computation time with GPU:',t2-t1
  !else
  !  print *,'Computation time with CPU:',t2-t1
  !endif

  ! stop all the processes and exit
  !call finalize_mpi()

  ! Kai added for restore back ibool, xstore, ystore zstore after smoothing
  if (myrank == 0) print *,' ..reading mesh back after smoothing for slice:',myrank

  ! neighbor database file
  call create_name_database(prname,myrank,LOCAL_PATH)
  prname_lp = prname(1:len_trim(prname))//'external_mesh.bin'

  ! gets number of point locations
  open(unit=IIN,file=trim(prname_lp),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error: could not open database file: ',trim(prname_lp)
    call exit_mpi(myrank, 'Error reading neighbors external mesh file')
  endif
    !===========BinHe==========! 
    !copy from src/tomography/post/sem_smooth.F90
    read(IIN) NSPEC_N
    read(IIN) NGLOB_N
    read(IIN) NSPEC_IRREGULAR_N
    ! allocates arrays
    allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1023')
    if (ier /= 0) stop 'Error allocating array ibool'
    allocate(xstore(NGLOB_N),ystore(NGLOB_N),zstore(NGLOB_N),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1024')
    if (ier /= 0) stop 'Error allocating array xstore etc.'

    if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
      if (NSPEC_IRREGULAR_N > 0) then
        allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR_N),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1025')
        if (ier /= 0) stop 'Error allocating array jacobian'
        allocate(dummy(NGLLX,NGLLY,NGLLZ,NSPEC_N),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1026')
        if (ier /= 0) stop 'Error allocating array dummy'
      else
        allocate(jacobian(1,1,1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1027')
        if (ier /= 0) stop 'Error allocating array jacobian'
        allocate(dummy(1,1,1,1),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 1028')
        if (ier /= 0) stop 'Error allocating array dummy'
      endif
      allocate(irregular_element_number(NSPEC_N),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1029')
      if (ier /= 0) stop 'Error allocating array irregular_element_number'
    endif

    ! ibool file
    read(IIN) ibool

    ! global point arrays
    read(IIN) xstore
    read(IIN) ystore
    read(IIN) zstore

    if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
    ! reads in jacobian
      read(IIN) irregular_element_number
      read(IIN) xix_regular
      read(IIN) jacobian_regular
      read(IIN) dummy ! xix
      read(IIN) dummy ! xiy
      read(IIN) dummy ! xiz
      read(IIN) dummy ! etax
      read(IIN) dummy ! etay
      read(IIN) dummy ! etaz
      read(IIN) dummy ! gammax
      read(IIN) dummy ! gammay
      read(IIN) dummy ! gammaz
      read(IIN) jacobian
    endif

  close(IIN)

end subroutine smooth_sem_fwat

!
! -----------------------------------------------------------------------------
!
  subroutine smoothing_weights_vec(x0,y0,z0,sigma_h2_inv,sigma_v2_inv,exp_val, &
                              xx_elem,yy_elem,zz_elem)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: exp_val
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: xx_elem, yy_elem, zz_elem
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,sigma_h2_inv,sigma_v2_inv

  ! local parameters
  real(kind=CUSTOM_REAL) :: dist_h_sq,dist_v_sq
  real(kind=CUSTOM_REAL) :: val
  real(kind=CUSTOM_REAL) :: x1,y1,z1

  integer :: i,j,k

  do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX

    ! point in second slice
    x1 = xx_elem(i,j,k)
    y1 = yy_elem(i,j,k)
    z1 = zz_elem(i,j,k)

    ! gets vertical and horizontal distance
    ! vertical distance (squared)
    dist_v_sq = (z0-z1)*(z0-z1)

    ! horizontal distance (squared)
    dist_h_sq = (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)

    ! Gaussian function
    val = exp(- dist_h_sq * sigma_h2_inv - dist_v_sq * sigma_v2_inv)

    ! stores values in element array
    exp_val(i,j,k) = val

  enddo; enddo; enddo ! NGLLZ,NGLLY,NGLLX

  end subroutine smoothing_weights_vec

!
! -----------------------------------------------------------------------------
!

  subroutine get_distance_vec(dist_h,dist_v,x0,y0,z0,x1,y1,z1)

! returns vector lengths as distances in radial and horizontal direction
! only for flat Earth with z in vertical direction

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: dist_h,dist_v
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,x1,y1,z1

  ! vertical distance
  !dist_v = sqrt( (z0-z1)*(z0-z1) )
  ! squared (to avoid costly square-root
  dist_v = (z0-z1)*(z0-z1)

  ! horizontal distance
  !dist_h = sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) )
  ! squared
  dist_h = (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)

  end subroutine get_distance_vec

!---------------------------------------------------
subroutine smooth_sem_pde_fwat(sigma_h,sigma_v,kernel_name,input_dir,&
                      output_dir,USE_GPU)

  use constants, only: HUGEVAL
  use postprocess_par, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,NGLLSQUARE, &
    MAX_STRING_LEN,IIN,IOUT,GAUSSALPHA,GAUSSBETA,PI,TWO_PI
  use specfem_par
  use specfem_par_elastic, only: nspec_inner_elastic, &
    nspec_outer_elastic,phase_ispec_inner_elastic
  ! use fullwave_adjoint_tomo_par, only: USE_SPH_SMOOTH, KERNEL_TAPER
  use fwat_input, only: tomo_par
  use taper3d
  !use specfem_par_acoustic, only: ispec_is_acoustic
  !use specfem_par_poroelastic, only: ispec_is_poroelastic

  implicit none 
  !integer, parameter :: NARGS = 6
  ! type(taper_cls) :: taper_3d
  integer, parameter :: PRINT_INFO_PER_STEP = 100000
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: rotate_r
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    dx_elem,dy_elem,dz_elem, stemp1,stemp2,stemp3, snewtemp1,snewtemp2,snewtemp3
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dat_glob, ddat_glob
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rvol
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: rvol_local
  integer :: i,j,k,l,iglob,ier,ispec,ispec_p,iphase,ispec_irreg

  !character(len=MAX_STRING_LEN) :: arg(6)
  character(len=MAX_STRING_LEN) :: input_dir, output_dir
  !character(len=MAX_STRING_LEN) :: prname_lp
  character(len=MAX_STRING_LEN*2) :: ks_file
  character(len=MAX_STRING_LEN*2) :: local_data_file


  !character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  !character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: kernel_name
  integer :: num_elements
  real t1,t2,tnow,tlast

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_v, ch, cv, cmax
  real(kind=CUSTOM_REAL) :: min_val, max_val, min_val_glob, max_val_glob
  
  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  real(kind=CUSTOM_REAL) :: xl,yl,zl,rl,rxl,ryl,rzl,&
    xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl 
  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3
  integer :: ntstep, istep
  double precision :: weight
  !real(kind=CUSTOM_REAL), dimension(MAX_KERNEL_NAMES) :: max_old,max_new,max_old_all,max_new_all
  !real(kind=CUSTOM_REAL), dimension(MAX_KERNEL_NAMES) :: min_old,min_new,min_old_all,min_new_all
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_vector_ext_mesh_smooth
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_recv_vector_ext_mesh_smooth
  logical :: BROADCAST_AFTER_READ, USE_GPU
 
  !call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) print *,"Running XSMOOTH_SEM_PDE"
  call synchronize_all()
  call cpu_time(t1)

  !! parse command line arguments
  !if (command_argument_count() /= NARGS) then
  !  if (myrank == 0) then
  !      print *,'USAGE:  mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR GPU_MODE'
  !    stop 'Please check command line arguments'
  !  endif
  !endif
  !call synchronize_all()

  !do i = 1, NARGS
  !  call get_command_argument(i,arg(i), status=ier)
  !enddo

  !read(arg(1),*) sigma_h
  !read(arg(2),*) sigma_v
  !kernel_name = arg(3)
  !input_dir= arg(4)
  !output_dir = arg(5)
  !read(arg(6),*) USE_GPU

  !call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker) 
  call synchronize_all()
  ! user output
  if (myrank == 0) then
    print *,'command line arguments:'
    print *,'  smoothing sigma_h , sigma_v                : ',sigma_h,sigma_v
    print *,'  input dir : ',trim(input_dir)
    print *,'  output dir: ',trim(output_dir)
    print *,"  GPU_MODE: ", USE_GPU
    print *
  endif
  
  ! reads the parameter file
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  if (ADIOS_ENABLED) stop 'Flag ADIOS_ENABLED not supported yet for smoothing, please rerun program...'

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *,'Error number of processors supposed to run on: ',NPROC
      print *,'Error number of MPI processors actually run on: ',sizeprocs
      print *
      print *,'Please rerun with: mpirun -np ',NPROC,' bin/xsmooth_sem .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif

  ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some
  ! array sizes below
  !call read_mesh_for_init()
  !
  !allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),irregular_element_number(NSPEC_AB),stat=ier)
  !if (ier /= 0) call exit_MPI_without_rank('error allocating array 980')

  !if (NSPEC_IRREGULAR > 0) then
  !  ! allocate arrays for storing the databases
  !  allocate(xix(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 981')
  !  allocate(xiy(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 982')
  !  allocate(xiz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 983')
  !  allocate(etax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 984')
  !  allocate(etay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 985')
  !  allocate(etaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 986')
  !  allocate(gammax(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 987')
  !  allocate(gammay(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 988')
  !  allocate(gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 989')
  !  allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 990')
  ! else
  !     ! allocate arrays for storing the databases
  !  allocate(xix(1,1,1,1),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 991')
  !  allocate(xiy(1,1,1,1),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 992')
  !  allocate(xiz(1,1,1,1),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 993')
  !  allocate(etax(1,1,1,1),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 994')
  !  allocate(etay(1,1,1,1),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 995')
  !  allocate(etaz(1,1,1,1),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 996')
  !  allocate(gammax(1,1,1,1),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 997')
  !  allocate(gammay(1,1,1,1),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 998')
  !  allocate(gammaz(1,1,1,1),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 999')
  !  allocate(jacobian(1,1,1,1),stat=ier)
  !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1000')
  !endif
  !if (ier /= 0) stop 'Error allocating arrays for databases'

  !! mesh node locations
  !allocate(xstore(NGLOB_AB),stat=ier)
  !if (ier /= 0) call exit_MPI_without_rank('error allocating array 1001')
  !allocate(ystore(NGLOB_AB),stat=ier)
  !if (ier /= 0) call exit_MPI_without_rank('error allocating array 1002')
  !allocate(zstore(NGLOB_AB),stat=ier)
  !if (ier /= 0) call exit_MPI_without_rank('error allocating array 1003')
  !if (ier /= 0) stop 'Error allocating arrays for mesh nodes'

  !! material properties
  !allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  !if (ier /= 0) call exit_MPI_without_rank('error allocating array 1004')
  !allocate(mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  !if (ier /= 0) call exit_MPI_without_rank('error allocating array 1005')
  !if (ier /= 0) stop 'Error allocating arrays for material properties'

  !! material flags
  !allocate(ispec_is_acoustic(NSPEC_AB),stat=ier)
  !if (ier /= 0) call exit_MPI_without_rank('error allocating array 1006')
  !allocate(ispec_is_elastic(NSPEC_AB),stat=ier)
  !if (ier /= 0) call exit_MPI_without_rank('error allocating array 1007')
  !allocate(ispec_is_poroelastic(NSPEC_AB),stat=ier)
  !if (ier /= 0) call exit_MPI_without_rank('error allocating array 1008')
  !if (ier /= 0) stop 'Error allocating arrays for material flags'
  !ispec_is_acoustic(:) = .false.
  !ispec_is_elastic(:) = .false.
  !ispec_is_poroelastic(:) = .false.


  !! reads in external mesh
  !call read_mesh_databases()

  !! gets mesh dimensions
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                            x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            distance_min_glob,distance_max_glob)

  !! outputs infos
  !if (myrank == 0) then
  !  print *,'mesh dimensions:'
  !  print *,'  Xmin and Xmax of the model = ',x_min_glob,x_max_glob
  !  print *,'  Ymin and Ymax of the model = ',y_min_glob,y_max_glob
  !  print *,'  Zmin and Zmax of the model = ',z_min_glob,z_max_glob
  !  print *
  !  print *,'  Max GLL point distance = ',distance_max_glob
  !  print *,'  Min GLL point distance = ',distance_min_glob
  !  print *,'  Max/min ratio = ',distance_max_glob/distance_min_glob
  !  print *
  !  print *,'  Max element size = ',elemsize_max_glob
  !  print *,'  Min element size = ',elemsize_min_glob
  !  print *,'  Max/min ratio = ',elemsize_max_glob/elemsize_min_glob
  !  print *
  !endif

  !! broadcast distance_min_glob to other processors
  ! check before broadcast
  !if (myrank == 1) print *, 'distance_min_glob = ', distance_min_glob, 'myrank=', myrank
  call bcast_all_singlecr(distance_min_glob)
  ! check after broadcast
  !if (myrank == 1) print *, 'distance_min_glob = ', distance_min_glob, 'myrank=', myrank

  if (tomo_par%USE_SPH_SMOOTH) then
    allocate(rotate_r(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1013')

    do ispec = 1, NSPEC_AB
      do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
        iglob = ibool(i,j,k,ispec)
        xl = xstore(iglob)
        yl = ystore(iglob)
        zl = zstore(iglob)
        rl = sqrt(xl*xl+yl*yl+zl*zl)
        rotate_r(1,i,j,k,ispec) = xl / rl
        rotate_r(2,i,j,k,ispec) = yl / rl
        rotate_r(3,i,j,k,ispec) = zl / rl
      enddo;enddo;enddo
    enddo
  endif
  !deallocate(xstore,ystore,zstore)
  !deallocate(ispec_is_acoustic, ispec_is_elastic, ispec_is_poroelastic)

  !! determine ch, cv, ntstep
  cmax = distance_min_glob ** 2 / 6.0
  if (sigma_v >= sigma_h) then
    cv = cmax
    ch = cv * (sigma_h ** 2) / (sigma_v ** 2)
  else
    ch = cmax
    cv = ch * (sigma_v ** 2) / (sigma_h ** 2)
  endif
  ntstep = int(ceiling((max(sigma_h,sigma_v)**2)/(2.0*cmax)))
  
  if (myrank == 0) print *, 'cv=', cv, 'ch=', ch, 'ntstep=', ntstep

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! initialize time iteration
   ! set up GLL points, weights and derivation matrices for reference element
   ! (between -1,1)
  call define_derivation_matrices(xigll,yigll,zigll,wxgll,wygll,wzgll, &
                                  hprime_xx,hprime_yy,hprime_zz, &
                                  hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                                  wgllwgll_xy,wgllwgll_xz,wgllwgll_yz)
  ! define transpose of derivation matrix
  do j = 1,NGLLY
    do i = 1,NGLLX
      hprime_xxT(j,i) = hprime_xx(i,j)
      hprime_yyT(j,i) = hprime_yy(i,j)
      hprime_zzT(j,i) = hprime_zz(i,j)

      hprimewgll_xxT(j,i) = hprimewgll_xx(i,j)
    enddo
  enddo

  allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  allocate(dat_glob(NGLOB_AB))
  allocate(ddat_glob(NGLOB_AB))
  allocate(buffer_send_vector_ext_mesh_smooth( &
              max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  allocate(buffer_recv_vector_ext_mesh_smooth( &
              max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
   ! prepare assemble array
  allocate(rvol(NGLOB_AB)) 
  rvol(:) = 0.0
  allocate(rvol_local(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
  do ispec = 1, NSPEC_AB
    ispec_irreg = irregular_element_number(ispec)
    do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
      weight =  wxgll(i)*wygll(j)*wzgll(k)
      jacobianl = jacobian(i,j,k,ispec_irreg)
      rvol_local(i,j,k,ispec) = real(dble(jacobianl)*weight,kind=CUSTOM_REAL)
      iglob = ibool(i,j,k,ispec)
      rvol(iglob) = rvol(iglob) + rvol_local(i,j,k,ispec)
    enddo;enddo;enddo
  enddo
  call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rvol, &
                       num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,&
                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                       my_neighbors_ext_mesh)
  rvol(:) = 1.0 / rvol(:)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! read in data to be smoothed
  ! data file
  write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'
  local_data_file = trim(prname) // trim(kernel_name) // '.bin'

  open(unit = IIN,file = trim(local_data_file),status='old',action='read',&
       form ='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening data file: ',trim(local_data_file)
    stop 'Error opening data file'
  endif

  read(IIN) dat
  close(IIN)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! project
  dat_glob(:) = 0.0
  do ispec = 1, NSPEC_AB
    do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
      iglob = ibool(i,j,k,ispec)
      dat_glob(iglob) = dat_glob(iglob) + dat(i,j,k,ispec) &
                               * rvol_local(i,j,k,ispec)
    enddo;enddo;enddo
  enddo
  call assemble_MPI_send_smooth(NPROC,NGLOB_AB,&
          dat_glob,buffer_send_vector_ext_mesh_smooth,&
          buffer_recv_vector_ext_mesh_smooth,&
          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
          my_neighbors_ext_mesh, &
          request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
  call assemble_MPI_w_ord_smooth(NPROC,NGLOB_AB,&
          dat_glob,buffer_recv_vector_ext_mesh_smooth,num_interfaces_ext_mesh,&
          max_nibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
          request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
          my_neighbors_ext_mesh,myrank)
  if (myrank == 0) print *, 'Before smoothing: '

  dat_glob(:) = dat_glob(:) * rvol(:)
  min_val = minval(dat_glob)
  max_val = maxval(dat_glob)
  call min_all_cr(min_val, min_val_glob)
  call max_all_cr(max_val, max_val_glob)
  if (myrank == 0) then
    print *, '  '//trim(kernel_name)
    print *, '    minval:', min_val_glob
    print *, '    maxval:', max_val_glob
    if (myrank == 0) call cpu_time(tlast)
  endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! broadcast glob array back to local array
  do ispec = 1, NSPEC_AB
    do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
      iglob = ibool(i,j,k,ispec)
      dat(i,j,k,ispec) = dat_glob(iglob)
    enddo;enddo;enddo
  enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do istep = 1, ntstep
    ddat_glob(:) = 0.0
    do iphase = 1,2
      if (iphase == 1) then
        num_elements = nspec_outer_elastic
      else
        num_elements = nspec_inner_elastic
      endif
      do ispec_p = 1,num_elements
        ispec = phase_ispec_inner_elastic(ispec_p,iphase)
        call get_gradient_element(dat(:,:,:,ispec), &
          dx_elem,dy_elem,dz_elem,&
          xix(:,:,:,ispec),xiy(:,:,:,ispec),xiz(:,:,:,ispec),&
          etax(:,:,:,ispec),etay(:,:,:,ispec),etaz(:,:,:,ispec),&
          gammax(:,:,:,ispec),gammay(:,:,:,ispec),gammaz(:,:,:,ispec))
        ispec_irreg = irregular_element_number(ispec)
        do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
          xixl = xix(i,j,k,ispec_irreg)
          xiyl = xiy(i,j,k,ispec_irreg)
          xizl = xiz(i,j,k,ispec_irreg)
          etaxl = etax(i,j,k,ispec_irreg)
          etayl = etay(i,j,k,ispec_irreg)
          etazl = etaz(i,j,k,ispec_irreg)
          gammaxl = gammax(i,j,k,ispec_irreg)
          gammayl = gammay(i,j,k,ispec_irreg)
          gammazl = gammaz(i,j,k,ispec_irreg)
          jacobianl = jacobian(i,j,k,ispec_irreg)
          if (tomo_par%USE_SPH_SMOOTH) then
            rxl = rotate_r(1,i,j,k,ispec)
            ryl = rotate_r(2,i,j,k,ispec)
            rzl = rotate_r(3,i,j,k,ispec)
            stemp1(i,j,k) = ((cv-ch) * (rxl*xixl+ryl*xiyl+rzl*xizl) * &
              (rxl*dx_elem(i,j,k)+ryl*dy_elem(i,j,k)+rzl*dz_elem(i,j,k)) +&
              ch * (xixl*dx_elem(i,j,k)+xiyl*dy_elem(i,j,k)&
                  +xizl*dz_elem(i,j,k))) * jacobianl
            stemp2(i,j,k) = ((cv-ch) * (rxl*etaxl+ryl*etayl+rzl*etazl) * &
              (rxl*dx_elem(i,j,k)+ryl*dy_elem(i,j,k)+rzl*dz_elem(i,j,k)) +&
              ch * (etaxl*dx_elem(i,j,k)+etayl*dy_elem(i,j,k)&
                  +etazl*dz_elem(i,j,k))) * jacobianl
            stemp3(i,j,k) = ((cv-ch) * (rxl*gammaxl+ryl*gammayl+rzl*gammazl) * &
              (rxl*dx_elem(i,j,k)+ryl*dy_elem(i,j,k)+rzl*dz_elem(i,j,k)) +&
              ch * (gammaxl*dx_elem(i,j,k)+gammayl*dy_elem(i,j,k)&
                  +gammazl*dz_elem(i,j,k))) * jacobianl
          else
            stemp1(i,j,k) = ((cv-ch) * xizl * dz_elem(i,j,k) +&
              ch * (xixl*dx_elem(i,j,k)+xiyl*dy_elem(i,j,k)&
                  +xizl*dz_elem(i,j,k))) * jacobianl
            stemp2(i,j,k) = ((cv-ch) * etazl * dz_elem(i,j,k) +&
              ch * (etaxl*dx_elem(i,j,k)+etayl*dy_elem(i,j,k)&
                  +etazl*dz_elem(i,j,k))) * jacobianl
            stemp3(i,j,k) = ((cv-ch) * gammazl * dz_elem(i,j,k) +&
              ch * (gammaxl*dx_elem(i,j,k)+gammayl*dy_elem(i,j,k)&
                  +gammazl*dz_elem(i,j,k))) * jacobianl
          endif
        enddo;enddo;enddo
        do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
          snewtemp1(i,j,k) = 0.0
          snewtemp2(i,j,k) = 0.0
          snewtemp3(i,j,k) = 0.0
          do l = 1, NGLLX
            fac1 = hprimewgll_xx(l,i)
            snewtemp1(i,j,k) = snewtemp1(i,j,k) + stemp1(l,j,k) * fac1
            fac2 = hprimewgll_yy(l,j)
            snewtemp2(i,j,k) = snewtemp2(i,j,k) + stemp2(i,l,k) * fac2
            fac3 = hprimewgll_zz(l,k)
            snewtemp3(i,j,k) = snewtemp3(i,j,k) + stemp3(i,j,l) * fac3
          enddo
          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)
          iglob = ibool(i,j,k,ispec)
          ddat_glob(iglob) = ddat_glob(iglob) - (fac1*snewtemp1(i,j,k)+&
                          fac2 * snewtemp2(i,j,k) + fac3 * snewtemp3(i,j,k))
        enddo;enddo;enddo 
      enddo  ! ispec_p = 1, num_elements
      !! assemble MPI
      if (iphase == 1) then
        call assemble_MPI_send_smooth(NPROC,NGLOB_AB,&
          ddat_glob,buffer_send_vector_ext_mesh_smooth,&
          buffer_recv_vector_ext_mesh_smooth,&
          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
          my_neighbors_ext_mesh, &
          request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
      else
        call assemble_MPI_w_ord_smooth(NPROC,NGLOB_AB,&
          ddat_glob,buffer_recv_vector_ext_mesh_smooth,num_interfaces_ext_mesh,&
          max_nibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
          request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
          my_neighbors_ext_mesh,myrank)
      endif
      !!!!!!!!!!!!!!!!!
    enddo !iphase = 1,2

    ddat_glob(:) = ddat_glob(:) * rvol(:)
    !! update
    dat_glob(:) = dat_glob(:) + ddat_glob(:)
    !! info
    if (mod(istep, PRINT_INFO_PER_STEP) == 0) then
      if (myrank == 0) print *, 'Step:', istep
      min_val = minval(dat_glob)
      max_val = maxval(dat_glob)
      call min_all_cr(min_val, min_val_glob)
      call max_all_cr(max_val, max_val_glob)
      if (myrank == 0) then
        print *, '  '//trim(kernel_name)
        print *, '    minval:', min_val_glob
        print *, '    maxval:', max_val_glob
        call cpu_time(tnow)
        print *, 'time since last message:', tnow-tlast
        call cpu_time(tlast)
      endif
    endif
    !!!!!!!!!!!!!
    !! broadcast glob array back to local array
    ! add taper here
    ! do ispec = 1, NSPEC_AB
    !   do k=1,NGLLZ;do j=1,NGLLY;do i=1,NGLLX
    !     iglob = ibool(i,j,k,ispec)
    !     if (tomo_par%KERNEL_TAPER) then
    !       taper_val = taper_3d%interp(xstore(iglob),ystore(iglob),zstore(iglob))
    !     else
    !       taper_val = 1.
    !     endif
    !     dat(i,j,k,ispec) = taper_val*dat_glob(iglob)
    !   enddo;enddo;enddo
    ! enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call synchronize_all()
  enddo
  call synchronize_all()
  call cpu_time(t2)
  if (myrank == 0) & 
    print *, 'Computation time with PDE-based smoothing on CPU:', t2-t1
  !! output
  ! file output
  ! smoothed kernel file name
  ! statistics
  min_val = minval(dat_glob)
  max_val = maxval(dat_glob)
  call min_all_cr(min_val, min_val_glob)
  call max_all_cr(max_val, max_val_glob)
  if (myrank == 0) then
    print *, 'After smoothing:'
    print *, '  '//trim(kernel_name)
    print *, '    minval:', min_val_glob
    print *, '    maxval:', max_val_glob
  endif
  write(ks_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_name)//'_smooth.bin'
  open(IOUT,file=trim(ks_file),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error opening smoothed kernel file'
  write(IOUT) dat
  close(IOUT)
  if (myrank == 0) print *,'written: ',trim(ks_file)

  !deallocate(ibool,irregular_element_number)
  !deallocate(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian)
  if (tomo_par%USE_SPH_SMOOTH) deallocate(rotate_r)
  deallocate(dat, dat_glob, ddat_glob)
  deallocate(buffer_send_vector_ext_mesh_smooth, &
             buffer_recv_vector_ext_mesh_smooth)
  deallocate(rvol, rvol_local)
  

  !call finalize_mpi()

end subroutine smooth_sem_pde_fwat

  subroutine assemble_MPI_send_smooth(NPROC,NGLOB_AB,&
          array_val,buffer_send_vector_ext_mesh_smooth,&
          buffer_recv_vector_ext_mesh_smooth,&
          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
          my_neighbors_ext_mesh, &
          request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

    ! sends data

  use constants, only: CUSTOM_REAL, itag

  implicit none

  integer :: NPROC
  integer :: NGLOB_AB

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), &
    dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_vector_ext_mesh_smooth,buffer_recv_vector_ext_mesh_smooth

  integer, dimension(num_interfaces_ext_mesh) :: &
    nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh):: &
    ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: &
    request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer ipoin,iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        buffer_send_vector_ext_mesh_smooth(ipoin,iinterface) = &
             array_val(ibool_interfaces_ext_mesh(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      call isend_cr(buffer_send_vector_ext_mesh_smooth(1,iinterface), &
                    nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_send_vector_ext_mesh(iinterface))
      call irecv_cr(buffer_recv_vector_ext_mesh_smooth(1,iinterface), &
                    nibool_interfaces_ext_mesh(iinterface), &
                    my_neighbors_ext_mesh(iinterface), &
                    itag, &
                    request_recv_vector_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_send_smooth


  subroutine assemble_MPI_w_ord_smooth(NPROC,NGLOB_AB, &
          array_val,buffer_recv_vector_ext_mesh_smooth,num_interfaces_ext_mesh,&
          max_nibool_interfaces_ext_mesh, &
          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
          request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
          my_neighbors_ext_mesh,myrank)

! waits for data to receive and assembles

! The goal of this version is to avoid different round-off errors in different
! processors.
! The contribution of each processor is added following the order of its rank.
! This guarantees that the sums are done in the same order on all processors.
!
! NOTE: this version assumes that the interfaces are ordered by increasing rank
! of the neighbor.
! That is currently done so in subroutine write_interfaces_database in
! decompose_mesh_SCOTCH/part_decompose_mesh_SCOTCH.f90
! A safety test could be added here.
!
! October 2012 - Surendra Somala and Jean-Paul Ampuero - Caltech Seismolab

  use constants, only: CUSTOM_REAL, itag

  implicit none

  integer :: NPROC
  integer :: NGLOB_AB
! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,myrank

  real(kind=CUSTOM_REAL), &
    dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh_smooth

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh)::&
    ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: &
    request_send_vector_ext_mesh,request_recv_vector_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbors_ext_mesh

  real(kind=CUSTOM_REAL), &
    dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
    mybuffer
  integer :: ipoin,iinterface,iglob
  logical :: need_add_my_contrib

! here we have to assemble all the contributions between partitions using MPI

! assemble only if more than one partition
  if (NPROC == 1) return

! move interface values of array_val to local buffers
  do iinterface = 1, num_interfaces_ext_mesh
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
      mybuffer(ipoin,iinterface) = array_val(iglob)
     ! set them to zero right away to avoid counting it more than once during
     ! assembly:
     ! buffers of higher rank get zeros on nodes shared with current buffer
      array_val(iglob) = 0._CUSTOM_REAL
    enddo
  enddo

! wait for communications completion (recv)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_recv_vector_ext_mesh(iinterface))
  enddo

! adding all contributions in order of processor rank
  need_add_my_contrib = .true.
  do iinterface = 1, num_interfaces_ext_mesh
    if (need_add_my_contrib .and. myrank < my_neighbors_ext_mesh(iinterface)) &
      call add_my_contrib()
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
      array_val(iglob) = array_val(iglob) + &
        buffer_recv_vector_ext_mesh_smooth(ipoin,iinterface)
    enddo
  enddo
  if (need_add_my_contrib) call add_my_contrib()

! wait for communications completion (send)
  do iinterface = 1, num_interfaces_ext_mesh
    call wait_req(request_send_vector_ext_mesh(iinterface))
  enddo

  contains

    subroutine add_my_contrib()

    integer :: my_iinterface,my_ipoin

    do my_iinterface = 1, num_interfaces_ext_mesh
      do my_ipoin = 1, nibool_interfaces_ext_mesh(my_iinterface)
        iglob = ibool_interfaces_ext_mesh(my_ipoin,my_iinterface)
        array_val(iglob) = array_val(iglob) + &
          mybuffer(my_ipoin,my_iinterface)
      enddo
    enddo
    need_add_my_contrib = .false.

    end subroutine add_my_contrib

  end subroutine assemble_MPI_w_ord_smooth


  subroutine get_gradient_element(s, dx_elem, dy_elem, dz_elem, &
     xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
  use constants, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ
  use specfem_par, only: hprime_xxT,hprime_yyT,hprime_zzT
  implicit none
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ), intent(in) :: s
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ), intent(out) :: &
          dx_elem, dy_elem, dz_elem
  integer :: i,j,k,l
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
  real(kind=CUSTOM_REAL) :: temp1l, temp2l, temp3l
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xix,xiy,xiz,etax,&
          etay,etaz,gammax,gammay,gammaz
  dx_elem(:,:,:) = 0.0
  dy_elem(:,:,:) = 0.0
  dz_elem(:,:,:) = 0.0
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        temp1l = 0.0
        temp2l = 0.0
        temp3l = 0.0
        do l = 1, NGLLX
          hp1 = hprime_xxT(l,i)
          temp1l = temp1l + s(l,j,k) * hp1
        enddo
        do l = 1, NGLLY
          hp2 = hprime_yyT(l,j)
          temp2l = temp2l + s(i,l,k) * hp2
        enddo
        do l = 1, NGLLZ
          hp3 = hprime_zzT(l,k)
          temp3l = temp3l + s(i,j,l) * hp3
        enddo
        dx_elem(i,j,k)=temp1l*xix(i,j,k)+&
                temp2l*etax(i,j,k)+temp3l*gammax(i,j,k)
        dy_elem(i,j,k)=temp1l*xiy(i,j,k)+&
                temp2l*etay(i,j,k)+temp3l*gammay(i,j,k)
        dz_elem(i,j,k)=temp1l*xiz(i,j,k)+&
                temp2l*etaz(i,j,k)+temp3l*gammaz(i,j,k)
      enddo
    enddo
  enddo
  end subroutine get_gradient_element
