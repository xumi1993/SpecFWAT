module postproc_sub
  use fullwave_adjoint_tomo_par
  use fwat_input
  use fwat_utils
  !! for reading database
  use specfem_par
  use taper3d
  use tomography_par, only: USE_ALPHA_BETA_RHO,USE_ISO_KERNELS, RHO_SCALING, &
                NGLOB, NSPEC,THRESHOLD_HESS
  
  implicit none


  character(len=MAX_STRING_LEN), public, dimension(2)  :: set_range
  character(len=MAX_STRING_LEN), public                :: type_name, sum_dir, model, model_prev,model_next 
  character(len=MAX_STRING_LEN), public                :: ekernel_dir_list='optimize/ekernel_dir.lst'
  real(kind=CUSTOM_REAL), public                       :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL), public                       :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL), public                       :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL), public                       :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL), public                       :: z_min_glob,z_max_glob
  type(taper_cls), public                              :: taper_3d
  character(len=MAX_STRING_LEN), public, dimension(:), allocatable :: fwat_kernel_names
  integer, public                                      :: kernel_num, isetb, isete, imod_current, imod_up,imod_down


  contains

  subroutine get_model_idx()
    read(model(2:),'(I2.2)') imod_current
    imod_up=imod_current+1
    imod_down=imod_current-1
    write(model_next,'(A1,I2.2)') 'M',imod_up
    if (imod_down < tomo_par%ITER_START) then
      model_prev='none'
    else
      write(model_prev,'(A1,I2.2)') 'M',imod_down
    endif
  end subroutine get_model_idx

  subroutine select_set_range()
    if (type_name=='noise') then
      set_range = tomo_par%NOISE_SET_RANGE(:)
      PRECOND_TYPE = noise_par%PRECOND_TYPE
    elseif (index(type_name,'tele')/=0 .or. type_name=='rf') then
      set_range = tomo_par%TELE_SET_RANGE(:)
      PRECOND_TYPE = tele_par%PRECOND_TYPE
    elseif (type_name=='leq') then
      set_range = tomo_par%LEQ_SET_RANGE(:)
      PRECOND_TYPE = leq_par%PRECOND_TYPE
    endif
    read(set_range(1)(4:),'(I3)') isetb 
    read(set_range(2)(4:),'(I3)') isete
    if (count(tomo_par%inv_type) > 1) then
      sum_dir = 'optimize/SUM_KERNELS_'//trim(model)//'_'//trim(type_name)
    else
      sum_dir = 'optimize/SUM_KERNELS_'//trim(model)
    endif

  end subroutine select_set_range

  subroutine set_evt_kernel_path()
    implicit none

    character(MAX_STRING_LEN)                          :: strset
    integer                                            :: ier,iset

    open(unit=400,file=trim(ekernel_dir_list),iostat=ier)
    do iset=isetb,isete
      write(strset,'(I0)') iset
      write(400,'(a)')'solver/'//trim(model)//'.set'//trim(strset)//'/GRADIENT'
    enddo
    close(400)
  end subroutine set_evt_kernel_path

  subroutine get_kernel_names()
    implicit none

    if (USE_ISO_KERNELS) then
      kernel_num = 3
      if (.not. allocated(fwat_kernel_names))allocate(fwat_kernel_names(kernel_num))
      fwat_kernel_names(1) = 'bulk_c_kernel'
      fwat_kernel_names(2) = 'bulk_beta_kernel'
      fwat_kernel_names(3) = 'rho_kernel'
    elseif (USE_ALPHA_BETA_RHO) then
      kernel_num = 3
      if (.not. allocated(fwat_kernel_names))allocate(fwat_kernel_names(kernel_num))
      fwat_kernel_names(1) = 'alpha_kernel'
      fwat_kernel_names(2) = 'beta_kernel'
      fwat_kernel_names(3) = 'rhop_kernel'
    else
      kernel_num = 4
      if (.not. allocated(fwat_kernel_names))allocate(fwat_kernel_names(kernel_num))
      fwat_kernel_names(1) = 'bulk_c_kernel'
      fwat_kernel_names(2) = 'bulk_betav_kernel'
      fwat_kernel_names(3) = 'bulk_betah_kernel'
      fwat_kernel_names(4) = 'eta_kernel'
    endif
  end subroutine get_kernel_names

  subroutine post_proc()
    use my_mpi

    implicit none

    ! integer, intent(in)                                :: inv_type
    ! character(MAX_STRING_LEN)                          :: model_name
    logical                                            :: norm_hess

    call select_set_range()
    if (myrank==0) then
      call set_evt_kernel_path()
      call system('mkdir -p '//trim(sum_dir))
      call system('mkdir -p optimize/SUM_KERNELS_'//trim(model))
    endif

    if(myrank==0) then
      ! write(OUT_FWAT_LOG,*) 'This is sum_precond ...'
      call write_timestamp_log(OUT_FWAT_LOG, 'This is sum_precond ...')
      if(PRECOND_TYPE=='default') then
        write(OUT_FWAT_LOG,*) 'USE THRESHOLD_HESS: ',THRESHOLD_HESS
      else
        write(OUT_FWAT_LOG,*) 'NO THRESHOLD_HESS'
      endif
      flush(OUT_FWAT_LOG)
    endif
    if (count(tomo_par%INV_TYPE)>1) then
      norm_hess = .true.
    else
      norm_hess = .false.
    endif
    call sum_preconditioned_kernels_fwat(ekernel_dir_list,sum_dir, norm_hess)

    if(myrank==0) then
      ! write(OUT_FWAT_LOG,*) 'This is smooth_kernel ...'
      call write_timestamp_log(OUT_FWAT_LOG, 'This is smooth_kernel ...')
      flush(OUT_FWAT_LOG)
    endif
    call smooth_kernel()
    call kernel_taper(sum_dir)
    if (USE_H5) then
      call grid_kernel()
    endif
  end subroutine post_proc

  subroutine grid_kernel()
    use hdf5_interface

    real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: kernel_data
    integer :: i

    do i = 1, kernel_num
      call rg%griddata(sum_dir, fwat_kernel_names(i), kernel_data)
      call h5write(trim(sum_dir)//trim(fwat_kernel_names(i))//'.h5', '/'//trim(fwat_kernel_names(i)), kernel_data)
    enddo
  end subroutine grid_kernel

  subroutine read_database()
    use specfem_par_elastic, only: ELASTIC_SIMULATION,ispec_is_elastic,rho_vp,rho_vs,min_resolved_period
    use specfem_par_acoustic, only: ACOUSTIC_SIMULATION,ispec_is_acoustic
    use specfem_par_poroelastic, only: POROELASTIC_SIMULATION,ispec_is_poroelastic,rho_vpI,rho_vpII,rho_vsI, &
                                      phistore,tortstore,rhoarraystore
    implicit none

    integer :: ier
  ! smooth
    ! Kai: Maybe we should move reding Database back to the smoothing subroutine.
    !!!  Read Database before smoothing !!! 
    ! read the value of NSPEC_AB and NGLOB_AB because we need it to define some array sizes below
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
    allocate(xstore(NGLOB_AB), &
            ystore(NGLOB_AB), &
            zstore(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays for mesh nodes'

    ! ! material properties
    ! allocate(kappastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
    !         mustore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    ! if (ier /= 0) stop 'Error allocating arrays for material properties'

    ! material flags
    allocate(ispec_is_acoustic(NSPEC_AB), &
            ispec_is_elastic(NSPEC_AB), &
            ispec_is_poroelastic(NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays for material flags'
    ispec_is_acoustic(:) = .false.
    ispec_is_elastic(:) = .false.
    ispec_is_poroelastic(:) = .false.

    NSPEC = NSPEC_AB
    NGLOB = NGLOB_AB

    ! reads in external mesh
    call read_mesh_databases()

    ! gets mesh dimensions
    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                              x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                              elemsize_min_glob,elemsize_max_glob, &
                              distance_min_glob,distance_max_glob)
    call synchronize_all()

  end subroutine read_database

  subroutine kernel_taper(input_dir)
    use interpolation_mod, only: trilin_interp
    implicit none

    character(len=MAX_STRING_LEN)                      :: input_dir
    integer                                            :: i, j, k, ispec, n,iglob
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_data, taper_val

    if (tomo_par%KERNEL_TAPER) then
      if (myrank==0) then
        ! write(OUT_FWAT_LOG,*) 'Tapering boudary area for smoothed kernel...'
        call write_timestamp_log(OUT_FWAT_LOG, 'Tapering boudary area for smoothed kernel...')
        flush(OUT_FWAT_LOG)
      endif

      ! gets mesh dimensions
      call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                              x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                              elemsize_min_glob,elemsize_max_glob, &
                              distance_min_glob,distance_max_glob)
      call synchronize_all()
      ! broadcast SEM size 
      call bcast_all_singlecr(x_min_glob)
      call bcast_all_singlecr(x_max_glob)
      call bcast_all_singlecr(y_min_glob)
      call bcast_all_singlecr(y_max_glob)
      call bcast_all_singlecr(z_min_glob)
      call bcast_all_singlecr(z_max_glob)
      call bcast_all_singlecr(x_min_glob)
      call bcast_all_singlecr(elemsize_min_glob)

      ! Initialize taper class
      call taper_3d%Destory()
      ! create regular grids for tapering mask 
      call taper_3d%create(x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                         elemsize_min_glob/2., elemsize_min_glob/2, elemsize_min_glob/2, &
                         tomo_par%TAPER_H_SUPPRESS, tomo_par%TAPER_H_BUFFER,&
                         tomo_par%TAPER_V_SUPPRESS, tomo_par%TAPER_V_BUFFER)
       
      allocate(taper_val(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
      do n = 1, kernel_num
        call read_smoothed_kernel(input_dir, fwat_kernel_names(n), kernel_data)
        do ispec = 1, NSPEC_AB
          do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            taper_val(i,j,k,ispec) = taper_3d%interp(xstore(iglob),ystore(iglob),zstore(iglob))
            kernel_data(i,j,k,ispec) = kernel_data(i,j,k,ispec) * taper_val(i,j,k,ispec)
          enddo;enddo;enddo
        enddo
        call write_smoothed_kernel(input_dir, fwat_kernel_names(n), kernel_data)
      enddo ! i = 1, kernel_num
      deallocate(taper_val)
    endif
  end subroutine kernel_taper

  subroutine kernel_norm(kernel_data)
    implicit none
    integer :: i, j, k, ispec
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable, intent(inout) :: kernel_data
    real(kind=CUSTOM_REAL) :: max_new,min_new,min_new_all,max_new_all,&
                              max_abs(2),max_abs_all
  
    ! statistics
    min_new = minval(kernel_data(:,:,:,:))
    max_new = maxval(kernel_data(:,:,:,:))
    call synchronize_all()
    call min_all_cr(min_new, min_new_all)
    call max_all_cr(max_new,max_new_all)
    call bcast_all_singlecr(min_new_all)
    call bcast_all_singlecr(max_new_all)

    ! Mijian added for normalization
    max_abs(1) = abs(min_new_all)
    max_abs(2) = abs(max_new_all)
    max_abs_all = maxval(max_abs)
    do ispec = 1, NSPEC_AB
      do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX
        kernel_data(i,j,k,ispec) = kernel_data(i,j,k,ispec)/max_abs_all
      enddo;enddo;enddo
    enddo
    call synchronize_all()
  end subroutine kernel_norm

  subroutine smooth_kernel()
    implicit none
    character(len=MAX_STRING_LEN)                      :: input_dir,output_dir
    logical                                            :: USE_GPU
    integer                                            :: i
    real(kind=CUSTOM_REAL)                             :: sigma_h, sigma_v

    if (type_name == 'noise') then
      sigma_h = tomo_par%NOISE_SIGMA_H
      sigma_v = tomo_par%NOISE_SIGMA_V
    elseif (index(type_name, 'tele') /= 0 .or. type_name=='rf') then
      sigma_h = tomo_par%TELE_SIGMA_H
      sigma_v = tomo_par%TELE_SIGMA_V
    elseif (type_name == 'leq') then
      sigma_h = tomo_par%LEQ_SIGMA_H
      sigma_v = tomo_par%LEQ_SIGMA_V
    else
      call exit_MPI('Please specify correct type_name for smooth_kernel')
    endif

    if (tomo_par%is_smooth) then
      USE_GPU = .false.
      input_dir = sum_dir
      output_dir = sum_dir
      do i=1,kernel_num
        if (myrank==0) then
          call write_timestamp_log(OUT_FWAT_LOG, 'Smoothing '//trim(fwat_kernel_names(i)))
          flush(OUT_FWAT_LOG)
        endif
        if (tomo_par%USE_SPH_SMOOTH) then
          call smooth_sem_pde_fwat(sigma_h,sigma_v,fwat_kernel_names(i),input_dir,&
                        output_dir,USE_GPU)
        else
          if (type_name == 'noise' .and. tomo_par%USE_RHO_SCALING_NOISE & 
              .and. fwat_kernel_names(i) == 'rhop_kernel') then
            call rho_scaling_fwat(input_dir,output_dir)
          else
            call smooth_sem_fwat(sigma_h,sigma_v,fwat_kernel_names(i),input_dir,&
                      output_dir,USE_GPU)
          endif
        endif
      enddo
    endif ! end if is_smooth
  end subroutine smooth_kernel

  subroutine sum_joint_kernels()
    ! use postprocess_par, only: FIOUT
    implicit none

    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel1, total_kernel
    integer :: ier,i,j,nchan
    double precision :: norm_val(NUM_INV_TYPE), norm, norm_sum
    character(len=MAX_STRING_LEN) :: output_dir

    allocate(total_kernel(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (myrank==0) then
      call write_timestamp_log(OUT_FWAT_LOG, 'This is summing kernels for joint inversion...')
      flush(OUT_FWAT_LOG)
    endif
    do i = 1, NUM_INV_TYPE
      if (.not. tomo_par%INV_TYPE(i)) cycle
      call calc_kernel0_std_weight(i, norm_val(i))
    enddo
    do i=1,kernel_num
      total_kernel = 0._CUSTOM_REAL
      do j=1,NUM_INV_TYPE
        if (.not. tomo_par%INV_TYPE(j)) cycle
        type_name = tomo_par%INV_TYPE_NAME(j)
        call select_set_range()
        output_dir = 'optimize/SUM_KERNELS_'//trim(model)//'_'//trim(type_name)
        if (USE_H5) then
          call rg%semdata(output_dir, fwat_kernel_names(i), kernel1)        
        else
          call read_smoothed_kernel(output_dir, fwat_kernel_names(i), kernel1)
        endif
        kernel1 = kernel1/norm_val(j)
        norm = maxval( abs(kernel1))
        call max_all_dp(norm, norm_sum)
        if (myrank == 0) then
          print *,'norm ',trim(fwat_kernel_names(i)),' kernel of ',trim(type_name),': ',norm_sum
          print *, 'max of gradient ',norm_val(j)
        endif
        total_kernel = total_kernel + tomo_par%JOINT_WEIGHT(j)*kernel1
        output_dir = 'optimize/SUM_KERNELS_'//trim(model)
        call write_smoothed_kernel(output_dir, fwat_kernel_names(i), total_kernel)
      enddo 
    enddo
    deallocate(kernel1,total_kernel)
  end subroutine sum_joint_kernels

  subroutine calc_kernel0_std_weight(itype, max_global)
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_data
    double precision :: max_local, max_global
    double precision :: max_ker(kernel_num)
    integer :: n, itype, iker, i
    character(len=MAX_STRING_LEN) :: type_name, output_dir, model0

    write(model0, '("M",i2.2)') tomo_par%ITER_START
    type_name = tomo_par%INV_TYPE_NAME(itype)
    output_dir = 'optimize/SUM_KERNELS_'//trim(model0)//'_'//trim(type_name)
    do iker=1,kernel_num
      call read_smoothed_kernel(output_dir, fwat_kernel_names(iker), kernel_data)
      max_local = maxval( abs(kernel_data))
      call max_all_dp(max_local, max_ker(iker))
      call bcast_all_singledp(max_ker(iker))
    enddo
    max_global = maxval(max_ker(:))
    call synchronize_all()
  end subroutine calc_kernel0_std_weight

  subroutine read_smoothed_kernel(input_dir, kernel_name, kernel_data, issmooth)

    logical, optional, intent(in) :: issmooth
    character(len=MAX_STRING_LEN) :: k_file,input_dir,kernel_name,this_kernel_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_data
    integer :: ier

    if(.not. allocated(kernel_data)) allocate(kernel_data(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    kernel_data = 0._CUSTOM_REAL

    if (.not. present(issmooth)) then
      this_kernel_name = trim(kernel_name)//'_smooth'
    elseif (issmooth) then
      this_kernel_name = trim(kernel_name)//'_smooth'
    else
      this_kernel_name = trim(kernel_name)
    endif

    write(k_file,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'//trim(this_kernel_name)//'.bin'
    open(FIIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ier)
    if (ier /= 0) then
      write(*,*) '  kernel not found: ',trim(k_file)
      stop 'Error kernel file not found'
    endif
    read(FIIN) kernel_data
    close(FIIN)
    call synchronize_all()
  end subroutine read_smoothed_kernel

  subroutine write_smoothed_kernel(output_dir, kernel_name, kernel_data, issmooth)

    logical, optional, intent(in) :: issmooth
    character(len=MAX_STRING_LEN) :: k_file,output_dir,kernel_name, this_kernel_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_data
    integer :: ier

    if (.not. present(issmooth)) then
      this_kernel_name = trim(kernel_name)//'_smooth'
    elseif (issmooth) then
      this_kernel_name = trim(kernel_name)//'_smooth'
    else
      this_kernel_name = trim(kernel_name)
    endif

    write(k_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(this_kernel_name)//'.bin'
    open(FIOUT,file=trim(k_file),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error writing smoothed kernel file'
    write(FIOUT) kernel_data(:,:,:,:)
    close(FIOUT)
    call synchronize_all()
  end subroutine write_smoothed_kernel

  subroutine rho_scaling_fwat(input_dir,output_dir)

    character(len=MAX_STRING_LEN)                          :: input_dir,output_dir
    character(len=MAX_STRING_LEN), parameter               :: beta_kernel_name='beta_kernel'
    character(len=MAX_STRING_LEN), parameter               :: rhop_kernel_name='rhop_kernel'
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: rho_kernel, beta_kernel
    integer                                                :: ier

    if (myrank == 0) then
      ! write(OUT_FWAT_LOG, *) 'Rho kernel scaling...'
      call write_timestamp_log(OUT_FWAT_LOG, 'Rho kernel scaling...')
      flush(OUT_FWAT_LOG)
    endif

    allocate(rho_kernel(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    rho_kernel = 0.

    call read_smoothed_kernel(input_dir, beta_kernel_name, beta_kernel)

    rho_kernel = RHO_SCALING * beta_kernel
    call write_smoothed_kernel(output_dir, rhop_kernel_name, rho_kernel)
    deallocate(rho_kernel, beta_kernel)
  end subroutine rho_scaling_fwat

  subroutine read_misfit(simu_type, model_name, sum_chi, nchan)
    implicit none

    integer, intent(out)                               :: nchan
    double precision, intent(out)                      :: sum_chi
    type(chi_table)                                    :: chi
    integer                                            :: setb,sete,iset,nflt,&
                                                          iband,i,nevt
    integer, parameter                                 :: col = 29, max_ndim=100000
    character(len=MAX_STRING_LEN), dimension(:), allocatable  :: evtnames
    character(len=MAX_STRING_LEN)                      :: fname,strset,bandname
    character(len=*),intent(in)                        :: simu_type, model_name
    double precision, dimension(:), allocatable        :: array
    double precision, dimension(max_ndim)              :: chi_non_zero
    real(kind=CUSTOM_REAL), dimension(:), allocatable  :: weight
    double precision                                   :: mean_chi,std_chi
    logical                                            :: exist

    call setup_common_variables(simu_type)
    type_name = simu_type
    call select_set_range()
    read(set_range(1)(4:),'(I3)') setb 
    read(set_range(2)(4:),'(I3)') sete
    if (simu_type == 'rf') then
      nflt = rf_par%NGAUSS
    else
      nflt = NUM_FILTER
    endif
    nchan = 0
    chi_non_zero = 0
    if (myrank == 0) then
      do iset=setb,sete
        write(strset,'("set",I0)') iset
        call read_evt_weight(strset, nevt, evtnames, weight)
        do iband=1, nflt
          if (simu_type == 'rf') then
            write(bandname,'(a1,f3.1)') 'F',rf_par%f0(iband)
          else
            call get_band_name(SHORT_P(iband),LONG_P(iband),bandname)
            ! write(bandname,'(a1,i3.3,a2,i3.3)') 'T',int(SHORT_P(iband)),'_T',int(LONG_P(iband))
          endif
          fname = 'misfits/'//trim(model_name)//'.'//trim(strset)//&
                  '_'//trim(bandname)//'_window_chi'
          inquire(file=trim(fname), exist=exist)
          if (.not. exist) cycle
          call chi%read_misfits(fname)
          allocate(array(chi%n_rows))
          ! fname = 'src_rec/sources'//trim(strset)//'.dat'
          ! call chi%apply_weight(col, nevt, evtnames, weight,array)
          ! print *, chi%tr_chi
          array = chi%get_column(col)
          do i = 1, chi%n_rows
            if (array(i) > 1.d-15) then
              chi_non_zero(nchan+1) = array(i)
              nchan = nchan + 1
            endif
          enddo
          deallocate(array)
          call chi%distory()
        enddo
        deallocate(evtnames, weight)
      enddo
      sum_chi = sum( chi_non_zero(1:nchan) )
      mean_chi = sum_chi / nchan
      std_chi = dsqrt( (sum( chi_non_zero(1:nchan)**2)/nchan - mean_chi**2 )*nchan)
    endif
    call bcast_all_singledp(sum_chi)
    call bcast_all_singlei(nchan)
  end subroutine

end module