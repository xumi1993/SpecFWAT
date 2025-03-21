module post_processing
  use config
  use fwat_mpi
  use fwat_constants
  use logger, only: log
  use input_params, fpar => fwat_par_global
  use utils, only: zeros
  use kernel_io

  implicit none

  type PostFlow
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: ker_names
    real(kind=cr), dimension(:,:,:,:,:), allocatable :: ker_data
    integer :: nker
    character(len=MAX_STRING_LEN), dimension(2) :: simu_types = [SIMU_TYPE_NOISE, SIMU_TYPE_TELE]
    contains
    procedure :: sum_kernel, init, sum_precond, write, generate_for_type, smooth_kernel, sum_joint_kernel
    procedure, private :: calc_kernel0_std_weight
  end type PostFlow
contains

  subroutine init(this)
    class(PostFlow), intent(inout) :: this
    integer :: itype

    call read_mesh_databases_minimum()

    if (fpar%update%model_type == 1) then
      this%nker = size(KERNEL_ISO)
      this%ker_names = KERNEL_ISO
    elseif (fpar%update%model_type == 2) then
      this%nker = size(KERNEL_AZI_ANI)
      this%ker_names = KERNEL_AZI_ANI
    else
      call exit_MPI(0, 'Unknown model type')
    endif

    call system('mkdir -p '//trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name))

    call log%init('output_fwd_post_processing.log')

  end subroutine init

  subroutine generate_for_type(this)
    class(PostFlow), intent(inout) :: this
    if (worldrank == 0) then
      if (is_joint) then
        call system('mkdir -p '//trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)//'_'//trim(simu_type))
      endif
    endif
    call synchronize_all()

    this%ker_data = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_AB, this%nker)
  end subroutine generate_for_type
  
  subroutine sum_kernel(this)
    class(PostFlow), intent(inout) :: this
    real(kind=cr), dimension(:,:,:,:), allocatable :: ker
    integer :: iker, ievt

    call log%write('This is sum kernel...', .true.)
    do iker = 1, this%nker
      do ievt = 1, fpar%acqui%nevents
        call read_event_kernel(ievt, trim(this%ker_names(iker))//'_kernel', ker)
        this%ker_data(:,:,:,:,iker) = this%ker_data(:,:,:,:,iker) + ker
      enddo
    enddo

  end subroutine sum_kernel

  subroutine sum_precond(this)
    class(PostFlow), intent(inout) :: this
    character(len=MAX_STRING_LEN), parameter :: hess_name = 'hess'
    real(kind=cr), dimension(:,:,:,:), allocatable :: total_hess, ker
    real(kind=cr) :: zdep
    integer :: iker, ievt, igllx, iglly, igllz, ispec, iglob

    if (simu_type /= SIMU_TYPE_TELE) then
      total_hess = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_AB)
      do ievt = 1, fpar%acqui%nevents
        call read_event_kernel(ievt, trim(hess_name)//'_kernel', ker)
        total_hess = total_hess + ker
      enddo
      call invert_hess(total_hess)

      do iker = 1, this%nker
        this%ker_data(:,:,:,:,iker) = this%ker_data(:,:,:,:,iker) * total_hess
      enddo
    else
      do ispec = 1, NSPEC_AB
        do igllx = 1, NGLLX
          do iglly = 1, NGLLY
            do igllz = 1, NGLLZ
              iglob = ibool(igllx,iglly,igllz,ispec)
              if (tele_par%PRECOND_TYPE == Z_PRECOND) then
                zdep = abs(zstore(iglob))/abs(z_min_glob)
              elseif (tele_par%PRECOND_TYPE == Z_SQRT_PRECOND) then
                zdep = sqrt(abs(zstore(iglob))/abs(z_min_glob))
              endif
              this%ker_data(igllx,iglly,igllz,ispec,:) = this%ker_data(igllx,iglly,igllz,ispec,:) / zdep
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine sum_precond

  subroutine smooth_kernel(this)
    class(PostFlow), intent(inout) :: this
    integer :: iker

    do iker = 1, this%nker
      call smooth_sem_pde(this%ker_data(:,:,:,:,iker), fpar%sim%SIGMA_H, fpar%sim%SIGMA_V)
    enddo

  end subroutine smooth_kernel

  subroutine invert_hess( hess_matrix )

    ! inverts the Hessian matrix
    ! the approximate Hessian is only defined for diagonal elements: like
    ! H_nn = \frac{ \partial^2 \chi }{ \partial \rho_n \partial \rho_n }
    ! on all GLL points, which are indexed (i,j,k,ispec)
        
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB), intent(inout) :: hess_matrix
    real(kind=CUSTOM_REAL) :: maxh,maxh_all
  
    ! maximum value of Hessian
    maxh = maxval( abs(hess_matrix) )
  
    ! determines maximum from all slices on main
    call max_all_all_cr(maxh, maxh_all)
  
    ! normalizes Hessian
    if (maxh_all < 1.e-18) then
      ! Hessian is zero, re-initializes
      hess_matrix = 1.0_cr
      !stop 'Error Hessian too small'
    else
      ! since Hessian has absolute values, this scales between [0,1]
      hess_matrix = hess_matrix / maxh_all
    endif

    ! inverts Hessian values
    where( abs(hess_matrix(:,:,:,:)) > THRESHOLD_HESS )
      hess_matrix = 1.0_cr / hess_matrix
    elsewhere
      hess_matrix = 1.0_cr / THRESHOLD_HESS
    endwhere
  
    maxh = maxval( abs(hess_matrix) )
    call max_all_all_cr(maxh, maxh_all)
    hess_matrix = hess_matrix / maxh_all
    ! rescales Hessian
    !hess_matrix = hess_matrix * maxh_all
  
  end subroutine invert_hess
  
  subroutine write(this, is_smooth)
    class(PostFlow), intent(inout) :: this
    logical, optional, intent(in) :: is_smooth
    integer :: iker, ievt
    character(len=MAX_STRING_LEN) :: suffix 
    logical :: is_simu_type

    if (present(is_smooth)) then
      suffix = '_kernel_smooth'
    else
      suffix = '_kernel'
    endif

    if (is_joint) then
      is_simu_type = .true.
    else
      is_simu_type = .false.
    endif
    do iker = 1, this%nker
      call write_kernel(trim(this%ker_names(iker))//trim(suffix), this%ker_data(:,:,:,:,iker), is_simu_type)
    enddo

  end subroutine write

  subroutine sum_joint_kernel(this)
    class(PostFlow), intent(inout) :: this
    real(kind=cr), dimension(:,:,:,:),allocatable :: kernel1, total_kernel
    real(kind=cr) :: norm, norm_glob
    integer :: iker, itype
    real(kind=cr), dimension(NUM_INV_TYPE) :: norm_val
    character(len=MAX_STRING_LEN) :: output_dir, type_name, msg

    call log%write('This is taking sum of kernels for joint inversion...', .true.)

    do itype = 1, NUM_INV_TYPE
      if (.not. fpar%postproc%INV_TYPE(itype)) cycle
      call this%calc_kernel0_std_weight(itype, norm_val(itype))
    enddo

    do iker = 1, this%nker
      total_kernel = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_AB)
      do itype = 1, NUM_INV_TYPE
        if (.not. fpar%postproc%INV_TYPE(itype)) cycle
        type_name = this%simu_types(itype)
        output_dir = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)//'_'//trim(type_name)
        call read_kernel(output_dir, trim(this%ker_names(iker))//'_kernel_smooth', kernel1)
        kernel1 = fpar%postproc%JOINT_WEIGHT(itype)*kernel1 / norm_val(itype)
        norm = maxval(abs(kernel1))
        call max_all_cr(norm, norm_glob)
        write(msg, '(a,e13.8)') 'Max '//trim(this%ker_names(iker))//' kernel of '//trim(type_name)//&
                               ': ', norm_glob
        call log%write(msg, .true.)
        total_kernel = total_kernel + kernel1
      enddo
      ! output_dir = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)
      call write_kernel(trim(this%ker_names(iker))//'_kernel_smooth', total_kernel, .false.)
    enddo

  end subroutine sum_joint_kernel

  subroutine calc_kernel0_std_weight(this, itype, max_global)
    class(PostFlow), intent(inout) :: this
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_data
    real(kind=cr) :: max_local, max_global
    real(kind=cr), dimension(:), allocatable :: max_ker
    integer :: n, itype, iker, i
    character(len=MAX_STRING_LEN) :: type_name, output_dir, model0

    max_ker = zeros(this%nker)
    type_name = this%simu_types(itype)
    output_dir = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)//'_'//trim(type_name)
    do iker = 1, this%nker
      call read_kernel(output_dir, trim(this%ker_names(iker))//'_kernel_smooth', kernel_data)
      max_local = maxval( abs(kernel_data))
      call max_all_dp(max_local, max_ker(iker))
      call bcast_all_singledp(max_ker(iker))
    enddo
    max_global = maxval(max_ker(:))
    call synchronize_all()
  end subroutine calc_kernel0_std_weight

end module post_processing