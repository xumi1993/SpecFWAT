module post_processing
  use config
  use fwat_mpi
  use fwat_constants
  use logger, only: log
  use input_params, fpar => fwat_par_global
  use utils, only: zeros
  use kernel_io

  implicit none

  type PostProc
    character(len=MAX_STRING_LEN), dimension(:), allocatable :: ker_names
    real(kind=cr), dimension(:,:,:,:,:), allocatable :: ker_data
    integer :: nker
    contains
    procedure :: sum_kernel, init, sum_precond, write
  end type PostProc
contains

  subroutine init(this)
    class(PostProc), intent(inout) :: this

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

    if (worldrank == 0) call system('mkdir -p '//trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name))

    this%ker_data = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_AB, this%nker)

    call log%init('output_fwd_post_processing.log')

  end subroutine init

  subroutine sum_kernel(this)
    class(PostProc), intent(inout) :: this
    real(kind=cr), dimension(:,:,:,:), allocatable :: ker
    integer :: iker, ievt

    call log%write('This is sum kernel...', .true.)
    do iker = 1, this%nker
      do ievt = 1, fpar%acqui%nevents
        call read_kernel(ievt, trim(this%ker_names(iker))//'_kernel', ker)
        this%ker_data(:,:,:,:,iker) = this%ker_data(:,:,:,:,iker) + ker
      enddo
    enddo

  end subroutine sum_kernel

  subroutine sum_precond(this)
    class(PostProc), intent(inout) :: this
    character(len=MAX_STRING_LEN), parameter :: hess_name = 'hess'
    real(kind=cr), dimension(:,:,:,:), allocatable :: total_hess, ker
    real(kind=cr) :: zdep
    integer :: iker, ievt, igllx, iglly, igllz, ispec

    if (simu_type /= SIMU_TYPE_TELE) then
      total_hess = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_AB)
      do ievt = 1, fpar%acqui%nevents
        call read_precond(ievt, trim(hess_name)//'_kernel', ker)
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
              if (tele_par%PRECOND_TYPE == Z_PRECOND) then
                zdep = abs(zstore(igllx,iglly,igllz,ispec))/abs(z_min_glob)
              elseif (tele_par%PRECOND_TYPE == Z_SQRT_PRECOND) then
                zdep = sqrt(abs(zstore(igllx,iglly,igllz,ispec))/abs(z_min_glob))
              endif
              this%ker_data(igllx,iglly,igllz,ispec,:) = this%ker_data(igllx,iglly,igllz,ispec,:) / zdep
            enddo
          enddo
        enddo
      enddo
    endif

  end subroutine sum_precond

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
    class(PostProc), intent(inout) :: this
    logical, optional, intent(in) :: is_smooth
    integer :: iker, ievt
    character(len=MAX_STRING_LEN) :: suffix 

    if (present(is_smooth)) then
      suffix = '_kernel_smooth'
    else
      suffix = '_kernel'
    endif

    do iker = 1, this%nker
      call write_kernel(trim(this%ker_names(iker))//trim(suffix), this%ker_data(:,:,:,:,iker))
    enddo

  end subroutine write

end module post_processing