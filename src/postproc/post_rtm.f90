module post_rtm
  use config
  use fwat_mpi
  use fwat_constants
  use logger, only: log
  use input_params, fpar => fwat_par_global
  use utils, only: zeros
  use kernel_io
  use taper3d
  use common_lib, only: get_dat_type, get_kernel_name_rtm, mkdir
  use smooth_mod
  use zprecond
  use model_grid_data, only: create_grid, write_grid, gll2grid

  implicit none
  character(len=MAX_STRING_LEN), private :: msg

  type :: PostRTM
    real(kind=cr), dimension(:,:,:,:), allocatable :: ker_data
    real(kind=cr), dimension(:,:,:), allocatable :: ker_data_grid
    character(len=MAX_STRING_LEN) :: kernel_path

    contains
    procedure :: init=>init_post_flow
    procedure :: sum_kernel, apply_precond, finalize
    procedure :: project_kernel_to_grid, write_rtm_kernel

  end type PostRTM

contains

  subroutine init_post_flow(this)
    class(PostRTM), intent(inout) :: this
    call get_kernel_name_rtm()

    if (.not. use_gll) call create_grid()

    call log%init('output_post_processing_'//trim(model_name)//'.log')
    call log%write('*******************************************', .false.)

    ! read src_rec for this data type
    call fpar%acqui%read()

    ! setup mesh
    call read_mesh_databases_for_init()

    this%kernel_path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(model_name)
    call mkdir(this%kernel_path)
    call synchronize_all()

    this%ker_data = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_FWAT)

  end subroutine init_post_flow

  subroutine sum_kernel(this)
    class(PostRTM), intent(inout) :: this
    real(kind=cr), dimension(:,:,:,:), allocatable :: ker
    integer :: ievt

    call log%write('This is writing sum of kernels...', .true.)
    do ievt = 1, fpar%acqui%nevents
      call read_event_kernel(ievt, trim(kernel_names(1))//'_kernel', ker)
      this%ker_data(:,:,:,:) = this%ker_data(:,:,:,:) + ker
    enddo
    call synchronize_all()

  end subroutine sum_kernel

  subroutine apply_precond(this)
    class(PostRTM), intent(inout) :: this
    real(kind=cr), dimension(:,:,:,:), allocatable :: total_hess

    call log%write('This is Z-preconditioning kernels...', .true.)
    call zprecond_gll(total_hess) 
    call synchronize_all()

    this%ker_data = this%ker_data * total_hess
  end subroutine apply_precond

  subroutine project_kernel_to_grid(this)
    class(PostRTM), intent(inout) :: this
    if (use_gll) return
    call log%write('This is projecting kernels to grid...', .true.)
    call gll2grid(this%ker_data, this%ker_data_grid)
  end subroutine project_kernel_to_grid

  subroutine write_rtm_kernel(this)
    class(PostRTM), intent(inout) :: this
    if (use_gll) then
      call write_kernel(this%kernel_path, kernel_names(1), this%ker_data)
    else
      call write_grid(trim(this%kernel_path)//'/'//trim(kernel_names(1))//'.h5', &
                      kernel_names(1), this%ker_data_grid)
    endif
  end subroutine write_rtm_kernel

  subroutine finalize(this)
    class(PostRTM), intent(inout) :: this
    if (allocated(this%ker_data)) deallocate(this%ker_data)
  end subroutine finalize

end module post_rtm