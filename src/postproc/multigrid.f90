module multigrid
  use config
  use fwat_mpi
  use fwat_constants
  use utils
  use input_params, only: fpar => fwat_par_global
  use external_model
  use hdf5_interface

  implicit none

  type :: InvGrid
    real(kind=cr), dimension(:,:), allocatable :: xinv, yinv, zinv
    integer :: ninvx, ninvy, ninvz
  contains
    procedure, public :: sem2inv, init=>init_inversion_grid, inv2grid
  end type InvGrid

contains

  subroutine init_inversion_grid(this)
    class(InvGrid), intent(inout) :: this
    integer :: nset
    real(kind=cr) :: dinvx, dinvy, dinvz, xadd, yadd, zadd, &
                     x_beg, y_beg, z_beg, x_end, y_end, z_end

    real(kind=cr), dimension(:), allocatable :: x_inv_1d, y_inv_1d, z_inv_1d
    integer :: i, j
    type(hdf5_file) :: h5file

    nset = fpar%postproc%n_inversion_grid
    this%ninvx = fpar%postproc%ninv(1) + 2
    this%ninvy = fpar%postproc%ninv(2) + 2
    this%ninvz = fpar%postproc%ninv(3) + 2

    this%xinv = zeros(this%ninvx, nset)
    this%yinv = zeros(this%ninvy, nset)
    this%zinv = zeros(this%ninvz, nset)

    dinvx = (ext_grid%x(ext_grid%nx) - ext_grid%x(1))/ real(fpar%postproc%ninv(1), kind=cr)
    xadd = dinvx / real(nset, kind=cr)
    x_beg = ext_grid%x(1) - xadd
    x_end = ext_grid%x(ext_grid%nx) + xadd

    dinvy = (ext_grid%y(ext_grid%ny) - ext_grid%y(1))/ real(fpar%postproc%ninv(2), kind=cr)
    yadd = dinvy / real(nset, kind=cr)
    y_beg = ext_grid%y(1) - yadd
    y_end = ext_grid%y(ext_grid%ny) + yadd

    dinvz = (ext_grid%z(ext_grid%nz) - ext_grid%z(1))/ real(fpar%postproc%ninv(3), kind=cr)
    zadd = dinvz / real(nset, kind=cr)
    z_beg = ext_grid%z(1) - zadd
    z_end = ext_grid%z(ext_grid%nz) + zadd

    x_inv_1d = zeros(this%ninvx)
    y_inv_1d = zeros(this%ninvy)
    z_inv_1d = zeros(this%ninvz)

    do i = 1, fpar%postproc%ninv(1)
      x_inv_1d(i) = ext_grid%x(1) + (i-1) * dinvx
    end do

    do i = 1, fpar%postproc%ninv(2)
      y_inv_1d(i) = ext_grid%y(1) + (i-1) * dinvy
    end do

    do i = 1, fpar%postproc%ninv(3)
      z_inv_1d(i) = ext_grid%z(1) + (i-1) * dinvz
    end do

    do i = 1, nset
      this%xinv(1, i) = x_beg
      this%xinv(this%ninvx, i) = x_end
      this%yinv(1, i) = y_beg
      this%yinv(this%ninvy, i) = y_end
      this%zinv(1, i) = z_beg
      this%zinv(this%ninvz, i) = z_end

      do j = 2, this%ninvx - 1
        this%xinv(j, i) = x_inv_1d(j-1) + xadd * (i-1)
      end do
      do j = 2, this%ninvy - 1
        this%yinv(j, i) = y_inv_1d(j-1) + yadd * (i-1)
      end do
      do j = 2, this%ninvz - 1
        this%zinv(j, i) = z_inv_1d(j-1) + zadd * (i-1)
      end do
    end do 
    call synchronize_all()

    if (worldrank == 0 .and. is_output_inv_grid) then
      call h5file%open(trim(OUTPUT_FILES)//'/inversion_grid.h5', status='new', action='write')
      call h5file%add('/xinv', this%xinv)
      call h5file%add('/yinv', this%yinv)
      call h5file%add('/zinv', this%zinv)
      call h5file%close(finalize=.true.)
    endif
    call synchronize_all()
  end subroutine init_inversion_grid

  subroutine sem2inv(this, data, data_inv)
    class(InvGrid), intent(inout) :: this
    real(kind=cr), dimension(:,:,:,:), intent(in) :: data
    real(kind=cr), dimension(:,:), allocatable, intent(out) :: data_inv
    real(kind=cr), dimension(:,:), allocatable :: tmp
    real(kind=cr) :: wx, wy, wz, wt
    integer idx, idy, idz, m, ngrid, igrid, i, j, k, ispec, iglob, n

    ngrid = this%ninvx * this%ninvy * this%ninvz * fpar%postproc%n_inversion_grid
    data_inv = zeros(this%ninvx * this%ninvy * this%ninvz, fpar%postproc%n_inversion_grid)
    tmp = zeros(this%ninvx * this%ninvy * this%ninvz, fpar%postproc%n_inversion_grid)
    
    do igrid = 1, fpar%postproc%n_inversion_grid
      do ispec = 1, NSPEC_FWAT
        do i = 1, NGLLX
          do j = 1, NGLLY
            do k = 1, NGLLZ
              iglob = ibool_fwat(i, j, k, ispec)
              call locate_bissection(dble(this%xinv(:, igrid)), this%ninvx, dble(xstore_fwat(iglob)), idx)
              if (idx == -1) call exit_mpi(worldrank, 'ERROR MULTIGRID: x is out of boundary')
              wx = (xstore_fwat(iglob) - this%xinv(idx, igrid)) / (this%xinv(idx+1, igrid) - this%xinv(idx, igrid))
              call locate_bissection(dble(this%yinv(:, igrid)), this%ninvy, dble(ystore_fwat(iglob)), idy)
              if (idy == -1) call exit_mpi(worldrank, 'ERROR MULTIGRID: y is out of boundary')
              wy = (ystore_fwat(iglob) - this%yinv(idy, igrid)) / (this%yinv(idy+1, igrid) - this%yinv(idy, igrid))
              call locate_bissection(dble(this%zinv(:, igrid)), this%ninvz, dble(zstore_fwat(iglob)), idz)
              if (idz == -1) call exit_mpi(worldrank, 'ERROR MULTIGRID: z is out of boundary')
              wz = (zstore_fwat(iglob) - this%zinv(idz, igrid)) / (this%zinv(idz+1, igrid) - this%zinv(idz, igrid))
              do n = 1, 8
                select case (n)
                  case (1)
                    m = this%ninvx * this%ninvy * (idz-1) + this%ninvx * (idy-1) + idx
                    wt = (1.0_cr - wx) * (1.0_cr - wy) * (1.0_cr - wz)
                  case (2)
                    m = this%ninvx * this%ninvy * (idz-1) + this%ninvx * idy + idx
                    wt = (1.0_cr - wx) * wy * (1.0_cr - wz)
                  case (3)
                    m = this%ninvx * this%ninvy * (idz-1) + this%ninvx * idy + idx + 1
                    wt = wx * wy * (1.0_cr - wz)
                  case (4)
                    m = this%ninvx * this%ninvy * (idz-1) + this%ninvx * (idy-1) + idx + 1
                    wt = wx * (1.0_cr - wy) * (1.0_cr - wz)
                  case (5)
                    m = this%ninvx * this%ninvy * idz + this%ninvx * (idy-1) + idx
                    wt = (1.0_cr - wx) * (1.0_cr - wy) * wz
                  case (6)
                    m = this%ninvx * this%ninvy * idz + this%ninvx * idy + idx
                    wt = (1.0_cr - wx) * wy * wz
                  case (7)
                    m = this%ninvx * this%ninvy * idz + this%ninvx * idy + idx + 1
                    wt = wx * wy * wz
                  case (8)
                    m = this%ninvx * this%ninvy * idz + this%ninvx * (idy-1) + idx + 1
                    wt = wx * (1.0_cr - wy) * wz
                end select
                tmp(m, igrid) = tmp(m, igrid) + wt * data(i, j, k, ispec)
              end do
            end do
          end do
        end do
      end do
    end do
    call synchronize_all()
    call sum_all_1Darray_cr(tmp, data_inv, ngrid)
    call bcast_all_cr(data_inv, ngrid)
  end subroutine sem2inv

  subroutine inv2grid(this, data_inv, data_grid)
    use external_model
    class(InvGrid), intent(inout) :: this
    real(kind=cr), dimension(:,:), intent(in) :: data_inv
    real(kind=cr), dimension(:,:,:), allocatable, intent(out) :: data_grid
    integer :: idx, idy, idz, m, ngrid, igrid, i, j, k, n
    real(kind=cr) :: wx, wy, wz, wt, val

    ngrid = this%ninvx * this%ninvy * this%ninvz * fpar%postproc%n_inversion_grid
    data_grid = zeros(ext_grid%nx, ext_grid%ny, ext_grid%nz)

    do igrid = 1, fpar%postproc%n_inversion_grid
      do i = 1, ext_grid%nx
        do j = 1, ext_grid%ny
          do k = 1, ext_grid%nz
            call locate_bissection(dble(this%xinv(:, igrid)), this%ninvx, dble(ext_grid%x(i)), idx)
            if (idx == -1) call exit_mpi(worldrank, 'ERROR MULTIGRID: x is out of boundary')
            wx = (ext_grid%x(i) - this%xinv(idx, igrid)) / (this%xinv(idx+1, igrid) - this%xinv(idx, igrid))
            call locate_bissection(dble(this%yinv(:, igrid)), this%ninvy, dble(ext_grid%y(j)), idy)
            if (idy == -1) call exit_mpi(worldrank, 'ERROR MULTIGRID: y is out of boundary')
            wy = (ext_grid%y(j) - this%yinv(idy, igrid)) / (this%yinv(idy+1, igrid) - this%yinv(idy, igrid))
            call locate_bissection(dble(this%zinv(:, igrid)), this%ninvz, dble(ext_grid%z(k)), idz)
            if (idz == -1) call exit_mpi(worldrank, 'ERROR MULTIGRID: z is out of boundary')
            wz = (ext_grid%z(k) - this%zinv(idz, igrid)) / (this%zinv(idz+1, igrid) - this%zinv(idz, igrid))
            val = 0.0_cr
            do n = 1, 8
              select case (n)
                case (1)
                  m = this%ninvx * this%ninvy * (idz-1) + this%ninvx * (idy-1) + idx
                  wt = (1.0_cr - wx) * (1.0_cr - wy) * (1.0_cr - wz)
                case (2)
                  m = this%ninvx * this%ninvy * (idz-1) + this%ninvx * idy + idx
                  wt = (1.0_cr - wx) * wy * (1.0_cr - wz)
                case (3)
                  m = this%ninvx * this%ninvy * (idz-1) + this%ninvx * idy + idx + 1
                  wt = wx * wy * (1.0_cr - wz)
                case (4)
                  m = this%ninvx * this%ninvy * (idz-1) + this%ninvx * (idy-1) + idx + 1
                  wt = wx * (1.0_cr - wy) * (1.0_cr - wz)
                case (5)
                  m = this%ninvx * this%ninvy * idz + this%ninvx * (idy-1) + idx
                  wt = (1.0_cr - wx) * (1.0_cr - wy) * wz
                case (6)
                  m = this%ninvx * this%ninvy * idz + this%ninvx * idy + idx
                  wt = (1.0_cr - wx) * wy * wz
                case (7)
                  m = this%ninvx * this%ninvy * idz + this%ninvx * idy + idx + 1
                  wt = wx * wy * wz
                case (8)
                  m = this%ninvx * this%ninvy * idz + this%ninvx * (idy-1) + idx + 1
                  wt = wx * (1.0_cr - wy) * wz
              end select
              val = val + wt * data_inv(m, igrid)
            end do
            data_grid(i, j, k) = data_grid(i, j, k) + val
          end do
        end do
      end do
    end do
    data_grid = data_grid / fpar%postproc%n_inversion_grid
    call synchronize_all()

  end subroutine inv2grid

end module multigrid