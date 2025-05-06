module zprecond
  use input_params, fpar => fwat_par_global
  use config
  use fwat_constants
  use utils, only: zeros, interp1
  use external_model
  implicit none

contains

  subroutine get_1d_precond(zl)
    integer :: iz
    real(kind=cr), dimension(:), allocatable, intent(out) :: zl

    zl = zeros(MEXT_V%nz)

    do iz = 1, MEXT_V%nz
      if (MEXT_V%z(iz) > 0.0) then
        zl(iz) = 1.e-8
      else
        zl(iz) = MEXT_V%z(iz)+1.e-8
      end if
      
      select case (fpar%sim%PRECOND_TYPE)
        case (2)
          zl(iz) = 1/sqrt(zl(iz)**2)
        case (3)
          zl(iz) = 1/sqrt(sqrt(zl(iz)**2))
      end select
    end do

  end subroutine get_1d_precond

  subroutine zprecond_gll(hess)
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: hess
    real(kind=cr), dimension(:), allocatable :: zl
    real(kind=dp), dimension(:), allocatable :: zl_dp
    real(kind=dp) :: zz
    real(kind=cr) :: maxh,maxh_all
    integer :: iz, ix, iy, ispec, iglob

    call get_1d_precond(zl)
    zl_dp = dble(zl)

    hess = zeros(NGLLX, NGLLY, NGLLZ, NSPEC_FWAT)
    
    do ispec = 1, NSPEC_FWAT
      do ix = 1, NGLLX; do iy = 1, NGLLY; do iz = 1, NGLLZ
        iglob = ibool_fwat(ix, iy, iz, ispec)
        zz = interp1(dble(MEXT_V%z), zl_dp, dble(zstore_fwat(iglob)))
        hess(ix, iy, iz, ispec) = sngl(zz)
      enddo; enddo; enddo
    enddo
    hess = hess * fpar%acqui%nevents
    
    ! scales between [0,1]
    maxh = maxval(abs(hess))
    call max_all_all_cr(maxh, maxh_all)
    hess = hess / maxh_all
    
    ! inverts Hessian values
    hess = 1.0_cr / hess

    ! normalizes Hessian
    maxh = maxval(abs(hess))
    call max_all_all_cr(maxh, maxh_all)
    hess = hess / maxh_all

    call synchronize_all()
  end subroutine zprecond_gll

  subroutine zprecond_grid(hess)
    real(kind=cr), dimension(:,:,:), allocatable, intent(out) :: hess
    real(kind=cr), dimension(:), allocatable :: zl
    real(kind=cr) :: maxh
    integer :: iz, ix, iy, ispec, iglob

    if (worldrank == 0) then
      call get_1d_precond(zl)

      hess = zeros(MEXT_V%nx, MEXT_V%ny, MEXT_V%nz)
      do ix = 1, MEXT_V%nx
        do iy = 1, MEXT_V%ny
          hess(ix, iy, :) = zl
        enddo
      enddo
      hess = hess * fpar%acqui%nevents

      maxh = maxval(abs(hess))
      hess = hess / maxh

      hess = 1.0_cr / hess

      maxh = maxval(abs(hess))
      hess = hess / maxh
    endif
    call synchronize_all()

  end subroutine zprecond_grid

end module zprecond