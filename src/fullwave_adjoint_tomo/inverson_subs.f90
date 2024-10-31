
subroutine init_inversion()

use fullwave_adjoint_tomo_par
use fwat_input
use generate_databases_par, only: IMODEL
use specfem_par
use fwat_utils, only: get_mesh_file_path

implicit none

integer :: itype
logical :: exist

if (myrank == 0) then
    ! check IMODEL
    if (IMODEL /= 7 .and. IMODEL /= 6) then
    print *, 'ERROR: MODEL must be gll or external'
    stop
    endif

    ! check for joint
    if (count(tomo_par%INV_TYPE) > 1) then
    is_joint = .true.
    do itype = 1, NUM_INV_TYPE
        if (tomo_par%INV_TYPE(itype)) then
        inquire(file=get_mesh_file_path(itype), exist=exist)
        if (.not. exist) then
            print *, 'ERROR: No found mesh file of ', trim(get_mesh_file_path(itype))
            stop
        endif
        endif
    end do
    endif ! is_joint
endif ! myrank == 0
call bcast_all_singlel(is_joint)

if (is_joint) then 
    call read_mesh_parameter_file_fwat(get_mesh_file_path(2))
    call meshfem3d_fwat()
    call generate_database_fwat(USE_H5)
    block
    use meshfem3D_par, only: NEX_XI, NEX_ETA, NER
    use generate_databases_par, only: xstore, ystore, zstore
    real(kind=CUSTOM_REAL) :: x_min_loc, x_max_loc, y_min_loc, y_max_loc, z_min_loc, z_max_loc
    real(kind=CUSTOM_REAL) :: x_min, x_max, y_min, y_max, z_min, z_max, deltax, deltay, deltaz

    x_min_loc = minval(xstore)
    x_max_loc = maxval(xstore)
    y_min_loc = minval(ystore)
    y_max_loc = maxval(ystore)
    z_min_loc = minval(zstore)
    z_max_loc = maxval(zstore)
    call min_all_all_cr(x_min_loc, x_min)
    call max_all_all_cr(x_max_loc, x_max)
    call min_all_all_cr(y_min_loc, y_min)
    call max_all_all_cr(y_max_loc, y_max)
    call min_all_all_cr(z_min_loc, z_min)
    call max_all_all_cr(z_max_loc, z_max)

    deltax = (x_max - x_min) / NEX_XI / NGLLX
    deltay = (y_max - y_min) / NEX_ETA / NGLLY
    deltaz = (z_max - z_min) / NER / NGLLZ

    call rg%init(x_min, x_max, y_min, y_max, z_min, z_max, deltax, deltay, deltaz)
    
    end block
endif
call synchronize_all()
end subroutine init_inversion
