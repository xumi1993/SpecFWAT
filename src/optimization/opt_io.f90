module opt_io
  use config
  use fwat_mpi
  use specfem_par
  use hdf5_interface
  use external_model
  use utils
  implicit none
  
  integer, private :: ier

contains

  subroutine read_model(iter, model_data)
    integer, intent(in) :: iter
    character(len=MAX_STRING_LEN) :: path, this_model, fprname
    integer :: imod
    real(kind=cr), dimension(:,:,:,:,:), allocatable, intent(out) :: model_data

    allocate(model_data(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nkernel))
    write(this_model, '(A1,I2.2)') 'M', iter
    path = trim(OPT_DIR)//'/MODEL_'//trim(this_model)
    call create_name_database(fprname, worldrank, path)
    do imod = 1, nkernel
      open(unit=IIN, file=trim(fprname)//trim(parameter_names(imod))//'.bin',&
           status='old', action='read', form='unformatted', iostat=ier)
      if (ier /= 0) then
        write(0, *) 'Error could not open database file: ',trim(fprname)//trim(parameter_names(imod))//'.bin'
        call exit_mpi(worldrank, 'Error opening database file')
      endif

      read(IIN) model_data(:,:,:,:,imod)
      close(IIN)
    end do

  end subroutine read_model

  subroutine read_model_grid_iso(iter, model_data)
    integer, intent(in) :: iter
    character(len=MAX_STRING_LEN) :: path, this_model, fprname
    integer :: imod
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: model_data
    real(kind=cr), dimension(:,:,:), allocatable :: gm
    type(hdf5_file) :: h5file

    if (worldrank == 0) then
      allocate(model_data(ext_grid%nx,ext_grid%ny,ext_grid%nz,3))
      write(this_model, '(A1,I2.2)') 'M', iter
      fprname = trim(OPT_DIR)//'/model_'//trim(this_model)//'.h5'
      call h5file%open(fprname, status='old', action='read')
      do imod = 1, 3
        call h5file%get('/'//trim(MODEL_ISO(imod)), gm)
        model_data(:,:,:,imod) = transpose_3(gm)
      end do
      call h5file%close(finalize=.true.)
    endif

  end subroutine read_model_grid_iso

  subroutine read_model_grid(iter, model_data)
    integer, intent(in) :: iter
    character(len=MAX_STRING_LEN) :: path, this_model, fprname
    integer :: imod
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: model_data
    real(kind=cr), dimension(:,:,:), allocatable :: gm
    type(hdf5_file) :: h5file

    if (worldrank == 0) then
      allocate(model_data(ext_grid%nx,ext_grid%ny,ext_grid%nz,nkernel))
      write(this_model, '(A1,I2.2)') 'M', iter
      fprname = trim(OPT_DIR)//'/model_'//trim(this_model)//'.h5'
      call h5file%open(fprname, status='old', action='read')
      do imod = 1, nkernel
        call h5file%get('/'//trim(parameter_names(imod)), gm)
        model_data(:,:,:,imod) = transpose_3(gm)
      end do
      call h5file%close(finalize=.true.)
    endif
  end subroutine read_model_grid

  subroutine read_gradient(iter, gradient_data)
    integer, intent(in) :: iter
    character(len=MAX_STRING_LEN) :: path, this_model, fprname
    integer :: imod
    real(kind=cr), dimension(:,:,:,:,:), allocatable, intent(out) :: gradient_data

    allocate(gradient_data(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nkernel))
    write(this_model, '(A1,I2.2)') 'M', iter
    path = trim(OPT_DIR)//'/SUM_KERNELS_'//trim(this_model)
    call create_name_database(fprname, worldrank, path)
    do imod = 1, nkernel
      open(unit=IIN, file=trim(fprname)//trim(kernel_names(imod))//'_kernel_smooth.bin',&
           status='old', action='read', form='unformatted', iostat=ier)
      if (ier /= 0) then
        write(0, *) 'Error could not open database file: ',trim(fprname)//trim(kernel_names(imod))//'.bin'
        call exit_mpi(worldrank, 'Error opening database file')
      endif

      ! allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      read(IIN) gradient_data(:,:,:,:,imod)
      close(IIN)
    end do
  end subroutine read_gradient

  subroutine read_gradient_grid(iter, gradient_data)
    integer, intent(in) :: iter
    character(len=MAX_STRING_LEN) :: path, this_model, fprname
    integer :: imod
    real(kind=cr), dimension(:,:,:,:), allocatable, intent(out) :: gradient_data
    real(kind=cr), dimension(:,:,:), allocatable :: gm
    type(hdf5_file) :: h5file

    if (worldrank == 0) then
      allocate(gradient_data(ext_grid%nx,ext_grid%ny,ext_grid%nz,nkernel))
      write(this_model, '(A1,I2.2)') 'M', iter
      fprname = trim(OPT_DIR)//'/gradient_'//trim(this_model)//'.h5'
      call h5file%open(fprname, status='old', action='read')
      do imod = 1, nkernel
        call h5file%get('/'//trim(kernel_names(imod))//'_kernel_smooth', gm)
        gradient_data(:,:,:,imod) = transpose_3(gm)
      end do
      call h5file%close(finalize=.true.)
    endif
  end subroutine read_gradient_grid

  subroutine write_model(OUTPUT_MODEL_PATH, model_data)
    character(len=*), intent(in) :: OUTPUT_MODEL_PATH
    real(kind=cr), dimension(:,:,:,:,:), intent(in) :: model_data
    character(len=MAX_STRING_LEN) :: path, this_model, fprname
    integer :: imod

    call create_name_database(fprname, worldrank, OUTPUT_MODEL_PATH)
    do imod = 1, nkernel
      open(unit=IIN, file=trim(fprname)//trim(parameter_names(imod))//'.bin',&
           status='replace', action='write', form='unformatted', iostat=ier)
      if (ier /= 0) then
        write(0, *) 'Error could not open database file: ',trim(fprname)//trim(parameter_names(imod))//'.bin'
        call exit_mpi(myrank,'Error opening database file')
      endif

      ! allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      write(IIN) model_data(:,:,:,:,imod)
      close(IIN)
    end do
  end subroutine write_model

  subroutine Parallel_ComputeInnerProduct(vect1, vect2, qp)

    ! use specfem_par, only: NSPEC_AB,  jacobianstore, jacobian_regular, irregular_element_number, wxgll, wygll, wzgll

    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(in)    :: vect1, vect2
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable   :: wks_1n, wks_2n
    real(kind=CUSTOM_REAL),                       intent(inout) :: qp
    integer                                                     :: Niv
    ! local
    real(kind=CUSTOM_REAL)                                      :: coeff, coeff_n1, coeff_n2
    ! try double precision
    double precision                                            :: jacobianl, weight, qp_dp, qp_tmp
    double precision                                            :: coeff_n1_dp, coeff_n2_dp
    integer                                                     :: ipar, i, j, k, ispec, ispec_irreg

    ! initializes
    qp = 0._CUSTOM_REAL

    Niv = size(vect1, 5)

    !! try normalization to avoid numerical overflow errors
    ! wks_1n(:,:,:,:,:) = 0.d0
    ! wks_2n(:,:,:,:,:) = 0.d0

    !call Parallel_ComputeL2normSquare(vect1 , Niv, coeff_n1)
    !call Parallel_ComputeL2normSquare(vect2 , Niv, coeff_n2)

    coeff = maxval(abs(vect1(:,:,:,:,:)))
    call max_all_all_cr(coeff, coeff_n1)

    if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1 = 1._CUSTOM_REAL

    coeff = maxval(abs(vect2(:,:,:,:,:)))
    call max_all_all_cr(coeff, coeff_n2)

    if (coeff_n2 == 0._CUSTOM_REAL) coeff_n2 = 1._CUSTOM_REAL

    coeff_n1_dp = coeff_n1
    coeff_n2_dp = coeff_n2

    ! wks_1n(:,:,:,:,:) = vect1(:,:,:,:,:)
    ! wks_2n(:,:,:,:,:) = vect2(:,:,:,:,:)

    wks_1n(:,:,:,:,:) = vect1(:,:,:,:,:) / coeff_n1_dp
    wks_2n(:,:,:,:,:) = vect2(:,:,:,:,:) / coeff_n2_dp

    ! L2 of normalized gradient
    qp_dp = 0.d0
    ! inner product with quadrature weights
    do ipar = 1, Niv
      do ispec = 1, NSPEC_AB
        ispec_irreg = irregular_element_number(ispec)
        if (ispec_irreg == 0) jacobianl = jacobian_regular
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              weight = wxgll(i) * wygll(j) * wzgll(k)
              if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)
              ! integrated inner product
              qp_dp = qp_dp + jacobianl * weight * wks_1n(i,j,k,ispec,ipar) * wks_2n(i,j,k,ispec,ipar)
              !qp = qp + jacobianl * weight * vect1(i,j,k,ispec,ipar) * vect2(i,j,k,ispec,ipar)
            enddo
          enddo
        enddo
      enddo
    enddo

    ! rescales inner product value
    qp_tmp = qp_dp * coeff_n1_dp * coeff_n2_dp

    ! due to overflow issues for larger meshes, uses double precision to sum over all slices
    call sum_all_all_dp(qp_tmp, qp_dp)

    ! converts to single precision
    if (abs(qp_dp) > HUGE(qp)) then
      write(*,*) "WARNING: Inner product numerical overflow, exceeds floating precision!!!"
      write(*,*) "         This will affect the Wolfe conditions, but will continue for now... "
      write(*,*) "         Please check your inputs (e.g.,source strengths) to have proper FWI setup."
      qp = sign(1.d0,qp_dp) * HUGE(qp) ! would overflow single precision
    else
      qp = real(qp_dp,kind=CUSTOM_REAL)
    endif

  end subroutine Parallel_ComputeInnerProduct

  subroutine model_interpolation(x, y, z, model_data, model)
    ! It is a placeholder for future implementation of model interpolation
    real(kind=cr), dimension(:), allocatable, intent(in) :: x, y, z
    real(kind=cr), dimension(:,:,:), allocatable, intent(in) :: model_data
    real(kind=cr), dimension(:,:,:), allocatable, intent(out) :: model
    integer :: i, j, k, nx, ny, nz
    
    nx = size(x)
    ny = size(y)
    nz = size(z)

    allocate(model(ext_grid%nx, ext_grid%ny, ext_grid%nz))
    do i = 1, ext_grid%nx
      do j = 1, ext_grid%ny
        do k = 1, ext_grid%nz
          if (ext_grid%x(i) < x(1) .or. ext_grid%x(i) > x(nx) .or. &
              ext_grid%y(j) < y(1) .or. ext_grid%y(j) > y(ny) .or. &
              ext_grid%z(k) < z(1) .or. ext_grid%z(k) > z(nz)) then
            model(i,j,k) = interp3_nearest_simple(x, y, z, model_data, ext_grid%x(i), ext_grid%y(j), ext_grid%z(k))
          else
            model(i,j,k) = interp3(x, y, z, model_data, ext_grid%x(i), ext_grid%y(j), ext_grid%z(k))
          endif
        enddo
      enddo
    enddo

  end subroutine model_interpolation
end module opt_io