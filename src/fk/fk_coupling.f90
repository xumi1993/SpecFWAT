module fk_coupling
  use specfem_par
  use specfem_par_coupling
  use config, only: worldrank, local_path_backup, compress_level
  use fwat_constants, only: cr, FKMODEL_PREFIX, SRC_REC_DIR, DEG2RAD
  use common_lib, only: mkdir
  use logger, only: log
#ifdef USE_CUDA
  use fk_gpu_backend, only: compute_fk_gpu
#endif

  implicit none
  integer, private :: ierr
#ifdef USE_CUDA
  logical, parameter :: fk_gpu_available = .true.
#else
  logical, parameter :: fk_gpu_available = .false.
#endif
contains

  subroutine read_fk_model(evtid)
    character(len=*), intent(in) :: evtid
    integer :: ier
    real(kind=cr) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box

    FKMODEL_FILE = trim(SRC_REC_DIR)//'/'//trim(FKMODEL_PREFIX)//'_'//trim(evtid)

    call FindBoundaryBox(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)
    call ReadFKModelInput(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)

    ! send FK parameters to others MPI slices
    call bcast_all_singlei(type_kpsv_fk)
    call bcast_all_singlei(nlayer)
    
    if (myrank > 0) then
      allocate(alpha_FK(nlayer), &
                beta_FK(nlayer), &
                rho_FK(nlayer), &
                mu_FK(nlayer), &
                h_FK(nlayer),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating arrays 2206')
      alpha_FK(:) = 0._cr; beta_FK(:) = 0._cr; rho_FK(:) = 0._cr
      mu_FK(:) = 0._cr; h_FK(:) = 0._cr
    endif

    call bcast_all_cr(alpha_FK, nlayer)
    call bcast_all_cr(beta_FK, nlayer)
    call bcast_all_cr(rho_FK, nlayer)
    call bcast_all_cr(mu_FK, nlayer)
    call bcast_all_cr(h_FK, nlayer)

    call bcast_all_singlecr(phi_FK)
    call bcast_all_singlecr(theta_FK)

    call bcast_all_singlecr(ff0)
    call bcast_all_singlecr(freq_sampling_fk)
    call bcast_all_singlecr(amplitude_fk)

    call bcast_all_singlecr(xx0)
    call bcast_all_singlecr(yy0)
    call bcast_all_singlecr(zz0)
    call bcast_all_singlecr(Z_REF_for_FK)

    call bcast_all_singlecr(tt0)
    call bcast_all_singlecr(tmax_fk)
    call synchronize_all()

  end subroutine read_fk_model

  subroutine free_fk_arrays()
    if (allocated(alpha_FK)) deallocate(alpha_FK)
    if (allocated(beta_FK)) deallocate(beta_FK)
    if (allocated(rho_FK)) deallocate(rho_FK)
    if (allocated(mu_FK)) deallocate(mu_FK)
    if (allocated(h_FK)) deallocate(h_FK)

  end subroutine free_fk_arrays


  subroutine compute_fk_wavefield(evtid)
    ! Shared by standalone FK and preproc; leave wavefield arrays available for SEM injection.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use input_params, only: fpar => fwat_par_global
    use hdf5, only: h5open_f
    character(len=*), intent(in) :: evtid
    real(kind=CUSTOM_REAL) :: p, tg, df
    integer :: status, layer, shift, saved_compression
    if (GPU_MODE .and. .not. fk_gpu_available) &
      call exit_MPI(myrank, 'GPU_MODE requires a USE_CUDA build for FK')
    call log%write('Computing FK wavefield for event '//trim(evtid), .true.)
    call read_fk_model(evtid)
    if (ff0 <= 0 .or. freq_sampling_fk <= 0 .or. tmax_fk <= 0) &
      call exit_MPI(myrank, 'FK frequencies and time window must be positive')
    if (any(alpha_FK <= 0) .or. any(beta_FK < 0) .or. any(rho_FK <= 0)) &
      call exit_MPI(myrank, 'Invalid FK layer properties')
    if (beta_FK(nlayer) <= 0) call exit_MPI(myrank, 'FK requires an elastic bottom half-space')
    zz0 = zz0 - Z_REF_for_FK
    phi_FK = phi_FK * PI/180.d0
    theta_FK = theta_FK * PI/180.d0
    if (type_kpsv_fk == 1) then
      p = sin(theta_FK)/alpha_FK(nlayer)
    else
      p = sin(theta_FK)/beta_FK(nlayer)
    endif
    if (abs(p) < 1.e-14_CUSTOM_REAL) p = sign(1.e-14_CUSTOM_REAL,p)
    ! The legacy FK takes real square roots; do not silently generate NaNs.
    if (any(1._CUSTOM_REAL/alpha_FK**2 - p**2 <= 0)) &
      call exit_MPI(myrank, 'Critical/evanescent P slowness is unsupported by the reference FK')
    do layer=1,nlayer
      if (beta_FK(layer) <= 0) cycle
      if (1._CUSTOM_REAL/beta_FK(layer)**2-p**2 <= 0) &
        call exit_MPI(myrank, 'Critical/evanescent S slowness is unsupported by the reference FK')
    enddo
    tg = 1.d0/ff0
    deltat = real(DT,CUSTOM_REAL)
    call count_num_boundary_points(num_abs_boundary_faces, abs_boundary_ispec, npt)
    call find_size_of_working_arrays(deltat, freq_sampling_fk, tmax_fk, NF_FOR_STORING, &
      NF_FOR_FFT, NPOW_FOR_INTERP, NP_RESAMP, df)
    if (NP_RESAMP < 1 .or. NP_RESAMP > 10000) call exit_MPI(myrank, 'Invalid FK resampling rate')
    shift=int(-tt0/deltat)
    if (tt0 > 0 .or. shift > min(NSTEP,NF_FOR_FFT)) &
      call exit_MPI(myrank, 'ORIGIN_TIME must be zero or negative, with shift <= min(NSTEP, FFT length)')
    if (NSTEP/NP_RESAMP > NF_FOR_STORING + NP_RESAMP) &
      call exit_MPI(myrank, 'FK time window is too short for NSTEP')
    ! Consecutive preproc events may leave either computed or cached wavefield arrays allocated.
    if (allocated(ipt_table)) deallocate(ipt_table)
    if (allocated(Veloc_FK)) deallocate(Veloc_FK)
    if (allocated(Tract_FK)) deallocate(Tract_FK)
    allocate(ipt_table(NGLLSQUARE,num_abs_boundary_faces), &
      Veloc_FK(NDIM,npt,-NP_RESAMP:NF_FOR_STORING+NP_RESAMP), &
      Tract_FK(NDIM,npt,-NP_RESAMP:NF_FOR_STORING+NP_RESAMP), stat=status)
    if (status /= 0) call exit_MPI(myrank, 'Cannot allocate FK arrays')
    ipt_table = 0
    Veloc_FK = 0
    Tract_FK = 0
    if (myrank == 0) write(*,'(a,i0,a,i0,a,i0)') &
      'FK: rank-0 boundary points=',npt, ', FFT length=',NF_FOR_FFT, ', resampling=',NP_RESAMP
    if (npt > 0) then
      if (GPU_MODE) then
#ifdef USE_CUDA
        call compute_fk_gpu(p, tg, df)
#else
        call exit_MPI(myrank, 'GPU_MODE requires a USE_CUDA build for FK')
#endif
      else
        ! Retain the reference kernel and all its numerical conventions.
        call FK3D_fwat(type_kpsv_fk, nlayer, NSTEP, npt, p, phi_FK, xx0, yy0, zz0, tg, &
          tt0, alpha_FK, beta_FK, rho_FK, h_FK, NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, df)
      endif
    endif
    if (.not. all(ieee_is_finite(Veloc_FK)) .or. .not. all(ieee_is_finite(Tract_FK))) &
      call exit_MPI(myrank, 'Non-finite FK output; nothing will be saved')
    call free_fk_arrays()
    call synchronize_all()
    call log%write('Finished FK wavefield computation', .true.)

    ! Keep the original cache writer and format; SAVE_FK controls saving, independently of GPU_MODE.
    if (fpar%sim%SAVE_FK) then
      call log%write('Saving FK wavefield', .true.)
      call h5open_f(status)
      if (status /= 0) call exit_MPI(myrank, 'Cannot initialize HDF5')
      ! The writer's fixed chunks require at least 16 points and 128 time samples.
      saved_compression = compress_level
      if (size(Veloc_FK,2) < 16 .or. size(Veloc_FK,3) < 128) compress_level = 0
      call write_fk_coupling_file(evtid)
      compress_level = saved_compression
      call log%write('Saved '//trim(local_path_backup)//'/FK_wavefield_'//trim(evtid), .true.)
    else
      call log%write('TELE.SAVE_FK is false; wavefield was computed without saving', .true.)
    endif
  end subroutine compute_fk_wavefield

  subroutine initialize_ipt_table()
    use specfem_par_elastic, only: ispec_is_elastic
    use specfem_par_acoustic, only: ispec_is_acoustic 
    integer :: ier, ipt, ispec, igll, iface

    if (allocated(ipt_table)) deallocate(ipt_table)
    if (npt > 0) then
      allocate(ipt_table(NGLLSQUARE,num_abs_boundary_faces), stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2202')
      ipt_table(:,:) = 0
    else
      ! dummy
      allocate(ipt_table(1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2204')
    endif
    ipt_table(:,:) = 0

    ipt = 0
    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
      if (ispec_is_elastic(ispec) .or. ispec_is_acoustic(ispec)) then
        ! reference GLL points on boundary face
        do igll = 1,NGLLSQUARE
          ipt = ipt + 1
          ipt_table(igll,iface) = ipt
        end do
      end if
    end do
  end subroutine initialize_ipt_table

  subroutine write_fk_coupling_file(evtid)
    use hdf5
    use hdf5_interface

    character(len=MAX_STRING_LEN), intent(in) :: evtid
    integer(HSIZE_T), dimension(3) :: dims, chunk_dims
    character(len=MAX_STRING_LEN) :: fname, out_dir
    integer(HID_T) :: file_id, dset_idv, dset_idt, dataspace_id, plist_id
    integer, parameter :: rank_fk=3
    type(hdf5_file) :: h5file

    dims = shape(Veloc_FK)
    chunk_dims = (/NDIM, int(dims(2)/16), 128/)

    out_dir = trim(local_path_backup)//'/FK_wavefield_'//trim(evtid)//'/'
    ! if (worldrank == 0) call system('mkdir -p '//trim(out_dir))
    call mkdir(out_dir)
    call synchronize_all()
    

    write(fname,'(a,i6.6,a)') trim(out_dir)//'proc', worldrank, '_fk_wavefield.h5'
    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, ierr)
    if (ierr /= 0) then
      call exit_MPI(worldrank, 'error creating file 2205')
    end if

    if (compress_level <= 0) then
      call h5file%open(fname)
      call h5file%add("/Veloc_FK", Veloc_FK)
      call h5file%add("/Tract_FK", Tract_FK)
      call h5file%close(finalize=.true.)
    else
      ! create HDF5 file
      call h5open_f(ierr)

      ! create dataspace for the dataset
      call h5screate_simple_f(rank_fk, dims, dataspace_id, ierr)

      ! create property list for chunking
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, ierr)
      call h5pset_chunk_f(plist_id, rank_fk, chunk_dims, ierr)
      call h5pset_deflate_f(plist_id, compress_level, ierr)  ! compression level

      ! create dataset for velocity
      call h5dcreate_f(file_id, "Veloc_FK", H5T_NATIVE_REAL, dataspace_id, &
                        dset_idv, ierr, plist_id)

      ! create dataset for traction
      call h5dcreate_f(file_id, "Tract_FK", H5T_NATIVE_REAL, dataspace_id, &
                        dset_idt, ierr, plist_id)

      ! write data to dataset
      call h5dwrite_f(dset_idv, H5T_NATIVE_REAL, Veloc_FK, dims, ierr)
      call h5dwrite_f(dset_idt, H5T_NATIVE_REAL, Tract_FK, dims, ierr)

      ! close dataset and file
      call h5dclose_f(dset_idv, ierr)
      call h5dclose_f(dset_idt, ierr)
      call h5sclose_f(dataspace_id, ierr)
      call h5pclose_f(plist_id, ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
    end if
    call synchronize_all()

  end subroutine write_fk_coupling_file

  subroutine read_fk_coupling_file(evtid)
    use hdf5_interface
    integer :: ier
    character(len=*) :: evtid
    character(len=MAX_STRING_LEN) :: out_dir, fkprname
    real(kind=cr) :: DF_FK
    type(hdf5_file) :: h5file
    real(kind=cr), dimension(:,:,:), allocatable :: tmp_array

    call count_num_boundary_points(num_abs_boundary_faces,abs_boundary_ispec,npt)
    call read_fk_model(evtid)
    deltat = real(DT, cr)
    call find_size_of_working_arrays(deltat, freq_sampling_fk, tmax_fk, NF_FOR_STORING, &
                                     NF_FOR_FFT, NPOW_FOR_INTERP, NP_RESAMP, DF_FK)
    
    call initialize_ipt_table()
    if (allocated(Veloc_FK)) deallocate(Veloc_FK)
    if (allocated(Tract_FK)) deallocate(Tract_FK)

    allocate(Veloc_FK(NDIM, npt, -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
    if (ier /= 0) call exit_MPI(worldrank, 'error allocating array 2210')

    allocate(Tract_FK(NDIM, npt, -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
    if (ier /= 0) call exit_MPI(worldrank, 'error allocating array 2210')

    out_dir = trim(local_path_backup)//'/FK_wavefield_'//trim(evtid)//'/'
    write(fkprname,'(a,i6.6,a)') trim(out_dir)//'proc', worldrank, '_fk_wavefield.h5'
    ! read FK solution from file
    ! open(FID, file=fkprname, form='unformatted', status='old', action='read', iostat=ier)
    ! if (ier /= 0) call exit_MPI(worldrank, 'error opening file 2205')
    ! read(FID) Veloc_FK, Tract_FK
    ! close(FID)
    call h5file%open(fkprname)
    call h5file%get("/Veloc_FK", tmp_array)
    Veloc_FK = tmp_array
    call h5file%get("/Tract_FK", tmp_array)
    Tract_FK = tmp_array
    call h5file%close(finalize=.true.)

    ! deallocate(alpha_FK, beta_FK, rho_FK, mu_FK, h_FK)
    call free_fk_arrays()

    call synchronize_all()
  end subroutine read_fk_coupling_file

  subroutine fktime(xx, yy, zz, tdelay)
    real(kind=cr), intent(in) :: xx, yy, zz
    real(kind=cr), intent(out) :: tdelay
    real(kind=cr), dimension(:), allocatable :: h, v_fk_input
    real(kind=cr) :: p, eta, z1, z0, theta_rad, phi_rad
    integer :: ilayer, j
    
    theta_rad = theta_FK * DEG2RAD
    phi_rad = phi_FK * DEG2RAD
    if (type_kpsv_fk == 1) then
      v_fk_input = alpha_FK
      p = sin(theta_rad) / v_fk_input(nlayer)
    else if (type_kpsv_fk == 2) then
      v_fk_input = beta_FK
      p = sin(theta_rad) / v_fk_input(nlayer)
    endif
    tdelay = p * (xx - xx0) * cos(phi_rad) + p * (yy - yy0) * sin(phi_rad)

    z0 = zz0 - Z_REF_for_FK
    z1 = zz - Z_REF_for_FK

    ! h = zeros(nlayer)
    allocate(h(nlayer))
    h = 0._cr
    ilayer = nlayer
    do j = nlayer - 1, 1, -1
      if (z1 <= sum(h_FK(j:nlayer))) then
        ilayer = j
        exit
      end if
    end do
    h(ilayer+1:nlayer) = h_FK(ilayer+1:nlayer)
    h(ilayer) = z1 - sum(h_FK(ilayer+1:nlayer))
    h(nlayer) = 0 - z0
    if (h(ilayer) < 0) then
      print *, 'Error setting layer thickness'
      stop
    end if
    do j = nlayer, ilayer, -1
      eta = sqrt(1 / v_fk_input(j)**2 - p**2)
      tdelay = tdelay + eta * h(j)
    end do

  end subroutine fktime
  
  function check_fk_files(evtid) result(res)
    use fwat_mpi, only: land_all_all_l

    character(len=*) :: evtid
    character(len=MAX_STRING_LEN) :: out_dir, fkprname
    logical :: findfile, res

    out_dir = './'//trim(local_path_backup)//'/FK_wavefield_'//trim(evtid)//'/'
    write(fkprname,'(a,i6.6,a)') trim(out_dir)//'proc', worldrank, '_fk_wavefield.h5'
    inquire(file=trim(fkprname), exist=findfile)
    call land_all_all_l(findfile, res)

  end function check_fk_files

end module