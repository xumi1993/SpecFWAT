module fk_coupling
  use specfem_par
  use specfem_par_coupling
  use config, only: worldrank, local_path_backup, compress_level
  use fwat_constants, only: cr, FKMODEL_PREFIX, SRC_REC_DIR, DEG2RAD
  use common_lib, only: mkdir

  implicit none
  integer, private :: ierr
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


  subroutine couple_with_injection_prepare_boundary_fwat(evtid)
    character(len=*), intent(in) :: evtid
    character(len=MAX_STRING_LEN)  :: out_dir, fkprname
    logical :: findfile
    integer :: ier, FID=858

    !! for FK point for intialization injected wavefield
    ! real(kind=cr) :: Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box
    real(kind=cr) :: ray_p,Tg,DF_FK

    real(kind=cr), parameter :: TOL_ZERO_TAKEOFF = 1.e-14
    
    out_dir = trim(local_path_backup)//'/FK_wavefield_'//trim(evtid)//'/'
    ! call system('mkdir -p '//trim(out_dir))
    call mkdir(out_dir)
    write(fkprname,'(a,i6.6,a)') trim(out_dir)//'proc', worldrank, '_fk_wavefield.bin'

    call read_fk_model(evtid)
    ! checks if anything to do
    ! for forward simulation only
    ! user output
    ! if (myrank == 0) then
    !   write(IMAIN,*) "preparing injection boundary"
    !   call flush_IMAIN()
    ! endif

    ! FK boundary
    ! initial setup for future FK3D calculations
    ! get MPI starting time for FK
    ! tstart = wtime()

    ! user output
    ! if (myrank == 0) then
    !   write(IMAIN,*) "  using FK injection technique"
    !   call flush_IMAIN()
    ! endif

    ! call FindBoundaryBox(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)
    ! call ReadFKModelInput(Xmin_box, Xmax_box, Ymin_box, Ymax_box, Zmin_box, Zmax_box)

    ! send FK parameters to others MPI slices
    ! call bcast_all_singlei(type_kpsv_fk)
    ! call bcast_all_singlei(nlayer)
    
    ! if (myrank > 0) then
    !   allocate(alpha_FK(nlayer), &
    !             beta_FK(nlayer), &
    !             rho_FK(nlayer), &
    !             mu_FK(nlayer), &
    !             h_FK(nlayer),stat=ier)
    !   if (ier /= 0) call exit_MPI_without_rank('error allocating arrays 2206')
    !   alpha_FK(:) = 0._cr; beta_FK(:) = 0._cr; rho_FK(:) = 0._cr
    !   mu_FK(:) = 0._cr; h_FK(:) = 0._cr
    ! endif

    ! call bcast_all_cr(alpha_FK, nlayer)
    ! call bcast_all_cr(beta_FK, nlayer)
    ! call bcast_all_cr(rho_FK, nlayer)
    ! call bcast_all_cr(mu_FK, nlayer)
    ! call bcast_all_cr(h_FK, nlayer)

    ! call bcast_all_singlecr(phi_FK)
    ! call bcast_all_singlecr(theta_FK)

    ! call bcast_all_singlecr(ff0)
    ! call bcast_all_singlecr(freq_sampling_fk)
    ! call bcast_all_singlecr(amplitude_fk)

    ! call bcast_all_singlecr(xx0)
    ! call bcast_all_singlecr(yy0)
    ! call bcast_all_singlecr(zz0)
    ! call bcast_all_singlecr(Z_REF_for_FK)

    ! call bcast_all_singlecr(tt0)
    ! call bcast_all_singlecr(tmax_fk)

    ! converts origin point Z to reference framework depth for FK,
    ! where top of lower half-space has to be at z==0
    zz0 = zz0 - Z_REF_for_FK

    ! converts to rad
    phi_FK   = phi_FK * PI/180.d0    ! azimuth
    theta_FK = theta_FK * PI/180.d0  ! take-off

    ! ray parameter p (according to Snell's law: sin(theta1)/v1 == sin(theta2)/v2)
    if (type_kpsv_fk == 1) then
      ! P-wave
      ray_p = sin(theta_FK)/alpha_FK(nlayer)    ! for vp (i.e., alpha)
    else if (type_kpsv_fk == 2) then
      ! SV-wave
      ray_p = sin(theta_FK)/beta_FK(nlayer)     ! for vs (i.e., beta)
    endif

    ! note: vertical incident (theta==0 -> p==0) is not handled.
    !       here, it limits ray parameter p to a very small value to handle the calculations
    if (abs(ray_p) < TOL_ZERO_TAKEOFF) ray_p = sign(TOL_ZERO_TAKEOFF,ray_p)

    ! maximum period
    Tg  = 1.d0 / ff0

    ! counts total number of (local) GLL points on absorbing boundary
    call count_num_boundary_points(num_abs_boundary_faces,abs_boundary_ispec,npt)

    !! compute the bottom midle point of the domain

    !! VM VM dealocate in case of severals runs occurs in inverse_problem program
    if (allocated(ipt_table)) deallocate(ipt_table)
    if (allocated(Veloc_FK))  deallocate(Veloc_FK)
    if (allocated(Tract_FK))  deallocate(Tract_FK)

    !! allocate memory for FK solution
    if (npt > 0) then
      allocate(ipt_table(NGLLSQUARE,num_abs_boundary_faces), stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2202')
    else
      ! dummy
      allocate(ipt_table(1,1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2204')
    endif
    ipt_table(:,:) = 0

    deltat = real(DT, cr)
    call find_size_of_working_arrays(deltat, freq_sampling_fk, tmax_fk, NF_FOR_STORING, &
                                      NF_FOR_FFT, NPOW_FOR_INTERP, NP_RESAMP, DF_FK)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  computed FK parameters:'
      write(IMAIN,*) '    frequency sampling rate        = ', freq_sampling_fk,"(Hz)"
      write(IMAIN,*) '    number of frequencies to store = ', NF_FOR_STORING
      write(IMAIN,*) '    number of frequencies for FFT  = ', NF_FOR_FFT
      write(IMAIN,*) '    power of 2 for FFT             = ', NPOW_FOR_INTERP
      write(IMAIN,*)
      write(IMAIN,*) '    simulation time step           = ', deltat,"(s)"
      write(IMAIN,*) '    total simulation length        = ', NSTEP*deltat,"(s)"
      write(IMAIN,*)
      write(IMAIN,*) '    FK time resampling rate        = ', NP_RESAMP
      write(IMAIN,*) '    new time step for F-K          = ', NP_RESAMP * deltat,"(s)"
      write(IMAIN,*) '    new time window length         = ', tmax_fk,"(s)"
      write(IMAIN,*)
      write(IMAIN,*) '    frequency step for F-K         = ', DF_FK,"(Hz)"
      write(IMAIN,*)
      write(IMAIN,*) '  total number of points on boundary = ',npt
      call flush_IMAIN()
    endif

    ! safety check with number of simulation time steps
    if (NSTEP/NP_RESAMP > NF_FOR_STORING + NP_RESAMP) then
      if (myrank == 0) then
        print *,'Error: FK time window length ',tmax_fk,' and NF_for_storing ',NF_FOR_STORING
        print *,'       are too small for chosen simulation length with NSTEP = ',NSTEP
        print *
        print *,'       you could use a smaller NSTEP <= ',NF_FOR_STORING*NP_RESAMP
        print *,'       or'
        print *,'       increase FK window length larger than ',(NSTEP/NP_RESAMP - NP_RESAMP) * NP_RESAMP * deltat
        print *,'       to have a NF for storing  larger than ',(NSTEP/NP_RESAMP - NP_RESAMP)
      endif
      stop 'Invalid FK setting'
    endif

    ! safety check
    if (NP_RESAMP == 0) then
      if (myrank == 0) then
        print *,'Error: FK resampling rate ',NP_RESAMP,' is invalid for frequency sampling rate ',freq_sampling_fk
        print *,'       and the chosen simulation DT = ',DT
        print *
        print *,'       you could use a higher frequency sampling rate>',1./(deltat)
        print *,'       (or increase the time stepping size DT if possible)'
      endif
      stop 'Invalid FK setting'
    endif

    ! limits resampling sizes
    if (NP_RESAMP > 10000) then
      if (myrank == 0) then
        print *,'Error: FK resampling rate ',NP_RESAMP,' is too high for frequency sampling rate ',freq_sampling_fk
        print *,'       and the chosen simulation DT = ',deltat
        print *
        print *,'       you could use a higher frequency sampling rate>',1./(10000*deltat)
        print *,'       (or increase the time stepping size DT if possible)'
      endif
      stop 'Invalid FK setting'
    endif

    if (npt > 0) then
      !! arrays for storing FK solution --------------------------------------------
      allocate(Veloc_FK(NDIM, npt, -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2210')
      if (ier /= 0) stop 'error while allocating Veloc_FK'
      Veloc_FK(:,:,:) = 0._CUSTOM_REAL

      allocate(Tract_FK(NDIM, npt, -NP_RESAMP:NF_FOR_STORING+NP_RESAMP),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2210')
      if (ier /= 0) stop 'error while allocating Veloc_FK'
      Tract_FK(:,:,:) = 0._CUSTOM_REAL

      call FK3D(type_kpsv_fk, nlayer, NSTEP, npt, &
                ray_p, phi_FK, xx0, yy0, zz0, Tg, &
                tt0, alpha_FK, beta_FK, rho_FK, h_FK, &
                NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, DF_FK)
    endif

    ! MX: write FK solution to file
    ! open(FID, file=fkprname, form='unformatted', status='unknown', action='write', iostat=ier)
    ! if (ier /= 0) call exit_MPI(worldrank, 'error opening file 2205')
    ! write(FID) Veloc_FK, Tract_FK
    ! close(FID)

    call write_fk_coupling_file(evtid)

    call synchronize_all()

    ! get MPI ending time for FK
    ! tCPU = wtime() - tstart

    ! ! user output
    ! if (myrank == 0) then
    !     write(IMAIN,'(a35,1x, f20.2, a7)')  " Elapsed time for FK computation : ",  tCPU, " sec. "
    !   write(IMAIN,*)
    !   call flush_IMAIN()
    ! endif

    call free_fk_arrays()

    ! * end of initial setup for future FK3D calculations *

  end subroutine couple_with_injection_prepare_boundary_fwat

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
    integer :: ier, FID=858
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
    real(kind=cr) :: p, z0, zi, eta, theta_rad, phi_rad
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

    z0 = zz0 - Z_ref_for_FK
    zi = zz - Z_ref_for_FK
    ! h = zeros(nlayer)
    allocate(h(nlayer))
    h = 0._cr
    ilayer = nlayer
    do j = nlayer - 1, 1, -1
      if (zi <= sum(h_FK(j:nlayer))) then
        ilayer = j
        exit
      end if
    end do
    h(ilayer+1:nlayer) = h_FK(ilayer+1:nlayer)
    h(ilayer) = zi - sum(h_FK(ilayer+1:nlayer))
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