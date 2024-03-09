!=====================================================================
!
!               Full Waveform Adjoint Tomography -v1.1
!               ---------------------------------------
!
!     Main historical authors: Kai Wang
!                              Macquarie Uni, Australia
!                            & University of Toronto, Canada
!                           (c) Martch 2020
!
!=====================================================================
!
module specfem_interface


  !! IMPORT SPECFEM ALL VARAIBLES
  use specfem_par
  use specfem_par_coupling
  use kdtree_search, only: kdtree_delete,kdtree_nodes_location,kdtree_nodes_index
  use shared_parameters
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par

  !! IMPORT inverse_problem VARIABLES
  use fullwave_adjoint_tomo_par
  use fwat_input

  !! IMPORT inverse problem modules
  !use adjoint_source
  !use input_output
  !use signal_processing


  implicit none


  !! Arrays for saving GPU kernels because at each GPU run the kernels are set to 0 and to perform
  !! summation over events we need to use those temporarry arrays
  ! real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable :: rho_ac_kl_GPU, kappa_ac_kl_GPU
  ! real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable :: rho_kl_GPU, kappa_kl_GPU, mu_kl_GPU
  ! real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable :: hess_rho_ac_kl_GPU, hess_kappa_ac_kl_GPU
  ! real(kind=CUSTOM_REAL), private, dimension(:,:,:,:),   allocatable :: hess_rho_kl_GPU, hess_kappa_kl_GPU, hess_mu_kl_GPU
  ! real(kind=CUSTOM_REAL), private, dimension(:,:,:,:,:), allocatable :: cijkl_kl_GPU

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------------------------------------------------------------------
!> initialize specfem before each call for the event ievent
!----------------------------------------------------------------------------------------------------------------------------------
  subroutine InitSpecfem() !(acqui_simu, ievent, inversion_param, iter_inverse)

    !integer,                                        intent(in)    ::  ievent, iter_inverse
    !type(inver),                                    intent(in)    ::  inversion_param
    !type(acqui),  dimension(:), allocatable,        intent(in)    ::  acqui_simu

    !integer                                                       :: irec,isrc,ier
    !integer                                                       :: icomp, it, irec_local
    !real(kind=CUSTOM_REAL)                                        :: DT, lw_tap
    !double precision                                              :: DT_dble
    !character(len=256)                                            :: name_file
    !character(len=MAX_STRING_LEN)                                 :: TRAC_PATH, dsname
    !integer(kind=8)                                               :: filesize
    !real(kind=CUSTOM_REAL), dimension(:), allocatable             :: raw_stf, filt_stf
    !character(len=MAX_LEN_STRING)                                 :: name_file_tmp, ch_to_add


    ! de-allocation for setup_search_kdtree() 
    if(allocated(kdtree_nodes_location)) deallocate(kdtree_nodes_location) 
    if(allocated(kdtree_nodes_index)) deallocate(kdtree_nodes_index)

    ! de-allocation for setup_sources_wk(source_fname) 
    if(allocated(islice_selected_source)) deallocate(islice_selected_source)
    if(allocated(ispec_selected_source)) deallocate(ispec_selected_source)
    if(allocated(Mxx)) deallocate(Mxx)
    if(allocated(Myy)) deallocate(Myy)
    if(allocated(Mzz)) deallocate(Mzz)
    if(allocated(Mxy)) deallocate(Mxy)
    if(allocated(Mxz)) deallocate(Mxz)
    if(allocated(Myz)) deallocate(Myz)
    if(allocated(xi_source)) deallocate(xi_source)
    if(allocated(eta_source)) deallocate(eta_source)
    if(allocated(gamma_source)) deallocate(gamma_source)
    if(allocated(tshift_src)) deallocate(tshift_src)
    if(allocated(hdur)) deallocate(hdur)
    if(allocated(hdur_Gaussian)) deallocate(hdur_Gaussian)
    if(allocated(utm_x_source)) deallocate(utm_x_source)
    if(allocated(utm_y_source)) deallocate(utm_y_source)
    if(allocated(nu_source)) deallocate(nu_source)

    if(allocated(factor_force_source)) deallocate(factor_force_source)
    if(allocated(comp_dir_vect_source_E)) deallocate(comp_dir_vect_source_E)
    if(allocated(comp_dir_vect_source_N)) deallocate(comp_dir_vect_source_N)
    if(allocated(comp_dir_vect_source_Z_UP)) deallocate(comp_dir_vect_source_Z_UP)
    if(allocated(user_source_time_function)) deallocate(user_source_time_function)

    if(allocated(run_number_of_the_source)) deallocate(run_number_of_the_source)

    ! de-allocation for setup_sources_precompute_arrays()
    if(allocated(pm1_source_encoding)) deallocate(pm1_source_encoding)
    if(allocated(sourcearray)) deallocate(sourcearray)
    if(allocated(sourcearrays)) deallocate(sourcearrays)
    if(allocated(source_adjoint)) deallocate(source_adjoint)

    !! de-allocation for setup_receivers_wk()
    if (allocated(islice_selected_rec)) deallocate(islice_selected_rec)
    if (allocated(ispec_selected_rec))  deallocate(ispec_selected_rec)
    if (allocated(xi_receiver)) deallocate(xi_receiver)
    if (allocated(eta_receiver)) deallocate(eta_receiver)
    if (allocated(gamma_receiver)) deallocate(gamma_receiver)
    if (allocated(station_name)) deallocate(station_name)
    if (allocated(network_name)) deallocate(network_name)
    if (allocated(nu)) deallocate(nu)

    ! de-allocation for setup_receivers_precompute_intp()
    if (allocated(hxir_store)) deallocate(hxir_store)
    if (allocated(hetar_store)) deallocate(hetar_store)
    if (allocated(hgammar_store)) deallocate(hgammar_store)
    if (allocated(hpxir_store)) deallocate(hpxir_store)
    if (allocated(hpetar_store)) deallocate(hpetar_store)
    if (allocated(hpgammar_store)) deallocate(hpgammar_store)
    if (allocated(number_receiver_global)) deallocate(number_receiver_global)

    if (allocated(seismograms_d)) deallocate(seismograms_d)
    if (allocated(seismograms_v)) deallocate(seismograms_v)
    if (allocated(seismograms_a)) deallocate(seismograms_a)
    if (allocated(seismograms_p)) deallocate(seismograms_p)

    ! de-allocation for prepare_attenuation()
    ! de-allocation for prepare_gravity()
    if(allocated(minus_deriv_gravity)) deallocate(minus_deriv_gravity)
    if(allocated(minus_g)) deallocate(minus_g)
    
    ! de-allocation for prepare_timerun_lddrk() 
    if(allocated(potential_acoustic_lddrk)) deallocate(potential_acoustic_lddrk)
    if(allocated(potential_dot_acoustic_lddrk)) deallocate(potential_dot_acoustic_lddrk)
    if(allocated(displ_lddrk)) deallocate(displ_lddrk)
    if(allocated(veloc_lddrk)) deallocate(veloc_lddrk)
    if(allocated(R_xx_lddrk)) deallocate(R_xx_lddrk)
    if(allocated(R_yy_lddrk)) deallocate(R_yy_lddrk)
    if(allocated(R_xy_lddrk)) deallocate(R_xy_lddrk)
    if(allocated(R_xz_lddrk)) deallocate(R_xz_lddrk)
    if(allocated(R_yz_lddrk)) deallocate(R_yz_lddrk)
    if(allocated(R_trace_lddrk)) deallocate(R_trace_lddrk)
    if(allocated(b_R_xx_lddrk)) deallocate(b_R_xx_lddrk)
    if(allocated(b_R_yy_lddrk)) deallocate(b_R_yy_lddrk)
    if(allocated(b_R_xy_lddrk)) deallocate(b_R_xy_lddrk)
    if(allocated(b_R_xz_lddrk)) deallocate(b_R_xz_lddrk)
    if(allocated(b_R_yz_lddrk)) deallocate(b_R_yz_lddrk)
    if(allocated(b_R_trace_lddrk)) deallocate(b_R_trace_lddrk)
   
    ! de-allocation for prepare_timerun_pml()
    if(allocated(spec_to_CPML)) deallocate(spec_to_CPML)
    if(allocated(CPML_type)) deallocate(CPML_type)
    if(allocated(PML_displ_old)) deallocate(PML_displ_old)
    if(allocated(PML_displ_new)) deallocate(PML_displ_new)
    if(allocated(PML_dux_dxl)) deallocate(PML_dux_dxl)
    if(allocated(PML_dux_dyl)) deallocate(PML_dux_dyl)
    if(allocated(PML_dux_dzl)) deallocate(PML_dux_dzl)
    if(allocated(PML_duy_dxl)) deallocate(PML_duy_dxl)
    if(allocated(PML_duy_dyl)) deallocate(PML_duy_dyl)
    if(allocated(PML_duy_dzl)) deallocate(PML_duy_dzl)
    if(allocated(PML_duz_dxl)) deallocate(PML_duz_dxl)
    if(allocated(PML_duz_dyl)) deallocate(PML_duz_dyl)
    if(allocated(PML_duz_dzl)) deallocate(PML_duz_dzl)
    if(allocated(PML_dux_dxl_old)) deallocate(PML_dux_dxl_old)
    if(allocated(PML_dux_dyl_old)) deallocate(PML_dux_dyl_old)
    if(allocated(PML_dux_dzl_old)) deallocate(PML_dux_dzl_old)
    if(allocated(PML_duy_dxl_old)) deallocate(PML_duy_dxl_old)
    if(allocated(PML_duy_dyl_old)) deallocate(PML_duy_dyl_old)
    if(allocated(PML_duy_dzl_old)) deallocate(PML_duy_dzl_old)
    if(allocated(PML_duz_dxl_old)) deallocate(PML_duz_dxl_old)
    if(allocated(PML_duz_dyl_old)) deallocate(PML_duz_dyl_old)
    if(allocated(PML_duz_dzl_old)) deallocate(PML_duz_dzl_old)
    if(allocated(PML_dux_dxl_new)) deallocate(PML_dux_dxl_new)
    if(allocated(PML_dux_dyl_new)) deallocate(PML_dux_dyl_new)
    if(allocated(PML_dux_dzl_new)) deallocate(PML_dux_dzl_new)
    if(allocated(PML_duy_dxl_new)) deallocate(PML_duy_dxl_new)
    if(allocated(PML_duy_dyl_new)) deallocate(PML_duy_dyl_new)
    if(allocated(PML_duy_dzl_new)) deallocate(PML_duy_dzl_new)
    if(allocated(PML_duz_dxl_new)) deallocate(PML_duz_dxl_new)
    if(allocated(PML_duz_dyl_new)) deallocate(PML_duz_dyl_new)
    if(allocated(PML_duz_dzl_new)) deallocate(PML_duz_dzl_new)
    if(allocated(rmemory_dux_dxl_x)) deallocate(rmemory_dux_dxl_x)
    if(allocated(rmemory_dux_dyl_x)) deallocate(rmemory_dux_dyl_x)
    if(allocated(rmemory_dux_dzl_x)) deallocate(rmemory_dux_dzl_x)
    if(allocated(rmemory_duy_dxl_x)) deallocate(rmemory_duy_dxl_x)
    if(allocated(rmemory_duy_dyl_x)) deallocate(rmemory_duy_dyl_x)
    if(allocated(rmemory_duz_dxl_x)) deallocate(rmemory_duz_dxl_x)
    if(allocated(rmemory_duz_dzl_x)) deallocate(rmemory_duz_dzl_x)
    if(allocated(rmemory_dux_dxl_y)) deallocate(rmemory_dux_dxl_y)
    if(allocated(rmemory_dux_dyl_y)) deallocate(rmemory_dux_dyl_y)
    if(allocated(rmemory_duy_dxl_y)) deallocate(rmemory_duy_dxl_y)
    if(allocated(rmemory_duy_dyl_y)) deallocate(rmemory_duy_dyl_y)
    if(allocated(rmemory_duy_dzl_y)) deallocate(rmemory_duy_dzl_y)
    if(allocated(rmemory_duz_dyl_y)) deallocate(rmemory_duz_dyl_y)
    if(allocated(rmemory_duz_dzl_y)) deallocate(rmemory_duz_dzl_y)
    if(allocated(rmemory_dux_dxl_z)) deallocate(rmemory_dux_dxl_z)
    if(allocated(rmemory_dux_dzl_z)) deallocate(rmemory_dux_dzl_z)
    if(allocated(rmemory_duy_dyl_z)) deallocate(rmemory_duy_dyl_z)
    if(allocated(rmemory_duy_dzl_z)) deallocate(rmemory_duy_dzl_z)
    if(allocated(rmemory_duz_dxl_z)) deallocate(rmemory_duz_dxl_z)
    if(allocated(rmemory_duz_dyl_z)) deallocate(rmemory_duz_dyl_z)
    if(allocated(rmemory_duz_dzl_z)) deallocate(rmemory_duz_dzl_z)
    if(allocated(rmemory_displ_elastic)) deallocate(rmemory_displ_elastic)
    if(allocated(accel_elastic_CPML)) deallocate(accel_elastic_CPML)
    if(allocated(PML_potential_acoustic_old)) deallocate(PML_potential_acoustic_old)
    if(allocated(PML_potential_acoustic_new)) deallocate(PML_potential_acoustic_new)
    if(allocated(rmemory_dpotential_dxl)) deallocate(rmemory_dpotential_dxl)
    if(allocated(rmemory_dpotential_dyl)) deallocate(rmemory_dpotential_dyl)
    if(allocated(rmemory_dpotential_dzl)) deallocate(rmemory_dpotential_dzl)
    if(allocated(rmemory_potential_acoustic)) deallocate(rmemory_potential_acoustic)
    if(allocated(potential_dot_dot_acoustic_CPML)) deallocate(potential_dot_dot_acoustic_CPML)
    if(allocated(rmemory_coupling_ac_el_displ)) deallocate(rmemory_coupling_ac_el_displ)
    if(allocated(rmemory_coupling_el_ac_potential_dot_dot)) deallocate(rmemory_coupling_el_ac_potential_dot_dot)
    if(allocated(rmemory_coupling_el_ac_potential)) deallocate(rmemory_coupling_el_ac_potential)
    if(allocated(b_PML_field)) deallocate(b_PML_field)
    if(allocated(b_PML_potential)) deallocate(b_PML_potential)

    ! de-allocation for prepare_timerun_adjoint()
    if(allocated(Mxx_der)) deallocate(Mxx_der)
    if(allocated(Myy_der)) deallocate(Myy_der)
    if(allocated(Mzz_der)) deallocate(Mzz_der)
    if(allocated(Mxy_der)) deallocate(Mxy_der)
    if(allocated(Mxz_der)) deallocate(Mxz_der)
    if(allocated(Myz_der)) deallocate(Myz_der)
    if(allocated(sloc_der)) deallocate(sloc_der)
    if(allocated(seismograms_eps)) deallocate(seismograms_eps)
    if(allocated(b_absorb_field)) deallocate(b_absorb_field)
    if(allocated(b_absorb_potential)) deallocate(b_absorb_potential)
    if(allocated(b_absorb_fields)) deallocate(b_absorb_fields)
    if(allocated(b_absorb_fieldw)) deallocate(b_absorb_fieldw)

    ! de-allocation for prepare_timerun_noise()
    if(allocated(noise_sourcearray)) deallocate(noise_sourcearray)
    if(allocated(normal_x_noise)) deallocate(normal_x_noise)
    if(allocated(normal_y_noise)) deallocate(normal_y_noise)
    if(allocated(normal_z_noise)) deallocate(normal_z_noise)
    if(allocated(mask_noise)) deallocate(mask_noise)
    if(allocated(noise_surface_movie)) deallocate(noise_surface_movie)


    !! this is to skip writing seismogram on disk by specffem (both forward and ajoint)
    !! do not change
    NTSTEP_BETWEEN_OUTPUT_SEISMOS=NSTEP
    !INVERSE_FWI_FULL_PROBLEM= .true.
    NTSTEP_BETWEEN_READ_ADJSRC=NSTEP

    ! initializes adjoint sources --------------------------------------------------------------------------------------------------
    !if (allocated(source_adjoint)) deallocate(source_adjoint)
    !allocate(source_adjoint(nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,NDIM),stat=ier)
    !if (ier /= 0) stop 'error allocating array source_adjoint'
    !source_adjoint(:,:,:) = 0._CUSTOM_REAL
    !if (SIMULATION_TYPE == 3) then
    !   do icomp=1,NDIM
    !      do it=1,NTSTEP_BETWEEN_OUTPUT_SEISMOS
    !         do irec_local=1, nadj_rec_local
    !            source_adjoint(irec_local, it, icomp) = acqui_simu(ievent)%adjoint_sources_input(icomp,irec_local,it)
    !         enddo
    !      enddo
    !   enddo
    !endif


    ! manage the storage in memory and in disk the simulaed wavefields
    !! with GPU we have to be carrefull of warning for seismogram not all things are allowed and this make the code crashes
    if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION ) then

       !! todo recuperer inversion_paral%get_synthetics_**
       SAVE_SEISMOGRAMS_DISPLACEMENT =.true.
       SAVE_SEISMOGRAMS_VELOCITY     =.false.
       SAVE_SEISMOGRAMS_ACCELERATION =.false.
       SAVE_SEISMOGRAMS_PRESSURE     =.true.

    else

       if (ACOUSTIC_SIMULATION) then
          SAVE_SEISMOGRAMS_PRESSURE     =.true.
          SAVE_SEISMOGRAMS_DISPLACEMENT =.false.
          SAVE_SEISMOGRAMS_VELOCITY     =.false.
          SAVE_SEISMOGRAMS_ACCELERATION =.false.
       endif

       if (ELASTIC_SIMULATION) then
          SAVE_SEISMOGRAMS_PRESSURE      =.false.
          SAVE_SEISMOGRAMS_DISPLACEMENT  =.true.
          SAVE_SEISMOGRAMS_VELOCITY      =.false.
          SAVE_SEISMOGRAMS_ACCELERATION  =.false.
       endif

    endif

    !! reallocate all GPU memory according the Fortran arrays
    !if (GPU_MODE) call prepare_GPU()

    !! open new log file for specfem -----------------------------------------------------------------------------------------------
    !if (myrank == 0 .and. SIMULATION_TYPE == 1) then
    !   close(IMAIN)
    !   write(name_file,'(a15,i5.5,a1,i5.5,a4)') '/output_solver_',ievent,'_',iter_inverse,'.txt'
    !   open(unit=IMAIN,file=trim(OUTPUT_FILES)//trim(name_file),status='unknown')
    !endif

    !! info on mesh and parameters ---------------
    !if (SIMULATION_TYPE == 1) then
    !   if (ELASTIC_SIMULATION) then
    !      call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
    !           ibool,xstore,ystore,zstore, &
    !           kappastore,mustore,rho_vp,rho_vs, &
    !           DT_dble,model_speed_max,min_resolved_period, &
    !           LOCAL_PATH,SAVE_MESH_FILES)
    !   else if (ACOUSTIC_SIMULATION) then
    !      allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    !      if (ier /= 0) stop 'error allocating array rho_vp'
    !      allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    !      if (ier /= 0) stop 'error allocating array rho_vs'
    !      rho_vp = sqrt( kappastore / rhostore ) * rhostore
    !      rho_vs = 0.0_CUSTOM_REAL
    !      call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
    !           ibool,xstore,ystore,zstore, &
    !           kappastore,mustore,rho_vp,rho_vs, &
    !           DT_dble,model_speed_max,min_resolved_period, &
    !           LOCAL_PATH,SAVE_MESH_FILES)
    !      deallocate(rho_vp,rho_vs)
    !   endif
    !endif

  end subroutine InitSpecfem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------------------------------------------------------------------
!> finalize, specfem after running each event
!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine  FinalizeSpecfem() !(acqui_simu, ievent)

    !integer,                                        intent(in)    :: ievent
    !type(acqui),  dimension(:), allocatable,        intent(in)    :: acqui_simu
    integer                                                       :: ier



    !! finalize specfem run --------------------------------------------------------------------------
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
       open(unit=IOUT,file=prname(1:len_trim(prname))//'save_forward_arrays.bin', &
            status='unknown',form='unformatted',iostat=ier)
       if (ier /= 0) then
          print *,'error: opening save_forward_arrays.bin'
          print *,'path: ',prname(1:len_trim(prname))//'save_forward_arrays.bin'
          call exit_mpi(myrank,'error opening file save_forward_arrays.bin')
       endif

       if (ACOUSTIC_SIMULATION) then
          write(IOUT) potential_acoustic
          write(IOUT) potential_dot_acoustic
          write(IOUT) potential_dot_dot_acoustic
       endif

       if (ELASTIC_SIMULATION) then
          write(IOUT) displ
          write(IOUT) veloc
          write(IOUT) accel

          if (ATTENUATION) then
             write(IOUT) R_trace
             write(IOUT) R_xx
             write(IOUT) R_yy
             write(IOUT) R_xy
             write(IOUT) R_xz
             write(IOUT) R_yz
             write(IOUT) epsilondev_trace
             write(IOUT) epsilondev_xx
             write(IOUT) epsilondev_yy
             write(IOUT) epsilondev_xy
             write(IOUT) epsilondev_xz
             write(IOUT) epsilondev_yz
          endif
       endif

       if (POROELASTIC_SIMULATION) then
          write(IOUT) displs_poroelastic
          write(IOUT) velocs_poroelastic
          write(IOUT) accels_poroelastic
          write(IOUT) displw_poroelastic
          write(IOUT) velocw_poroelastic
          write(IOUT) accelw_poroelastic
       endif

       close(IOUT)
    endif

  ! adjoint simulations
  if (SIMULATION_TYPE == 3 .and. tomo_par%SAVE_OUTPUT_EACH_EVENT) then
    ! adjoint kernels
    if (myrank ==0) then 
      write(OUT_FWAT_LOG,*) 'Save event kernel to '//trim(prname)
      flush(OUT_FWAT_LOG)
    endif
    call save_adjoint_kernels_fwat()
  endif

  ! seismograms and source parameter gradients for (pure type=2) adjoint simulation runs
  if (SIMULATION_TYPE == 2) then
    if (nrec_local > 0) then
      ! seismograms (strain)
      call write_adj_seismograms2_to_file(myrank,seismograms_eps,number_receiver_global,nrec_local,it,DT,NSTEP,t0)
      ! source gradients  (for sources in elastic domains)
      call save_kernels_source_derivatives()
    endif
  endif

  ! stacey absorbing fields will be reconstructed for adjoint simulations
  ! using snapshot files of wavefields
  if (STACEY_ABSORBING_CONDITIONS) then
    ! closes absorbing wavefield saved/to-be-saved by forward simulations
    if (num_abs_boundary_faces > 0 .and. (SIMULATION_TYPE == 3 .or. &
          (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

      if (ELASTIC_SIMULATION) call close_file_abs(IOABS)
      if (ACOUSTIC_SIMULATION) call close_file_abs(IOABS_AC)

    endif
  endif

    !if (GPU_MODE .and. SIMULATION_TYPE == 3) then

    !   if (ACOUSTIC_SIMULATION) then
    !      rho_ac_kl_GPU(:,:,:,:)=rho_ac_kl_GPU(:,:,:,:)+rho_ac_kl(:,:,:,:)
    !      kappa_ac_kl_GPU(:,:,:,:)=kappa_ac_kl_GPU(:,:,:,:)+kappa_ac_kl(:,:,:,:)
    !      if (APPROXIMATE_HESS_KL) then
    !         hess_rho_ac_kl_GPU(:,:,:,:)=hess_rho_ac_kl_GPU(:,:,:,:)+hess_rho_ac_kl(:,:,:,:)
    !         hess_kappa_ac_kl_GPU(:,:,:,:)=hess_kappa_ac_kl_GPU(:,:,:,:)+hess_kappa_ac_kl(:,:,:,:)
    !      endif
    !   endif

    !   if (ELASTIC_SIMULATION) then
    !      rho_kl_GPU(:,:,:,:)=rho_kl_GPU(:,:,:,:)+rho_kl(:,:,:,:)
    !      if (ANISOTROPIC_KL) then
    !         cijkl_kl_GPU(:,:,:,:,:)=cijkl_kl_GPU(:,:,:,:,:)+cijkl_kl(:,:,:,:,:)
    !      else
    !         kappa_kl_GPU(:,:,:,:)=kappa_kl_GPU(:,:,:,:)+kappa_kl(:,:,:,:)
    !         mu_kl_GPU(:,:,:,:)= mu_kl_GPU(:,:,:,:)+ mu_kl(:,:,:,:)
    !         if (APPROXIMATE_HESS_KL) then
    !            hess_rho_kl_GPU(:,:,:,:)=hess_rho_kl_GPU(:,:,:,:)+hess_rho_kl(:,:,:,:)
    !            hess_kappa_kl_GPU(:,:,:,:)=hess_kappa_kl_GPU(:,:,:,:)+hess_kappa_kl(:,:,:,:)
    !            hess_mu_kl_GPU(:,:,:,:)= hess_mu_kl_GPU(:,:,:,:)+ hess_mu_kl(:,:,:,:)
    !         endif
    !      endif

    !   endif

    !endif
    ! synchronize all the processes to make sure everybody has finished
    call synchronize_all()
  end subroutine FinalizeSpecfem

end module specfem_interface
