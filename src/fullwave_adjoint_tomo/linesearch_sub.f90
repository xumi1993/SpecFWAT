module linesearch_sub
  use constants, only: MAX_STRING_LEN
  use fwat_input
  use fwat_utils
  use postproc_sub, only: read_database,model,set_range,read_misfit,type_name
  implicit none

  
contains
  subroutine get_wolfe_amijo(qp, q0)
    use tomography_kernels_iso

    real(CUSTOM_REAL) :: ,qp
    character(len=MAX_STRING_LEN) :: fname, model
    real(kind=CUSTOM_REAL),dimension(:),allocatable :: direction_gll,&
        direction,d_vector,gradient_gll, gradient,g_vector
    double precision :: mean_chi
    integer :: i, j,k,iker,nchan

    integer :: iker

    call read_database() ! from postproc_sub, read_mesh_databases

    allocate(gradient_gll(NKERNEL,NGLOB),gradient(NKERNEL*NGLOB))
    allocate(direction_gll(NKERNEL,NGLOB),direction(NKERNEL*NGLOB))
    allocate(g_vector(NGLLX,NGLLY,NGLLZ,NSPEC))
    allocate(d_vector(NGLLX,NGLLY,NGLLZ,NSPEC))

    ! read gradient 
    call read_kernels_iso()

    ! reads in smoothed kernels: bulk, beta, rho
    if (trim(tomo_par%OPT_METHOD)=='CG'.and. iter_current-iter_start>0) then
      call read_kernels_cg_iso_old()
    endif

    ! calculates gradient
    ! steepest descent method
    if (trim(tomo_par%OPT_METHOD)=='SD' .or. iter_current-iter_start==0) then
      call get_sd_direction_iso()
    elseif (trim(tomo_par%OPT_METHOD)=='CG') then 
      call get_cg_direction_iso()
    else
      call get_lbfgs_direction_iso()
    endif
    do iker=1,NKERNEL
      if (iker == 1) then
        d_vector = model_dbulk
        g_vector = kernel_bulk
      endif
      if (iker == 2) then
        d_vector = model_dbeta
        g_vector = kernel_beta
      endif
      if (iker == 3) then
        d_vector = model_drho
        g_vector = kernel_rho
      endif
      do ispec=1,NSPEC
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              iglob=ibool(i,j,k,ispec)
              direction_gll(iker,iglob)=d_vector(i,j,k,ispec)
              gradient_gll(iker,iglob)=g_vector(i,j,k,ispec)
            enddo
          enddo
        enddo
      enddo
    enddo 
    direction(1:NGLOB)=direction_gll(1,1:NGLOB)
    direction(NGLOB+1:2*NGLOB)=direction_gll(2,1:NGLOB)
    direction(2*NGLOB+1:3*NGLOB)=direction_gll(3,1:NGLOB)
    gradient(1:NGLOB)=gradient_gll(1,1:NGLOB)
    gradient(NGLOB+1:2*NGLOB)=gradient_gll(2,1:NGLOB)
    gradient(2*NGLOB+1:3*NGLOB)=gradient_gll(3,1:NGLOB)
    call Parallel_ComputeInnerProduct(gradient,direction,NKERNEL,qp)
    call read_misfit(type_name, model, mean_chi, nchan)
    q0 = mean_chi*nchan

  end subroutine 

end module