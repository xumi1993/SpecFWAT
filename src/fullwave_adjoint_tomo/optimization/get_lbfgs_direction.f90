! This program is used to compute lbfgs update direction
! Author: Kai Wang, wangkaim8@gmail.com
! University of Toronto, Ontario, CA
! Last modified: Dec 23, 2018
! 
! 2018/12/23, Kai, revised from compute_direction_lbfgs.f90
!     written by Hejun Zhu
!=====================================================================

  subroutine get_lbfgs_direction_iso()

! calculates TI gradient based on a conjugate gradient method
!
! based on: Tarantola, Inverse problem theory, 2005.
!                  section 6.22.7 conjugate directions, page 217.
!                  formula for alpha_n based on Polak & Ribiere (1969)
!
! note: we use a preconditioner F_0 = 1, thus lambda_n = gamma_n in (6.322)
!          and use gamma_n as the smoothed kernel (for bulk_c, bulk_betav,..).
!
!          however, one could see smoothing as preconditioner F_0, thus
!          gamma_n would be un-smoothed kernel and lambda_n would be smoothed one...
!          i'm not sure if this makes a difference.

  use tomography_par
  use fwat_input, only: tomo_par
  use tomography_kernels_iso
  !!! Kai added
  use specfem_par, only: jacobian,wxgll, wygll, wzgll
  !!!
  implicit none
  real(kind=CUSTOM_REAL) :: minmax(4),depthmax(2),depthmax_radius(2),max
  real(kind=CUSTOM_REAL) :: r,rmax_vs,depthmax_depth
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_beta,norm_rho
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_beta_sum,norm_rho_sum
  real(kind=CUSTOM_REAL) :: min_vp,max_vp,min_vs,max_vs,min_rho,max_rho
  integer :: maxindex(1)
  !!! Kai added !!!
  ! L-BFGS
  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB) :: q_vector,r_vector
  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB) :: gradient1,gradient0,model1,model0
  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB) :: gradient_diff,model_diff
  real(kind=CUSTOM_REAL),dimension(128) :: p(0:1000),a(0:1000)
  real(kind=CUSTOM_REAL) :: p_tmp,p_sum,a_tmp,a_sum,b_tmp,b_sum
  real(kind=CUSTOM_REAL) :: b,p_k_up,p_k_down,p_k_up_sum,p_k_down_sum,p_k
  integer :: iglob
  integer :: i,j,k,ispec,ier
  integer :: iter_store,istore
  real(kind=CUSTOM_REAL),dimension(NKERNEL,NGLOB)::vector_gll
  
  ! jacobian
  ! See use sepcfem_par
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: jacobian
  real(kind=CUSTOM_REAL) :: jacobian_regular!BinHe
  integer, dimension(:), allocatable :: irregular_element_number!BinHe
  !real(kind=CUSTOM_REAL) :: volumel
  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll
  !real(kind=CUSTOM_REAL) :: jacobianl
  integer :: ival,NSPEC_IRREGULAR,s1_jac,s2_jac,s3_jac,s4_jac!BinHe
  character(len=MAX_STRING_LEN) :: m_file!BinHe
  ! if normalization is applied 
  logical :: USE_NORM=.false.
  !!!
  
 

  ! allocate arrays for storing gradient
  allocate(model_dbulk(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_dbeta(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_drho(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)

  if (ier /= 0) stop 'error allocating gradient arrays'

  ! initializes arrays
  model_dbulk = 0.0_CUSTOM_REAL
  model_dbeta = 0.0_CUSTOM_REAL
  model_drho = 0.0_CUSTOM_REAL
  !******* The following is modified from compute_kernel_integral.f90 ***********
  if(allocated(jacobian)) deallocate(jacobian)
  !allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  !if (ier /= 0) stop 'Error allocating jacobian array'
  ! reads NSPEC_IRREGULAR
  write(m_file,'(a,i6.6,a)') trim(INPUT_DATABASES_DIR)//'proc',myrank,trim(REG)//'external_mesh.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif

  read(IIN) ival !nspec
  if (ival /= nspec) call exit_mpi(myrank,'Error invalid nspec value in external_mesh.bin')
  read(IIN) ival !nglob
  if (ival /= nglob) call exit_mpi(myrank,'Error invalid nspec value in external_mesh.bin')

  read(IIN) NSPEC_IRREGULAR
  close(IIN)

  !allocations
  if (NSPEC_IRREGULAR > 0) then
    allocate(jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1064')
    s1_jac = NGLLX
    s2_jac = NGLLY
    s3_jac = NGLLZ
    s4_jac = NSPEC_IRREGULAR
  else
    allocate(jacobian(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1065')
    s1_jac = 1
    s2_jac = 1
    s3_jac = 1
    s4_jac = 1
  endif
  if (ier /= 0) then
    print *,'Error allocating array jacobian'
    call exit_mpi(myrank,'error allocation jacobian')
  endif

  allocate(irregular_element_number(NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1066')


  ! GLL points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  ! builds jacobian
  !call compute_jacobian(jacobian)
  call compute_jacobian(jacobian,irregular_element_number,jacobian_regular,s1_jac,s2_jac,s3_jac,s4_jac)
  !******* The following is modified from compute_direction_lbfgs.f90 ************
  iter_store = iter_current-tomo_par%LBFGS_M_STORE
  if ( iter_store <= iter_start ) then
        iter_store = iter_start
  endif
  if (myrank == 0) print *, 'stored iteration:',iter_store

  !!! TODO: read kernel name from postproc
  nkernel_name(1) = "alpha_kernel_smooth"
  nkernel_name(2) = "beta_kernel_smooth"
  nkernel_name(3) = "rhop_kernel_smooth"
  nmodel_name(1) = "vp"
  nmodel_name(2) = "vs"
  nmodel_name(3) = "rho"


  ! initialize arrays
  a(:)=0.0
  p(:)=0.0
  gradient1(:)=0.0
  gradient0(:)=0.0
  model1(:)=0.0
  model0(:)=0.0
  gradient_diff(:)=0.0
  model_diff(:)=0.0
  q_vector(:)=0.0
  r_vector(:)=0.0

  call get_gradient(iter_current,q_vector)  ! q<--- g_k

  if (myrank == 0) then
     print *,'************************************************'
     print *,'*******starting backward store *****************'
     print *,'************************************************'
  endif

  do istore=iter_current-1,iter_store,-1
     call get_gradient(istore+1,gradient1)
     call get_gradient(istore,gradient0)
     call get_gll_model(istore+1,model1)
     call get_gll_model(istore,model0)
     gradient_diff=gradient1-gradient0
     model_diff=model1-model0

     if (USE_NORM) then
       call Parallel_ComputeInnerProduct(gradient_diff, model_diff, NKERNEL, p_sum) 
     else
       p_tmp=sum(gradient_diff*model_diff)
       call sum_all_all_cr(p_tmp,p_sum)
     endif
     if (abs(p_sum) < 1.e-22) call exit_mpi(myrank,'p vector is zero')
     p(istore) = 1.0_CUSTOM_REAL / p_sum !1.0/p_sum

     if (USE_NORM) then
       call Parallel_ComputeInnerProduct(model_diff, q_vector, NKERNEL, a_sum)
     else
       a_tmp=sum(model_diff*q_vector)
       call sum_all_all_cr(a_tmp,a_sum)
     endif
     a(istore)=p(istore)*a_sum    !a_i <--- p_i*s_i^T*q

     if (myrank == 0) print *,'a,p:',a(istore),p(istore)
     q_vector=q_vector-a(istore)*gradient_diff  !q<--- q-a_i*y_i
  enddo
  istore=iter_current-1
  call get_gradient(istore+1,gradient1)
  call get_gradient(istore,gradient0)
  call get_gll_model(istore+1,model1)
  call get_gll_model(istore,model0)
  gradient_diff=gradient1-gradient0
  model_diff=model1-model0

! this implements Algorithm equation (9.6) on page 226 of the book of
! Jorge Nocedal and Stephen Wright, "Numerical Optimization", Springer (2006)
  
  if (USE_NORM) then
    call Parallel_ComputeInnerProduct(gradient_diff, model_diff, NKERNEL, p_k_up_sum)  
    call Parallel_ComputeL2normSquare(gradient_diff, NKERNEL, p_k_down_sum)  
  else
    p_k_up=sum(gradient_diff*model_diff)
    p_k_down=sum(gradient_diff*gradient_diff)
    call sum_all_all_cr(p_k_up,p_k_up_sum)
    call sum_all_all_cr(p_k_down,p_k_down_sum)
  endif
  p_k=p_k_up_sum/p_k_down_sum

  if ( myrank == 0) print *,'p_k:',p_k
  r_vector=p_k*q_vector
  !r_vector=1.0*q_vector

  if (myrank == 0) then
     print *,'******************************************'
     print *,'********starting forward store ***********'
     print *,'******************************************'
  endif

  do istore=iter_store,iter_current-1,1
     call get_gradient(istore+1,gradient1)
     call get_gradient(istore,gradient0)
     call get_gll_model(istore+1,model1)
     call get_gll_model(istore,model0)

     gradient_diff=gradient1-gradient0
     model_diff=model1-model0

     
     if (USE_NORM) then
       call Parallel_ComputeInnerProduct(gradient_diff, r_vector, NKERNEL, b_sum)
     else
       b_tmp=sum(gradient_diff*r_vector)
       call sum_all_all_cr(b_tmp,b_sum)
     endif
     b=p(istore)*b_sum

     if (myrank == 0) print *,'a,b:',a(istore),b

     r_vector=r_vector+model_diff*(a(istore)-b)

  enddo
  r_vector=-1.0*r_vector


   ! initializes kernel maximum
  depthmax(:) = 0._CUSTOM_REAL

  vector_gll(1,1:NGLOB) = r_vector(1:NGLOB)
  vector_gll(2,1:NGLOB) = r_vector(NGLOB+1:2*NGLOB)
  vector_gll(3,1:NGLOB) = r_vector(2*NGLOB+1:3*NGLOB)
   
  do ispec=1,NSPEC
     do k=1,NGLLZ
        do j=1,NGLLY
           do i=1,NGLLX
              iglob=ibool(i,j,k,ispec)
              model_dbulk(i,j,k,ispec) = vector_gll(1,iglob)
              model_dbeta(i,j,k,ispec) = vector_gll(2,iglob)
              model_drho(i,j,k,ispec)  = vector_gll(3,iglob)
             ! determines maximum kernel betav value within given radius
              if (USE_DEPTH_RANGE_MAXIMUM) then
                r = z(iglob)

                ! stores maximum kernel betav/betah value in this depth slice,
                ! since betav/betah are most likely dominating
                if (r < R_TOP .and. r > R_BOTTOM) then
                  ! kernel betav value
                  max_vs = abs( model_dbeta(i,j,k,ispec) )
                  if (depthmax(1) < max_vs) then
                    depthmax(1) = max_vs
                    depthmax_radius(1) = r
                  endif
                endif
              endif

           enddo
         enddo
      enddo
  enddo

!********************************************************
  ! stores model_dbulk, ... arrays
  ! note: stores these new gradient before we scale them with the step length
  call write_gradient_iso()

  ! statistics
  call min_all_cr(minval(model_dbulk),min_vp)
  call max_all_cr(maxval(model_dbulk),max_vp)

  call min_all_cr(minval(model_dbeta),min_vs)
  call max_all_cr(maxval(model_dbeta),max_vs)

  call min_all_cr(minval(model_drho),min_rho)
  call max_all_cr(maxval(model_drho),max_rho)

  if (myrank == 0) then
    print *,'initial gradient updates:'
    print *,'  vp min/max : ',min_vp,max_vp
    print *,'  vs min/max: ',min_vs,max_vs
    print *,'  rho min/max  : ',min_rho,max_rho
    print *
  endif

  ! determines maximum kernel betav value within given radius
  if (USE_DEPTH_RANGE_MAXIMUM) then
    ! maximum of all processes stored in max_vs
    call max_all_cr(depthmax(1),max_vs)
    call max_all_cr(depthmax_radius(1),rmax_vs)
  endif

  ! determines step length
  ! based on maximum gradient value (either vsv or vsh)
  if (myrank == 0) then

    ! determines maximum kernel betav value within given radius
    if (USE_DEPTH_RANGE_MAXIMUM) then
      depthmax(1) = max_vs
      depthmax_radius(1) = rmax_vs

      max = maxval(depthmax)
      maxindex = maxloc(depthmax)
      depthmax_depth = depthmax_radius(maxindex(1))
      ! maximum in given depth range
      print *,'  using depth maximum: '
      print *,'  between depths (top/bottom)   : ',R_TOP,R_BOTTOM
      print *,'  maximum kernel value          : ',max
      print *,'  depth of maximum kernel value : ',depthmax_depth
      print *
    else
      ! maximum gradient values
      minmax(1) = abs(min_vs)
      minmax(2) = abs(max_vs)
      minmax(3) = abs(max_vp)
      minmax(4) = abs(max_vp)

      ! maximum value of all kernel maxima
      max = maxval(minmax)
    endif
    print *,'step length:'
    print *,'  using kernel maximum: ',max

    ! checks maximum value
    if (max < 1.e-25) stop 'Error maximum kernel value too small for update'

    ! chooses step length such that it becomes the desired, given step factor as inputted
    step_length = step_fac/max

    print *,'  step length : ',step_length
    print *

  endif
  call bcast_all_singlecr(step_length)


  ! gradient length sqrt( v^T * v )
  norm_bulk = sum( model_dbulk * model_dbulk )
  norm_beta = sum( model_dbeta * model_dbeta )
  norm_rho = sum( model_drho * model_drho )

  call sum_all_cr(norm_bulk,norm_bulk_sum)
  call sum_all_cr(norm_beta,norm_beta_sum)
  call sum_all_cr(norm_rho,norm_rho_sum)

  if (myrank == 0) then
    norm_bulk = sqrt(norm_bulk_sum)
    norm_beta = sqrt(norm_beta_sum)
    norm_rho = sqrt(norm_rho_sum)

    print *,'norm model updates:'
    print *,'  bulk : ',norm_bulk
    print *,'  beta: ',norm_beta
    print *,'  rho  : ',norm_rho
    print *
  endif

  ! multiply model updates by a subjective factor that will change the step
  model_dbulk(:,:,:,:) = step_length * model_dbulk(:,:,:,:)
  model_dbeta(:,:,:,:) = step_length * model_dbeta(:,:,:,:)
  model_drho(:,:,:,:) = step_length * model_drho(:,:,:,:)

  ! statistics
  call min_all_cr(minval(model_dbulk),min_vp)
  call max_all_cr(maxval(model_dbulk),max_vs)

  call min_all_cr(minval(model_dbeta),min_vs)
  call max_all_cr(maxval(model_dbeta),max_vs)

  call min_all_cr(minval(model_drho),min_rho)
  call max_all_cr(maxval(model_drho),max_rho)


  deallocate(jacobian)
  deallocate(irregular_element_number)!BinHe

  if (myrank == 0) then
    print *,'scaled gradient:'
    print *,'  vp min/max : ',min_vp,max_vp
    print *,'  vs min/max: ',min_vs,max_vs
    print *,'  eta min/max  : ',min_rho,max_rho
    print *
  endif
  call synchronize_all()

  end subroutine get_lbfgs_direction_iso

! The following subroutine modified from compute_direction_lbfgs.f90 
!=================================================================
subroutine get_gradient(iter,gradient)
  !use globe_parameter
  use tomography_par ! Kai modified 
  ! use fullwave_adjoint_tomo_par, only: USE_RHO_SCALING_FWAT
  use fwat_input, only: tomo_par

  implicit none
  integer::iter,iglob
  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB)::gradient
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC)::vector
  real(kind=CUSTOM_REAL),dimension(NKERNEL,NGLOB)::vector_gll
  integer::ispec,i,j,k,ier ! Kai modefied

  do iker=1,NKERNEL
     write(dirname,'(a,i2.2)') 'optimize/SUM_KERNELS_M',iter
     !write(dirname,'(a,i2.2)') 'SUM_KERNELS_M',iter
     write(filename,'(a,i6.6,a)') trim(dirname)//'/proc',myrank,'_'//trim(nkernel_name(iker))//'.bin'
     open(1001,file=trim(filename),status='old',form='unformatted',iostat=ier)
     if ( myrank == 0) print *,'reading gradient:',trim(filename)
     if (ier /= 0 ) then
        print *,'error reading:',trim(filename)
        call exit_mpi(myrank,'file not found')
     endif
     read(1001) vector(:,:,:,1:NSPEC)
     close(1001)
     do ispec=1,NSPEC
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 iglob=ibool(i,j,k,ispec)
                 vector_gll(iker,iglob)=vector(i,j,k,ispec)
              enddo
            enddo
         enddo
     enddo
  enddo
  if (tomo_par%USE_RHO_SCALING_FWAT) then
     vector_gll(3,:)=RHO_SCALING * vector_gll(2,:)
  endif
  gradient(1:NGLOB)=vector_gll(1,1:NGLOB)
  gradient(NGLOB+1:2*NGLOB)=vector_gll(2,1:NGLOB)
  gradient(2*NGLOB+1:3*NGLOB)=vector_gll(3,1:NGLOB)
  !!!
end subroutine get_gradient


subroutine get_gll_model(iter,model)
  !use globe_parameter
  use tomography_par ! Kai modified 
  implicit none
  integer::iter,iglob
  real(kind=CUSTOM_REAL),dimension(NKERNEL*NGLOB):: model
  real(kind=CUSTOM_REAL),dimension(NKERNEL,NGLOB):: vector_gll
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC):: vector
  integer::ispec,i,j,k,ier ! Kai modefied
  
  do iker=1,NKERNEL
     write(dirname,'(a,i2.2)') 'optimize/MODEL_M',iter
     !write(dirname,'(a,i2.2)') 'MODEL_M',iter
     write(filename,'(a,i6.6,a)') trim(dirname)//'/proc',myrank,'_'//trim(nmodel_name(iker))//'.bin'
     open(1001,file=trim(filename),status='old',form='unformatted',iostat=ier)
     if ( myrank == 0) print *,'reading model:',trim(filename)
     if ( ier /= 0) then
        print *,'error reading:',trim(filename)
        call exit_mpi(myrank,'file not found')
     endif
     read(1001) vector(:,:,:,1:NSPEC)
     close(1001)
     do ispec=1,NSPEC
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 iglob=ibool(i,j,k,ispec)
                 vector_gll(iker,iglob)=vector(i,j,k,ispec)
              enddo
           enddo
        enddo
     enddo
  enddo
  model(1:NGLOB)=log(vector_gll(1,1:NGLOB))
  model(NGLOB+1:2*NGLOB)=log(vector_gll(2,1:NGLOB))
  model(2*NGLOB+1:3*NGLOB)=log(vector_gll(3,1:NGLOB))
end subroutine get_gll_model

!====================================================
! The following codes are modified from
! src/inverse_problem_for_model/inversion_scheme/inversion_scheme_mod.f90 
!---------------------------------------------------------------------------------------------

  subroutine Parallel_ComputeInnerProduct(vect1, vect2, Niv, qp)

    use tomography_par
    use specfem_par, only: jacobian, wxgll, wygll, wzgll

    real(kind=CUSTOM_REAL), dimension(NKERNEL*NGLOB)     :: vect1, vect2
    real(kind=CUSTOM_REAL)                               :: qp
    integer                                              :: Niv
    real(kind=CUSTOM_REAL)                                                   :: qp_tmp_single
    ! try double precision
    real(kind=8)                                                             :: jacobianl, weight, qp_tmp
    integer                                                                  :: ipar, i, j, k, ispec
    real(kind=CUSTOM_REAL)                                                   :: coeff, coeff_n1, coeff_n2
    real(kind=8)                                                             :: coeff_n1_dp, coeff_n2_dp
    integer :: iglob

    !! try normalization to avoid numerical errors
    !call Parallel_ComputeL2normSquare(vect1 , Niv, coeff_n1)
    !call Parallel_ComputeL2normSquare(vect2 , Niv, coeff_n2)

    coeff=maxval(abs(vect1(:)))
    call max_all_all_cr(coeff, coeff_n1)
    if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1=1._CUSTOM_REAL
    vect1(:) = vect1(:) / coeff_n1

    coeff=maxval(abs(vect2(:)))
    call max_all_all_cr(coeff, coeff_n2)
    if (coeff_n2 == 0._CUSTOM_REAL) coeff_n2=1._CUSTOM_REAL
    vect2(:) = vect2(:) / coeff_n2

    coeff_n1_dp = coeff_n1
    coeff_n2_dp = coeff_n2

    qp_tmp=0._CUSTOM_REAL

    do ipar=1, Niv
       do ispec = 1, NSPEC
          do k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX
                   iglob=ibool(i,j,k,ispec)
                   weight = wxgll(i)*wygll(j)*wzgll(k)
                   jacobianl = jacobian(i,j,k,ispec)
                   qp_tmp = qp_tmp + jacobianl * weight * vect1(iglob+(ipar-1)*NGLOB) * vect2(iglob+(ipar-1)*NGLOB)
                   !qp = qp + jacobianl * weight * vect1(i,j,k,ispec,ipar) * vect2(i,j,k,ispec,ipar)
                enddo
             enddo
          enddo
       enddo
    enddo

    qp_tmp_single = qp_tmp * coeff_n1_dp * coeff_n2_dp
    qp=0.
    call sum_all_all_cr(qp_tmp_single, qp)

  end subroutine Parallel_ComputeInnerProduct
  !---------------------------------------------------------------------------------------------


   subroutine Parallel_ComputeL2normSquare(vect1 , Niv, qp)
 
     use tomography_par
     use specfem_par, only: jacobian, wxgll, wygll, wzgll

     real(kind=CUSTOM_REAL), dimension(NKERNEL*NGLOB)           :: vect1
     real(kind=CUSTOM_REAL)                                     :: qp
     integer                                                    :: Niv
     real(kind=CUSTOM_REAL) :: coeff, coeff_n1
     real(kind=CUSTOM_REAL) :: qp_tmp
     real(kind=8) :: jacobianl, weight, qp_dp, coeff_n1_dp
     integer :: ipar, i, j, k, ispec,iglob

     qp=0.d0
     qp_dp=0.d0

     coeff=maxval(abs(vect1(:)))
     call max_all_all_cr(coeff, coeff_n1)

     if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1=1._CUSTOM_REAL

     vect1(:) = vect1(:) / coeff_n1
     coeff_n1_dp=coeff_n1

     do ipar=1,Niv
        do ispec = 1, NSPEC
           do k=1,NGLLZ
              do j=1,NGLLY
                 do i=1,NGLLX
                    iglob=ibool(i,j,k,ispec)
                    weight = wxgll(i)*wygll(j)*wzgll(k)
                    jacobianl = jacobian(i,j,k,ispec)
                    qp_dp = qp_dp + jacobianl * weight * vect1(iglob+(ipar-1)*NGLOB) **2
                 enddo
              enddo
           enddo
        enddo
     enddo

     qp_tmp = qp_dp * coeff_n1_dp * coeff_n1_dp
     qp=0.
     call sum_all_all_cr(qp_tmp, qp)

  end subroutine Parallel_ComputeL2normSquare
  !---------------------------------------------------------------------------------------------
