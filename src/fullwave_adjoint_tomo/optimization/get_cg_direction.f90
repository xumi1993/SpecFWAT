!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine get_cg_direction_iso()

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

  use tomography_kernels_iso
  use tomography_kernels_iso_cg
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha_bulk,alpha_beta,alpha_rho,alpha_all
  real(kind=CUSTOM_REAL) :: minmax(4),depthmax(2),depthmax_radius(2),max
  real(kind=CUSTOM_REAL) :: r,rmax_vs,depthmax_depth
  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_beta,norm_rho
  real(kind=CUSTOM_REAL) :: norm_bulk_old,norm_beta_old,norm_rho_old
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_beta_sum,norm_rho_sum
  real(kind=CUSTOM_REAL) :: min_vp,max_vp,min_vs,max_vs,min_rho,max_rho
  integer :: maxindex(1)
  real(kind=CUSTOM_REAL) :: ratio_bulk,ratio_beta,ratio_rho
  integer :: iglob
  integer :: i,j,k,ispec,ier

  ! allocate arrays for storing gradient
  allocate(model_dbulk(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_dbeta(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_drho(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)

  if (ier /= 0) stop 'error allocating gradient arrays'

  ! initializes arrays
  model_dbulk = 0.0_CUSTOM_REAL
  model_dbeta = 0.0_CUSTOM_REAL
  model_drho = 0.0_CUSTOM_REAL

  ! old kernel/gradient
  ! length ( gamma_(n-1)^T * lambda_(n-1) )
  norm_bulk_old = sum( kernel_bulk_old * kernel_bulk_old )
  norm_beta_old = sum( kernel_beta_old * kernel_beta_old )
  norm_rho_old = sum( kernel_rho_old * kernel_rho_old )

  call sum_all_cr(norm_bulk_old,norm_bulk_sum)
  call sum_all_cr(norm_beta_old,norm_beta_sum)
  call sum_all_cr(norm_rho_old,norm_rho_sum)

  ! don't use square root, just take gamma^T * gamma
  norm_bulk_old = norm_bulk_sum
  norm_beta_old = norm_beta_sum
  norm_rho_old = norm_rho_sum

  if (myrank == 0) then
    print *,'norm squared old gradient:'
    print *,'  bulk : ',norm_bulk_old
    print *,'  beta: ',norm_beta_old
    print *,'  rho  : ',norm_rho_old
    print *

    ! checks lengths
    if (norm_bulk_old < 1.e-22) call exit_mpi(myrank,'norm old gradient bulk is zero')
    if (norm_beta_old < 1.e-22) call exit_mpi(myrank,'norm old gradient beta is zero')
    if (norm_rho_old < 1.e-22) call exit_mpi(myrank,'norm old gradient rho is zero')
  endif

  ! Powell, 1977: checks orthogonality between old and new gradient
  ! gets length of ( gamma_(n-1)^T * gamma_n )
  norm_bulk = sum( kernel_bulk_old * kernel_bulk )
  norm_beta = sum( kernel_beta_old * kernel_beta )
  norm_rho = sum( kernel_rho_old * kernel_rho )

  call sum_all_cr(norm_bulk,norm_bulk_sum)
  call sum_all_cr(norm_beta,norm_beta_sum)
  call sum_all_cr(norm_rho,norm_rho_sum)

  if (myrank == 0) then
    ! ratio:  ( g_n * g_n-1) / ( g_n-1 * g_n-1)
    ratio_bulk = norm_bulk_sum / norm_bulk_old
    ratio_beta = norm_beta_sum / norm_beta_old
    ratio_rho = norm_rho_sum / norm_rho_old

    ! if ratio > 0.2 (empirical threshold value), then one should restart with a steepest descent
    print *,'Powell ratio: ( > 0.2 then restart with steepest descent)'
    print *,'  bulk : ',ratio_bulk
    print *,'  beta: ',ratio_beta
    print *,'  rho  : ',ratio_rho
    print *

    if (ratio_bulk > 0.2 .and. ratio_beta > 0.2 .and. ratio_rho > 0.2) then
      print *,'  critical ratio found!'
      print *
      print *,'****************'
      print *
      print *,'  Please consider doing a steepest descent instead cg...'
      print *
      print *,'****************'
    endif
  endif


  ! difference kernel/gradient
  ! length ( ( gamma_n - gamma_(n-1))^T * lambda_n )
  norm_bulk = sum( (kernel_bulk - kernel_bulk_old) * kernel_bulk )
  norm_beta = sum( (kernel_beta - kernel_beta_old) * kernel_beta )
  norm_rho = sum( (kernel_rho - kernel_rho_old) * kernel_rho )

  call sum_all_cr(norm_bulk,norm_bulk_sum)
  call sum_all_cr(norm_beta,norm_beta_sum)
  call sum_all_cr(norm_rho,norm_rho_sum)

  ! don't take square root, since norm_bulk_sum could be negative
  ! just use (gamma_n - gamma_n-1)^T * lambda_n
  norm_bulk = norm_bulk_sum
  norm_beta = norm_beta_sum
  norm_rho = norm_rho_sum

  if (myrank == 0) then
    print *,'norm squared difference gradient:'
    print *,'  bulk : ',norm_bulk
    print *,'  beta: ',norm_beta
    print *,'  rho  : ',norm_rho
    print *
  endif

  ! calculates ratio based on Polak & Ribiere (1969)
  if (myrank == 0) then
    if (USE_SEPARATE_CG_STEPLENGTHS) then
      ! calculates steplength alpha for each parameter
      alpha_bulk = norm_bulk / norm_bulk_old
      alpha_beta = norm_beta / norm_beta_old
      alpha_rho = norm_rho / norm_rho_old

      ! only if contribution is positive it will be considered, otherwise
      ! we set it to zero so that it becomes a steepest descent update
      if (alpha_bulk < 0.0) then
        alpha_bulk = 0.0
      endif
      if (alpha_beta < 0.0) then
        alpha_beta = 0.0
      endif
      if (alpha_rho < 0.0) then
        alpha_rho = 0.0
      endif

    else
      ! calculates only a single steplength applied to all
      alpha_all = (norm_bulk + norm_beta + norm_rho) &
                  / (norm_bulk_old + norm_beta_old + norm_rho_old)
      ! only if contribution is positive it will be considered, otherwise
      ! we set it to zero so that it becomes a steepest descent update
      if (alpha_all < 0.0) then
        alpha_all = 0.0
      endif
      !!! Kai added !!!
      if (ratio_bulk > 0.2 .and. ratio_beta > 0.2 .and. ratio_rho > 0.2) then
        print *,'Powell ratio > 0.2, restart with steepest descent'
        print *,'Set alpha to zero for bulk, beta, and rho.'
        alpha_all = 0.0
      endif
      !!! 

      ! sets each steplength to same single one
      alpha_bulk = alpha_all
      alpha_beta = alpha_all
      alpha_rho = alpha_all
    endif
    ! user output
    print *,'alpha gradient:'
    print *,'  bulk : ',alpha_bulk
    print *,'  beta: ',alpha_beta
    print *,'  rho  : ',alpha_rho
    print *
  endif
  ! broadcast values from rank 0 to all others
  call bcast_all_singlecr(alpha_bulk)
  call bcast_all_singlecr(alpha_beta)
  call bcast_all_singlecr(alpha_rho)

  ! initializes kernel maximum
  depthmax(:) = 0._CUSTOM_REAL

  ! gradient in negative direction
  if (USE_OLD_GRADIENT) then
    ! uses old kernel/gradient updates ( phi_n-1)
    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

              ! note: uses old gradient update (phi_(n-1) as model_bulk_old, but
              !       given in negative gradient direction

              ! for bulk
              model_dbulk(i,j,k,ispec) = - kernel_bulk(i,j,k,ispec) &
                                         + alpha_bulk * model_dbulk_old(i,j,k,ispec)
              ! for shear
              model_dbeta(i,j,k,ispec) = - kernel_beta(i,j,k,ispec) &
                                          + alpha_beta * model_dbeta_old(i,j,k,ispec)

              ! for eta
              model_drho(i,j,k,ispec) = - kernel_rho(i,j,k,ispec) &
                                        + alpha_rho * model_drho_old(i,j,k,ispec)

              ! determines maximum kernel betav value within given radius
              if (USE_DEPTH_RANGE_MAXIMUM) then
                ! get radius of point
                iglob = ibool(i,j,k,ispec)
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
  else
    ! uses only old kernel/gradient
    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

              ! note: uses old kernels (lambda_(n-1) ) in negative gradient direction
              ! for bulk
              model_dbulk(i,j,k,ispec) = - kernel_bulk(i,j,k,ispec) &
                                         - alpha_bulk * kernel_bulk_old(i,j,k,ispec)
              ! for shear
              model_dbeta(i,j,k,ispec) = - kernel_beta(i,j,k,ispec) &
                                          - alpha_beta * kernel_beta_old(i,j,k,ispec)

              ! for eta
              model_drho(i,j,k,ispec) = - kernel_rho(i,j,k,ispec) &
                                        - alpha_rho * kernel_rho_old(i,j,k,ispec)


              ! determines maximum kernel betav value within given radius
              if (USE_DEPTH_RANGE_MAXIMUM) then
                ! get radius of point
                iglob = ibool(i,j,k,ispec)
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
  endif

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

  if (myrank == 0) then
    print *,'scaled gradient:'
    print *,'  vp min/max : ',min_vp,max_vp
    print *,'  vs min/max: ',min_vs,max_vs
    print *,'  eta min/max  : ',min_rho,max_rho
    print *
  endif
  call synchronize_all()

  end subroutine get_cg_direction_iso

  subroutine get_cg_direction_tiso()

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

  use tomography_kernels_tiso
  use tomography_kernels_tiso_cg
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha_bulk,alpha_betav,alpha_betah,alpha_eta,alpha_all
  real(kind=CUSTOM_REAL) :: alpha_alphav,alpha_alphah
  real(kind=CUSTOM_REAL) :: minmax(4),depthmax(2),depthmax_radius(2),max
  real(kind=CUSTOM_REAL) :: r,rmax_vsv,rmax_vsh,depthmax_depth
  ! gradient vector norm ( v^T * v )
  real(kind=CUSTOM_REAL) :: norm_bulk,norm_betav,norm_betah,norm_eta
  real(kind=CUSTOM_REAL) :: norm_alphav,norm_alphah
  real(kind=CUSTOM_REAL) :: norm_bulk_old,norm_betav_old,norm_betah_old,norm_eta_old
  real(kind=CUSTOM_REAL) :: norm_alphav_old,norm_alphah_old
  real(kind=CUSTOM_REAL) :: norm_bulk_sum,norm_betav_sum, &
    norm_betah_sum,norm_eta_sum
  real(kind=CUSTOM_REAL) :: norm_alphav_sum, norm_alphah_sum
  real(kind=CUSTOM_REAL) :: min_vsv,min_vsh,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk
  real(kind=CUSTOM_REAL) :: min_vpv,min_vph,max_vpv,max_vph
  integer :: maxindex(1)
  real(kind=CUSTOM_REAL) :: ratio_bulk,ratio_betav,ratio_betah,ratio_eta
  real(kind=CUSTOM_REAL) :: ratio_alphav,ratio_alphah
  integer :: iglob
  integer :: i,j,k,ispec,ier

  ! allocate arrays for storing gradient
  ! transversely isotropic arrays
  if (USE_ALPHA_BETA_RHO_TISO) then
     allocate(model_dalphav(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_dalphah(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_dbetav(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_dbetah(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_deta(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  else
     allocate(model_dbulk(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_dbetav(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_dbetah(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_deta(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  endif
  if (ier /= 0) stop 'error allocating gradient arrays'

  ! initializes arrays
  if (USE_ALPHA_BETA_RHO_TISO) then
  model_dalphav = 0.0_CUSTOM_REAL
  model_dalphah = 0.0_CUSTOM_REAL
  else
  model_dbulk = 0.0_CUSTOM_REAL
  endif
  model_dbetav = 0.0_CUSTOM_REAL
  model_dbetah = 0.0_CUSTOM_REAL
  model_deta = 0.0_CUSTOM_REAL

  ! old kernel/gradient
  ! length ( gamma_(n-1)^T * lambda_(n-1) )
  if (USE_ALPHA_BETA_RHO_TISO) then
  norm_alphav_old = sum( kernel_alphav_old * kernel_alphav_old )
  norm_alphah_old = sum( kernel_alphah_old * kernel_alphah_old )
  else
  norm_bulk_old = sum( kernel_bulk_old * kernel_bulk_old )
  endif
  norm_betav_old = sum( kernel_betav_old * kernel_betav_old )
  norm_betah_old = sum( kernel_betah_old * kernel_betah_old )
  norm_eta_old = sum( kernel_eta_old * kernel_eta_old )

  if (USE_ALPHA_BETA_RHO_TISO) then
  call sum_all_cr(norm_alphav_old,norm_alphav_sum)
  call sum_all_cr(norm_alphah_old,norm_alphah_sum)
  else
  call sum_all_cr(norm_bulk_old,norm_bulk_sum)
  endif 
  call sum_all_cr(norm_betav_old,norm_betav_sum)
  call sum_all_cr(norm_betah_old,norm_betah_sum)
  call sum_all_cr(norm_eta_old,norm_eta_sum)

  ! don't use square root, just take gamma^T * gamma
  if (USE_ALPHA_BETA_RHO_TISO) then
  norm_alphav_old = norm_alphav_sum
  norm_alphah_old = norm_alphah_sum
  else 
  norm_bulk_old = norm_bulk_sum
  endif
  norm_betav_old = norm_betav_sum
  norm_betah_old = norm_betah_sum
  norm_eta_old = norm_eta_sum

  if (myrank == 0) then
    print *,'norm squared old gradient:'
    if (USE_ALPHA_BETA_RHO_TISO) then
    print *,'  alphav: ',norm_alphav_old
    print *,'  alphah: ',norm_alphah_old
    else
    print *,'  bulk : ',norm_bulk_old
    endif
    print *,'  betav: ',norm_betav_old
    print *,'  betah: ',norm_betah_old
    print *,'  eta  : ',norm_eta_old
    print *

    ! checks lengths
    if (USE_ALPHA_BETA_RHO_TISO) then
    if (norm_alphav_old < 1.e-22) call exit_mpi(myrank,'norm old gradient alphav is zero')
    if (norm_alphah_old < 1.e-22) call exit_mpi(myrank,'norm old gradient alphah is zero')
    else  
    if (norm_bulk_old < 1.e-22) call exit_mpi(myrank,'norm old gradient bulk is zero')
    endif
    if (norm_betav_old < 1.e-22) call exit_mpi(myrank,'norm old gradient betav is zero')
    if (norm_betah_old < 1.e-22) call exit_mpi(myrank,'norm old gradient betah is zero')
    if (norm_eta_old < 1.e-22) call exit_mpi(myrank,'norm old gradient eta is zero')
  endif

  ! Powell, 1977: checks orthogonality between old and new gradient
  ! gets length of ( gamma_(n-1)^T * gamma_n )
  if (USE_ALPHA_BETA_RHO_TISO) then
  norm_alphav = sum( kernel_alphav_old * kernel_alphav )
  norm_alphah = sum( kernel_alphah_old * kernel_alphah )
  else 
  norm_bulk = sum( kernel_bulk_old * kernel_bulk )
  endif
  norm_betav = sum( kernel_betav_old * kernel_betav )
  norm_betah = sum( kernel_betah_old * kernel_betah )
  norm_eta = sum( kernel_eta_old * kernel_eta )

  if (USE_ALPHA_BETA_RHO_TISO) then
  call sum_all_cr(norm_alphav,norm_alphav_sum)
  call sum_all_cr(norm_alphah,norm_alphah_sum)
  else 
  call sum_all_cr(norm_bulk,norm_bulk_sum)
  endif 
  call sum_all_cr(norm_betav,norm_betav_sum)
  call sum_all_cr(norm_betah,norm_betah_sum)
  call sum_all_cr(norm_eta,norm_eta_sum)

  if (myrank == 0) then
    ! ratio:  ( g_n * g_n-1) / ( g_n-1 * g_n-1)
    if (USE_ALPHA_BETA_RHO_TISO) then
    ratio_alphav = norm_alphav_sum / norm_alphav_old
    ratio_alphah = norm_alphah_sum / norm_alphah_old
    else 
    ratio_bulk = norm_bulk_sum / norm_bulk_old
    endif
    ratio_betav = norm_betav_sum / norm_betav_old
    ratio_betah = norm_betah_sum / norm_betah_old
    ratio_eta = norm_eta_sum / norm_eta_old

    ! if ratio > 0.2 (empirical threshold value), then one should restart with a steepest descent
    print *,'Powell ratio: ( > 0.2 then restart with steepest descent)'
    if (USE_ALPHA_BETA_RHO_TISO) then
    print *,'  alphav: ',ratio_alphav
    print *,'  alphah: ',ratio_alphah
    else 
    print *,'  bulk : ',ratio_bulk
    endif
    print *,'  betav: ',ratio_betav
    print *,'  betah: ',ratio_betah
    print *,'  eta  : ',ratio_eta
    print *
    if (USE_ALPHA_BETA_RHO_TISO) then
    if (ratio_alphav > 0.2 .and. ratio_alphah > 0.2 .and. ratio_betav > 0.2  & 
      .and. ratio_betah > 0.2 .and. ratio_eta > 0.2) then
      print *,'  critical ratio found!'
      print *
      print *,'****************'
      print *
      print *,'  Please consider doing a steepest descent instead cg...'
      print *
      print *,'****************'
    endif
    else
    if (ratio_bulk > 0.2 .and. ratio_betav > 0.2 .and. ratio_betah > 0.2 &
      .and. ratio_eta > 0.2) then
      print *,'  critical ratio found!'
      print *
      print *,'****************'
      print *
      print *,'  Please consider doing a steepest descent instead cg...'
      print *
      print *,'****************'
    endif
    endif
  endif


  ! difference kernel/gradient
  ! length ( ( gamma_n - gamma_(n-1))^T * lambda_n )
  if (USE_ALPHA_BETA_RHO_TISO) then
  norm_alphav = sum( (kernel_alphav - kernel_alphav_old) * kernel_alphav )
  norm_alphah = sum( (kernel_alphah - kernel_alphah_old) * kernel_alphah )
  else 
  norm_bulk = sum( (kernel_bulk - kernel_bulk_old) * kernel_bulk )
  endif
  norm_betav = sum( (kernel_betav - kernel_betav_old) * kernel_betav )
  norm_betah = sum( (kernel_betah - kernel_betah_old) * kernel_betah )
  norm_eta = sum( (kernel_eta - kernel_eta_old) * kernel_eta )

  if (USE_ALPHA_BETA_RHO_TISO) then
  call sum_all_cr(norm_alphav,norm_alphav_sum)
  call sum_all_cr(norm_alphah,norm_alphah_sum)
  else 
  call sum_all_cr(norm_bulk,norm_bulk_sum)
  endif
  call sum_all_cr(norm_betav,norm_betav_sum)
  call sum_all_cr(norm_betah,norm_betah_sum)
  call sum_all_cr(norm_eta,norm_eta_sum)

  ! don't take square root, since norm_bulk_sum could be negative
  ! just use (gamma_n - gamma_n-1)^T * lambda_n
  if (USE_ALPHA_BETA_RHO_TISO) then
  norm_alphav = norm_alphav_sum
  norm_alphah = norm_alphah_sum
  else 
  norm_bulk = norm_bulk_sum
  endif
  norm_betav = norm_betav_sum
  norm_betah = norm_betah_sum
  norm_eta = norm_eta_sum

  if (myrank == 0) then
    print *,'norm squared difference gradient:'
    if (USE_ALPHA_BETA_RHO_TISO) then
    print *,'  alphav: ',norm_alphav
    print *,'  alphah: ',norm_alphah
    else 
    print *,'  bulk : ',norm_bulk
    endif
    print *,'  betav: ',norm_betav
    print *,'  betah: ',norm_betah
    print *,'  eta  : ',norm_eta
    print *
  endif

  ! calculates ratio based on Polak & Ribiere (1969)
  if (myrank == 0) then
    if (USE_SEPARATE_CG_STEPLENGTHS) then
      ! calculates steplength alpha for each parameter
      if (USE_ALPHA_BETA_RHO_TISO) then
      alpha_alphav = norm_alphav / norm_alphav_old
      alpha_alphah = norm_alphah / norm_alphah_old
      else 
      alpha_bulk = norm_bulk / norm_bulk_old
      endif
      alpha_betav = norm_betav / norm_betav_old
      alpha_betah = norm_betah / norm_betah_old
      alpha_eta = norm_eta / norm_eta_old

      ! only if contribution is positive it will be considered, otherwise
      ! we set it to zero so that it becomes a steepest descent update
      if (USE_ALPHA_BETA_RHO_TISO) then
      if (alpha_alphav < 0.0) then
        alpha_alphav = 0.0
      endif
      if (alpha_alphah < 0.0) then
        alpha_alphah = 0.0
      endif
      else 
      if (alpha_bulk < 0.0) then
        alpha_bulk = 0.0
      endif
      endif
      if (alpha_betav < 0.0) then
        alpha_betav = 0.0
      endif
      if (alpha_betah < 0.0) then
        alpha_betah = 0.0
      endif
      if (alpha_eta < 0.0) then
        alpha_eta = 0.0
      endif

    else
      ! calculates only a single steplength applied to all
      if (USE_ALPHA_BETA_RHO_TISO) then
      alpha_all = (norm_alphav + norm_alphah + norm_betav + norm_betah + norm_eta) &
                  / (norm_alphav_old+ norm_alphah_old + norm_betav_old + norm_betah_old + norm_eta_old)
      else
      alpha_all = (norm_bulk + norm_betav + norm_betah + norm_eta) &
                  / (norm_bulk_old + norm_betav_old + norm_betah_old + norm_eta_old)
      endif

      ! only if contribution is positive it will be considered, otherwise
      ! we set it to zero so that it becomes a steepest descent update
      if (alpha_all < 0.0) then
        alpha_all = 0.0
      endif

      ! sets each steplength to same single one
      if (USE_ALPHA_BETA_RHO_TISO) then
      alpha_alphav = alpha_all
      alpha_alphah = alpha_all
      else 
      alpha_bulk = alpha_all
      endif
      alpha_betav = alpha_all
      alpha_betah = alpha_all
      alpha_eta = alpha_all
    endif
    ! user output
    print *,'alpha gradient:'
    if (USE_ALPHA_BETA_RHO_TISO) then
    print *,'  alphav: ',alpha_alphav
    print *,'  alphah: ',alpha_alphah
    else 
    print *,'  bulk : ',alpha_bulk
    endif
    print *,'  betav: ',alpha_betav
    print *,'  betah: ',alpha_betah
    print *,'  eta  : ',alpha_eta
    print *
  endif
  ! broadcast values from rank 0 to all others
  if (USE_ALPHA_BETA_RHO_TISO) then
  call bcast_all_singlecr(alpha_alphav)
  call bcast_all_singlecr(alpha_alphah)
  else 
  call bcast_all_singlecr(alpha_bulk)
  endif
  call bcast_all_singlecr(alpha_betav)
  call bcast_all_singlecr(alpha_betah)
  call bcast_all_singlecr(alpha_eta)

  ! initializes kernel maximum
  depthmax(:) = 0._CUSTOM_REAL

  ! gradient in negative direction
  if (USE_OLD_GRADIENT) then
    ! uses old kernel/gradient updates ( phi_n-1)
    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

              ! note: uses old gradient update (phi_(n-1) as model_bulk_old, but
              !       given in negative gradient direction

              if (USE_ALPHA_BETA_RHO_TISO) then
              ! for P
              model_dalphav(i,j,k,ispec) = - kernel_alphav(i,j,k,ispec) &
                                          + alpha_alphav * model_dalphav_old(i,j,k,ispec)

              model_dalphah(i,j,k,ispec) = - kernel_alphah(i,j,k,ispec) &
                                          + alpha_alphah * model_dalphah_old(i,j,k,ispec)
              else 
              ! for bulk
              model_dbulk(i,j,k,ispec) = - kernel_bulk(i,j,k,ispec) &
                                         + alpha_bulk * model_dbulk_old(i,j,k,ispec)
              endif

              ! for shear
              model_dbetav(i,j,k,ispec) = - kernel_betav(i,j,k,ispec) &
                                          + alpha_betav * model_dbetav_old(i,j,k,ispec)

              model_dbetah(i,j,k,ispec) = - kernel_betah(i,j,k,ispec) &
                                          + alpha_betah * model_dbetah_old(i,j,k,ispec)

              ! for eta
              model_deta(i,j,k,ispec) = - kernel_eta(i,j,k,ispec) &
                                        + alpha_eta * model_deta_old(i,j,k,ispec)

              ! determines maximum kernel betav value within given radius
              if (USE_DEPTH_RANGE_MAXIMUM) then
                ! get radius of point
                iglob = ibool(i,j,k,ispec)
                r = z(iglob)

                ! stores maximum kernel betav/betah value in this depth slice,
                ! since betav/betah are most likely dominating
                if (r < R_TOP .and. r > R_BOTTOM) then
                  ! kernel betav value
                  max_vsv = abs( model_dbetav(i,j,k,ispec) )
                  if (depthmax(1) < max_vsv) then
                    depthmax(1) = max_vsv
                    depthmax_radius(1) = r
                  endif
                  ! kernel betav value
                  max_vsh = abs( model_dbetah(i,j,k,ispec) )
                  if (depthmax(2) < max_vsh) then
                    depthmax(2) = max_vsh
                    depthmax_radius(2) = r
                  endif
                endif
              endif

          enddo
        enddo
      enddo
    enddo
  else
    ! uses only old kernel/gradient
    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

              ! note: uses old kernels (lambda_(n-1) ) in negative gradient direction
              if (USE_ALPHA_BETA_RHO_TISO) then
              ! for shear
              model_dalphav(i,j,k,ispec) = - kernel_alphav(i,j,k,ispec) &
                                          - alpha_alphav * kernel_alphav_old(i,j,k,ispec)

              model_dalphah(i,j,k,ispec) = - kernel_alphah(i,j,k,ispec) &
                                          - alpha_alphah * kernel_alphah_old(i,j,k,ispec)
              else

              ! for bulk
              model_dbulk(i,j,k,ispec) = - kernel_bulk(i,j,k,ispec) &
                                         - alpha_bulk * kernel_bulk_old(i,j,k,ispec)
              endif

              ! for shear
              model_dbetav(i,j,k,ispec) = - kernel_betav(i,j,k,ispec) &
                                          - alpha_betav * kernel_betav_old(i,j,k,ispec)

              model_dbetah(i,j,k,ispec) = - kernel_betah(i,j,k,ispec) &
                                          - alpha_betah * kernel_betah_old(i,j,k,ispec)

              ! for eta
              model_deta(i,j,k,ispec) = - kernel_eta(i,j,k,ispec) &
                                        - alpha_eta * kernel_eta_old(i,j,k,ispec)


              ! determines maximum kernel betav value within given radius
              if (USE_DEPTH_RANGE_MAXIMUM) then
                ! get radius of point
                iglob = ibool(i,j,k,ispec)
                r = z(iglob)

                ! stores maximum kernel betav/betah value in this depth slice,
                ! since betav/betah are most likely dominating
                if (r < R_TOP .and. r > R_BOTTOM) then
                  ! kernel betav value
                  max_vsv = abs( model_dbetav(i,j,k,ispec) )
                  if (depthmax(1) < max_vsv) then
                    depthmax(1) = max_vsv
                    depthmax_radius(1) = r
                  endif
                  ! kernel betav value
                  max_vsh = abs( model_dbetah(i,j,k,ispec) )
                  if (depthmax(2) < max_vsh) then
                    depthmax(2) = max_vsh
                    depthmax_radius(2) = r
                  endif
                endif
              endif

          enddo
        enddo
      enddo
    enddo
  endif

  ! stores model_dbulk, ... arrays
  ! note: stores these new gradient before we scale them with the step length
  call write_gradient_tiso()

  ! statistics
  if (USE_ALPHA_BETA_RHO_TISO) then 
  call min_all_cr(minval(model_dalphav),min_vpv)
  call max_all_cr(maxval(model_dalphav),max_vpv)

  call min_all_cr(minval(model_dalphah),min_vph)
  call max_all_cr(maxval(model_dalphah),max_vph)
  else

  call min_all_cr(minval(model_dbulk),min_bulk)
  call max_all_cr(maxval(model_dbulk),max_bulk)
  endif

  call min_all_cr(minval(model_dbetav),min_vsv)
  call max_all_cr(maxval(model_dbetav),max_vsv)

  call min_all_cr(minval(model_dbetah),min_vsh)
  call max_all_cr(maxval(model_dbetah),max_vsh)

  call min_all_cr(minval(model_deta),min_eta)
  call max_all_cr(maxval(model_deta),max_eta)

  if (myrank == 0) then
    print *,'initial gradient updates:'
    if (USE_ALPHA_BETA_RHO_TISO) then
    print *,'  alphav min/max: ',min_vpv,max_vpv
    print *,'  alphah min/max: ',min_vph,max_vph
    else 
    print *,'  bulk min/max : ',min_bulk,max_bulk
    endif
    print *,'  betav min/max: ',min_vsv,max_vsv
    print *,'  betah min/max: ',min_vsh,max_vsh
    print *,'  eta min/max  : ',min_eta,max_eta
    print *
  endif

  ! determines maximum kernel betav value within given radius
  if (USE_DEPTH_RANGE_MAXIMUM) then
    ! maximum of all processes stored in max_vsv
    call max_all_cr(depthmax(1),max_vsv)
    call max_all_cr(depthmax(2),max_vsh)
    call max_all_cr(depthmax_radius(1),rmax_vsv)
    call max_all_cr(depthmax_radius(2),rmax_vsh)
  endif

  ! determines step length
  ! based on maximum gradient value (either vsv or vsh)
  if (myrank == 0) then

    ! determines maximum kernel betav value within given radius
    if (USE_DEPTH_RANGE_MAXIMUM) then
      depthmax(1) = max_vsv
      depthmax(2) = max_vsh
      depthmax_radius(1) = rmax_vsv
      depthmax_radius(2) = rmax_vsh

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
      minmax(1) = abs(min_vsv)
      minmax(2) = abs(max_vsv)
      minmax(3) = abs(min_vsh)
      minmax(4) = abs(max_vsh)

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
  if (USE_ALPHA_BETA_RHO_TISO) then
  norm_alphav = sum( model_dalphav * model_dalphav )
  norm_alphah = sum( model_dalphah * model_dalphah )
  else 
  norm_bulk = sum( model_dbulk * model_dbulk )
  endif
  norm_betav = sum( model_dbetav * model_dbetav )
  norm_betah = sum( model_dbetah * model_dbetah )
  norm_eta = sum( model_deta * model_deta )

  if (USE_ALPHA_BETA_RHO_TISO) then
  call sum_all_cr(norm_alphav,norm_alphav_sum)
  call sum_all_cr(norm_alphah,norm_alphah_sum)
  else 
  call sum_all_cr(norm_bulk,norm_bulk_sum)
  endif
  call sum_all_cr(norm_betav,norm_betav_sum)
  call sum_all_cr(norm_betah,norm_betah_sum)
  call sum_all_cr(norm_eta,norm_eta_sum)

  if (myrank == 0) then
    if (USE_ALPHA_BETA_RHO_TISO) then
    norm_alphav = sqrt(norm_alphav_sum)
    norm_alphah = sqrt(norm_alphah_sum)
    else 
    norm_bulk = sqrt(norm_bulk_sum)
    endif
    norm_betav = sqrt(norm_betav_sum)
    norm_betah = sqrt(norm_betah_sum)
    norm_eta = sqrt(norm_eta_sum)

    print *,'norm model updates:'
    if (USE_ALPHA_BETA_RHO_TISO) then
    print *,'  alphav: ',norm_alphav
    print *,'  alphah: ',norm_alphah
    else 
    print *,'  bulk : ',norm_bulk
    endif
    print *,'  betav: ',norm_betav
    print *,'  betah: ',norm_betah
    print *,'  eta  : ',norm_eta
    print *
  endif

  ! multiply model updates by a subjective factor that will change the step
  if (USE_ALPHA_BETA_RHO_TISO) then
  model_dalphav(:,:,:,:) = step_length * model_dalphav(:,:,:,:)
  model_dalphah(:,:,:,:) = step_length * model_dalphah(:,:,:,:)
  else 
  model_dbulk(:,:,:,:) = step_length * model_dbulk(:,:,:,:)
  endif
  model_dbetav(:,:,:,:) = step_length * model_dbetav(:,:,:,:)
  model_dbetah(:,:,:,:) = step_length * model_dbetah(:,:,:,:)
  model_deta(:,:,:,:) = step_length * model_deta(:,:,:,:)

  ! statistics
  if (USE_ALPHA_BETA_RHO_TISO) then
  call min_all_cr(minval(model_dalphav),min_vpv)
  call max_all_cr(maxval(model_dalphav),max_vpv)

  call min_all_cr(minval(model_dalphah),min_vph)
  call max_all_cr(maxval(model_dalphah),max_vph)
  else

  call min_all_cr(minval(model_dbulk),min_bulk)
  call max_all_cr(maxval(model_dbulk),max_bulk)
  endif

  call min_all_cr(minval(model_dbetav),min_vsv)
  call max_all_cr(maxval(model_dbetav),max_vsv)

  call min_all_cr(minval(model_dbetah),min_vsh)
  call max_all_cr(maxval(model_dbetah),max_vsh)

  call min_all_cr(minval(model_deta),min_eta)
  call max_all_cr(maxval(model_deta),max_eta)

  if (myrank == 0) then
    print *,'scaled gradient:'
    if (USE_ALPHA_BETA_RHO_TISO) then
    print *,'  alphav min/max: ',min_vpv,max_vpv
    print *,'  alphah min/max: ',min_vph,max_vph
    else 
    print *,'  bulk min/max : ',min_bulk,max_bulk
    endif
    print *,'  betav min/max: ',min_vsv,max_vsv
    print *,'  betah min/max: ',min_vsh,max_vsh
    print *,'  eta min/max  : ',min_eta,max_eta
    print *
  endif
  call synchronize_all()

  end subroutine get_cg_direction_tiso

