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

module tomography_par

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN, &
    NGLLX,NGLLY,NGLLZ,IIN,IOUT,FOUR_THIRDS,R_EARTH_KM,GAUSSALPHA,GAUSSBETA
  use specfem_par, only: myrank,sizeprocs,ibool

  implicit none

  ! tomography parameter settings
  include "constants_tomography.h"

  ! mesh size
  integer :: NSPEC, NGLOB

  ! volume
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: x, y, z
  !integer, dimension(:,:,:,:),allocatable :: ibool

  ! model update length
  real(kind=CUSTOM_REAL) :: step_fac,step_length

  ! MPI process
  !integer :: myrank,sizeprocs

  !!! Kai added for lbfgs
  integer,parameter:: NKERNEL=3
  !!! Mijian move m_store to tomo_par
  integer,parameter:: m_store=5 ! stored model step 3 <= m_store <= 7

  integer:: iker
  character(len=MAX_STRING_LEN) :: filename,dirname

  character(len=MAX_STRING_LEN) :: nkernel_name(NKERNEL)
  character(len=MAX_STRING_LEN) :: nmodel_name(NKERNEL)

  integer:: iter_start,iter_current

end module tomography_par

!
!-------------------------------------------------------------------------------------------------
!

module tomography_kernels_iso

  use tomography_par

  implicit none

  ! kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_bulk,kernel_beta,kernel_rho,total_kernel

  ! gradients for model updates
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_dbulk,model_dbeta,model_drho

end module tomography_kernels_iso


!
!-------------------------------------------------------------------------------------------------
!

module tomography_kernels_tiso

  use tomography_par

  implicit none

  ! kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_bulk,kernel_betav,kernel_betah,kernel_eta
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_alphav,kernel_alphah

  ! gradients for model updates
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_dbulk,model_dbetah,model_dbetav,model_deta
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_dalphah,model_dalphav

end module tomography_kernels_tiso

!
!-------------------------------------------------------------------------------------------------
!

module tomography_kernels_iso_cg

  use tomography_par

  implicit none

  ! flags to determine wheter old gradients (model_dbulk,..) can be used or
  ! if update is based on old kernels only (kernel_bulk,..)
  logical :: USE_OLD_GRADIENT

  ! kernels from former iteration
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
        kernel_alpha_old,kernel_bulk_old,kernel_beta_old,kernel_rho_old

  ! gradients for model updates from former iteration
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
        model_dbulk_old,model_dbeta_old,model_drho_old

end module tomography_kernels_iso_cg


!
!-------------------------------------------------------------------------------------------------
!

module tomography_kernels_tiso_cg

  use tomography_par

  implicit none

  ! flags to determine wheter old gradients (model_dbulk,..) can be used or
  ! if update is based on old kernels only (kernel_bulk,..)
  logical :: USE_OLD_GRADIENT

  ! kernels from former iteration
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
        kernel_bulk_old,kernel_betav_old,kernel_betah_old,kernel_eta_old &
        ,kernel_alphav_old,kernel_alphah_old

  ! gradients for model updates from former iteration
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
        model_dbulk_old,model_dbetah_old,model_dbetav_old,model_deta_old &
        ,model_dalphah_old,model_dalphav_old

end module tomography_kernels_tiso_cg

!
!-------------------------------------------------------------------------------------------------
!

module tomography_model_iso

  use tomography_par

  implicit none

  ! isotropic model files
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_vp,model_vs,model_rho
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_vp_new,model_vs_new,model_rho_new

end module tomography_model_iso

!
!-------------------------------------------------------------------------------------------------
!

module tomography_model_tiso

  use tomography_par

  implicit none

  ! transverse isotropic model files
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
        model_vpv,model_vph,model_vsv,model_vsh,model_eta,model_rho
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
        model_vpv_new,model_vph_new,model_vsv_new,model_vsh_new,model_eta_new,model_rho_new

end module tomography_model_tiso


