!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

!--------------------------------------------------------------------------------------------------
!
! generic model file
!
! note: the idea is to super-impose velocity model values on the GLL points,
!          additional to the ones assigned on the CUBIT mesh
!
! most of the routines here are place-holders, please add/implement your own routines
!
!--------------------------------------------------------------------------------------------------


  module external_model
  use constants
  use config
!---
!
! ADD YOUR MODEL HERE
!
!---
  implicit none
  integer, parameter :: dummy_size = 1

  ! only here to illustrate an example
  type model_external_variables
      sequence
      real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: vp, vs, rho, Gcp, Gsp
      real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: qkappa_atten, qmu_atten
      real(kind=CUSTOM_REAL), dimension(:), allocatable :: x, y, z
      real(kind=CUSTOM_REAL) :: dx, dy, dz
      integer nx, ny, nz
    end type model_external_variables
  type (model_external_variables) MEXT_V, ext_grid

  end module external_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_broadcast()

! standard routine to setup model

  use external_model

  use constants

  ! use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH

  implicit none

  ! local parameters
  integer :: idummy

  ! dummy to ignore compiler warnings
  idummy = myrank

  ! safety check
  ! if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) then
  !   print *,'Coupling with injection technique or mesh a chunk of the earth requires in Par_file: MODEL = coupled'
  !   stop 'Error model external'
  ! endif


!---
!
! ADD YOUR MODEL HERE
!
!---

  ! the variables read are declared and stored in structure MEXT_V

  ! initialize the model
  if(allocated(MEXT_V%vp)) deallocate(MEXT_V%vp)
  if(allocated(MEXT_V%vs)) deallocate(MEXT_V%vs)
  if(allocated(MEXT_V%rho)) deallocate(MEXT_V%rho)
  if(allocated(MEXT_V%Gcp)) deallocate(MEXT_V%Gcp)
  if(allocated(MEXT_V%Gsp)) deallocate(MEXT_V%Gsp)
  if(allocated(MEXT_V%qkappa_atten)) deallocate(MEXT_V%qkappa_atten)
  if(allocated(MEXT_V%qmu_atten)) deallocate(MEXT_V%qmu_atten)
  if(allocated(MEXT_V%x)) deallocate(MEXT_V%x)
  if(allocated(MEXT_V%y)) deallocate(MEXT_V%y)
  if(allocated(MEXT_V%z)) deallocate(MEXT_V%z)

  if (myrank == 0) call read_external_model()

  ! broadcast the information read on the main to the nodes
  call bcast_all_singlei(MEXT_V%nx)
  call bcast_all_singlei(MEXT_V%ny)
  call bcast_all_singlei(MEXT_V%nz)
  if (myrank /= 0) then
    allocate(MEXT_V%vp(MEXT_V%nx,MEXT_V%ny,MEXT_V%nz))
    allocate(MEXT_V%vs(MEXT_V%nx,MEXT_V%ny,MEXT_V%nz))
    allocate(MEXT_V%rho(MEXT_V%nx,MEXT_V%ny,MEXT_V%nz))
    allocate(MEXT_V%x(MEXT_V%nx))
    allocate(MEXT_V%y(MEXT_V%ny))
    allocate(MEXT_V%z(MEXT_V%nz))
    if (parameter_type == 2) then
      allocate(MEXT_V%Gcp(MEXT_V%nx,MEXT_V%ny,MEXT_V%nz))
      allocate(MEXT_V%Gsp(MEXT_V%nx,MEXT_V%ny,MEXT_V%nz))
    endif 
  endif
  call bcast_all_cr(MEXT_V%vp, size(MEXT_V%vp))
  call bcast_all_cr(MEXT_V%vs, size(MEXT_V%vs))
  call bcast_all_cr(MEXT_V%rho, size(MEXT_V%rho))    
  if (parameter_type == 2) then
    call bcast_all_cr(MEXT_V%Gcp, size(MEXT_V%Gcp))
    call bcast_all_cr(MEXT_V%Gsp, size(MEXT_V%Gsp))
  endif
  call bcast_all_cr(MEXT_V%x, size(MEXT_V%x))
  call bcast_all_cr(MEXT_V%y, size(MEXT_V%y))
  call bcast_all_cr(MEXT_V%z, size(MEXT_V%z))


  end subroutine model_external_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_external_model()

  use external_model
  use hdf5_interface
  use shared_input_parameters, only: TOMOGRAPHY_PATH
  use utils, only: transpose_3

  implicit none

  character(len=MAX_STRING_LEN) :: fname

!---
!
! ADD YOUR MODEL HERE
!
!---
  fname = trim(TOMOGRAPHY_PATH)//'/tomography_model.h5'
  call h5read(fname, '/x', MEXT_V%x)
  call h5read(fname, '/y', MEXT_V%y)
  call h5read(fname, '/z', MEXT_V%z)
  MEXT_V%nx = size(MEXT_V%x)
  MEXT_V%ny = size(MEXT_V%y)
  MEXT_V%nz = size(MEXT_V%z)
  call h5read(fname, '/vp', MEXT_V%vp)
  call h5read(fname, '/vs', MEXT_V%vs)
  call h5read(fname, '/rho', MEXT_V%rho)
  MEXT_V%vp = transpose_3(MEXT_V%vp)
  MEXT_V%vs = transpose_3(MEXT_V%vs)
  MEXT_V%rho = transpose_3(MEXT_V%rho)    
  if (parameter_type == 2) then
    call h5read(fname, '/gcp', MEXT_V%Gcp)
    call h5read(fname, '/gsp', MEXT_V%Gsp)
    MEXT_V%Gcp = transpose_3(MEXT_V%Gcp)
    MEXT_V%Gsp = transpose_3(MEXT_V%Gsp)
  endif

  end subroutine read_external_model


!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_external_values(xmesh,ymesh,zmesh,ispec,rho,vp,vs,qkappa_atten,qmu_atten,iflag_aniso,idomain_id )

! given a GLL point, returns super-imposed velocity model values

  !use constants, only: MIDX,MIDY,MIDZ

  use external_model
  use utils, only: interp3, interp3_nearest_simple
  use generate_databases_par, only: HUGEVAL,TINYVAL,IDOMAIN_ELASTIC,CUSTOM_REAL

  ! use create_regions_mesh_ext_par, only: xstore_unique,ystore_unique,zstore_unique, &
    ! num_free_surface_faces,free_surface_ijk,free_surface_ispec,nglob_unique

  implicit none

  ! GLL point
  double precision, intent(in) :: xmesh,ymesh,zmesh
  real(kind=CUSTOM_REAL) :: xsem, ysem, zsem
  integer, intent(in) :: ispec

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL), intent(out) :: vp,vs,rho

  ! attenuation flag
  real(kind=CUSTOM_REAL), intent(out) :: qkappa_atten,qmu_atten

  ! anisotropy flag
  integer, intent(out) :: iflag_aniso

  ! acoustic/elastic/.. domain flag ( 1 = acoustic / 2 = elastic / ... )
  integer, intent(out) :: idomain_id

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! GLL point location converted to real
  xsem = xmesh
  ysem = ymesh
  zsem = zmesh

  if (xsem < MEXT_V%x(1) .or. xsem > MEXT_V%x(MEXT_V%nx) .or. &
      ysem < MEXT_V%y(1) .or. ysem > MEXT_V%y(MEXT_V%ny) .or. &
      zsem < MEXT_V%z(1) .or. zsem > MEXT_V%z(MEXT_V%nz)) then
    rho = interp3_nearest_simple(MEXT_V%x, MEXT_V%y, MEXT_V%z, MEXT_V%rho, xsem, ysem, zsem)
    vp = interp3_nearest_simple(MEXT_V%x, MEXT_V%y, MEXT_V%z, MEXT_V%vp, xsem, ysem, zsem)
    vs = interp3_nearest_simple(MEXT_V%x, MEXT_V%y, MEXT_V%z, MEXT_V%vs, xsem, ysem, zsem)
  else
    rho = interp3(MEXT_V%x, MEXT_V%y, MEXT_V%z, MEXT_V%rho, xsem, ysem, zsem)
    vp = interp3(MEXT_V%x, MEXT_V%y, MEXT_V%z, MEXT_V%vp, xsem, ysem, zsem)
    vs = interp3(MEXT_V%x, MEXT_V%y, MEXT_V%z, MEXT_V%vs, xsem, ysem, zsem)
  endif

  ! attenuation: PREM crust value
  qmu_atten = 600._CUSTOM_REAL

  ! no Q_Kappa in this model, use a dummy very high value of 9999. as a flag for no QKappa attenuation
  qkappa_atten = 9999._CUSTOM_REAL

  ! no anisotropy
  iflag_aniso = 0

  ! elastic material
  idomain_id = IDOMAIN_ELASTIC

  end subroutine model_external_values

  subroutine model_external_values_aniso(xmesh,ymesh,zmesh,rho,vp,vs,&
                                          c11,c12,c13,c14,c15,c16, &
                                          c22,c23,c24,c25,c26,c33, &
                                          c34,c35,c36,c44,c45,c46,c55,c56,c66,&
                                          iflag_aniso )
    use external_model
    use utils, only: interp3, interp3_nearest_simple
    use aniso, only: AnisoStruct
    use generate_databases_par, only: IDOMAIN_ELASTIC,CUSTOM_REAL
    double precision, intent(in) :: xmesh,ymesh,zmesh
    real(kind=CUSTOM_REAL) :: xsem, ysem, zsem
    ! density, Vp and Vs
    real(kind=CUSTOM_REAL), intent(in) :: vp,vs,rho
    real(kind=CUSTOM_REAL) :: Gcp, Gsp

    ! anisotropy parameters
    real(kind=CUSTOM_REAL), intent(out) :: c11,c12,c13,c14,c15,c16, &
                                          c22,c23,c24,c25,c26,c33, &
                                          c34,c35,c36,c44,c45,c46,c55,c56,c66
    ! anisotropy flag
    integer, intent(out) :: iflag_aniso

    type(AnisoStruct) :: anistruct

    xsem = xmesh
    ysem = ymesh
    zsem = zmesh

    ! interpolate the anisotropic parameters
    if (xsem < MEXT_V%x(1) .or. xsem > MEXT_V%x(MEXT_V%nx) .or. &
      ysem < MEXT_V%y(1) .or. ysem > MEXT_V%y(MEXT_V%ny) .or. &
      zsem < MEXT_V%z(1) .or. zsem > MEXT_V%z(MEXT_V%nz)) then
      Gcp = interp3_nearest_simple(MEXT_V%x, MEXT_V%y, MEXT_V%z, MEXT_V%Gcp, xsem, ysem, zsem)
      Gsp = interp3_nearest_simple(MEXT_V%x, MEXT_V%y, MEXT_V%z, MEXT_V%Gsp, xsem, ysem, zsem)
    else
      Gcp = interp3(MEXT_V%x, MEXT_V%y, MEXT_V%z, MEXT_V%Gcp, xsem, ysem, zsem)
      Gsp = interp3(MEXT_V%x, MEXT_V%y, MEXT_V%z, MEXT_V%Gsp, xsem, ysem, zsem)
    endif

    call anistruct%init_iso(vp, vs, rho)
    if (parameter_type == 2) then
      call anistruct%hti2aniso(Gcp, Gsp)
    else
      call exit_mpi(0, 'Error: anisotropic model not available')
    endif
    call anistruct%output(c11,c12,c13,c14,c15,c16, &
                          c22,c23,c24,c25,c26,c33, &
                          c34,c35,c36,c44,c45,c46,c55,c56,c66)

    iflag_aniso = 1

  end subroutine model_external_values_aniso



