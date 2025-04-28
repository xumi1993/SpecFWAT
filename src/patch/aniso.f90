module aniso
  use config
  use fwat_constants

  implicit none

  real(kind=CUSTOM_REAL), private :: d11,d12,d13,d14,d15,d16,d22,d23,&
                                      d24,d25,d26,d33,d34,d35,d36, &
                                      d44,d45,d46,d55,d56,d66

  ! for anisotropy simulations in a halfspace model

  ! only related to body waves
  ! one-zeta term
  real(kind=CUSTOM_REAL), parameter :: FACTOR_CS1p_A = 0.2_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_CS1sv_A = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_CS1sh_N = 0._CUSTOM_REAL
  ! three-zeta term
  real(kind=CUSTOM_REAL), parameter :: FACTOR_CS3_L = 0._CUSTOM_REAL

  ! Relative to Love wave
  ! four-zeta term
  real(kind=CUSTOM_REAL), parameter :: FACTOR_N = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_E_N = 0._CUSTOM_REAL

  ! Relative to Rayleigh wave
  ! two-zeta term
  real(kind=CUSTOM_REAL), parameter :: FACTOR_A = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_C = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_F = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_H_F = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_B_A = 0._CUSTOM_REAL

  ! Relative to both Love wave and Rayleigh wave
  ! two-zeta term
  real(kind=CUSTOM_REAL), parameter :: FACTOR_L = 0._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: FACTOR_G_L = 0._CUSTOM_REAL

  type :: AnisoStruct
    real(kind=CUSTOM_REAL) :: rho, vp, vs, A, C, F, AL, AN, eta_aniso
    real(kind=CUSTOM_REAL) :: Bc, Bs, Gc, Gs, Hc, Hs, Ec, Es
    real(kind=CUSTOM_REAL) :: C1p, C1sv, C1sh, C3, S1p, S1sv, S1sh, S3
    real(kind=CUSTOM_REAL) :: c11, c12, c13, c14, c15, c16
    real(kind=CUSTOM_REAL) :: c22, c23, c24, c25, c26, c33
    real(kind=CUSTOM_REAL) :: c34, c35, c36, c44, c45, c46
    real(kind=CUSTOM_REAL) :: c55, c56, c66
    contains
    procedure :: hti2aniso, glob2cart, init_iso => initialize_iso, output => output_cijkl
    procedure, private :: to_cijkl
  end type

contains

  subroutine initialize_iso(this, vp, vs, rho)
    class(AnisoStruct), intent(inout) :: this
    real(kind=CUSTOM_REAL), intent(in) :: vp, vs, rho

    this%eta_aniso = 1.0_CUSTOM_REAL
    this%A = rho*vp**2
    this%C = rho*vp**2
    this%AL = rho*vs**2
    this%AN = rho*vs**2
    this%F = this%eta_aniso*(this%A - 2.*this%AL)

    ! zeta-dependant terms
    this%C1p = 0._CUSTOM_REAL
    this%C1sv = 0._CUSTOM_REAL
    this%C1sh = 0._CUSTOM_REAL
    this%S1p = 0._CUSTOM_REAL
    this%S1sv = 0._CUSTOM_REAL
    this%S1sh = 0._CUSTOM_REAL

    ! two-zeta dependant terms
    this%Gc = 0._CUSTOM_REAL
    this%Gs = 0._CUSTOM_REAL

    this%Bc = 0._CUSTOM_REAL
    this%Bs = 0._CUSTOM_REAL

    this%Hc = 0._CUSTOM_REAL
    this%Hs = 0._CUSTOM_REAL

    ! three-zeta dependant terms
    this%C3 = 0._CUSTOM_REAL
    this%S3 = 0._CUSTOM_REAL

    ! four-zeta dependant terms
    this%Ec = 0._CUSTOM_REAL
    this%Es = 0._CUSTOM_REAL
  end subroutine initialize_iso

  subroutine hti2aniso(this, L, Gc, Gs)
    class(AnisoStruct), intent(inout) :: this
    real(kind=CUSTOM_REAL), intent(in) :: Gc, Gs, L
    real(kind=CUSTOM_REAL) :: AL, AN, C1p, C1sh, C1sv, C3, Hs, S1p, S1sh, S1sv

    ! This subroutine sets the HTI parameters based on the input values
    this%AL = L
    this%AN = L
    this%F = this%eta_aniso*(this%A - 2.*this%AL)

    this%Gc = Gc
    this%Gs = Gs

    call this%to_cijkl()
    call this%glob2cart()
    
  end subroutine hti2aniso

  subroutine to_cijkl(this)
    class(AnisoStruct), intent(inout) :: this
    ! The mapping to the Cijkl matrix used in the code
    ! (1---East, 2---North, 3---up)
    ! The mapping from the elastic coefficients to the elastic tensor elements
    ! in the local Cartesian coordinate system (classical geographic) used in the
    ! global code (1---South, 2---East, 3---up)
    ! Always keep the following part when you modify this subroutine
    d11 = this%A + this%Ec + this%Bc
    d12 = this%A - 2.*this%AN - this%Ec
    d13 = this%F + this%Hc
    d14 = this%S3 + 2.*this%S1sh + 2.*this%S1p
    d15 = 2.*this%C1p + this%C3
    d16 = -this%Bs/2. - this%Es
    d22 = this%A + this%Ec - this%Bc
    d23 = this%F - this%Hc
    d24 = 2.*this%S1p - this%S3
    d25 = 2.*this%C1p - 2.*this%C1sh - this%C3
    d26 = -this%Bs/2. + this%Es
    d33 = this%C
    d34 = 2.*(this%S1p - this%S1sv)
    d35 = 2.*(this%C1p - this%C1sv)
    d36 = -this%Hs
    d44 = this%AL - this%Gc
    d45 = -this%Gs
    d46 = this%C1sh - this%C3
    d55 = this%AL + this%Gc
    d56 = this%S3 - this%S1sh
    d66 = this%AN - this%Ec
  end subroutine

  subroutine glob2cart(this)
    class(AnisoStruct), intent(inout) :: this
    ! The mapping to the global Cartesian coordinate system used in the code
    ! (1---East, 2---North, 3---up)
    this%c11 = d22
    this%c12 = d12
    this%c13 = d23
    this%c14 = - d25
    this%c15 = d24
    this%c16 = - d26
    this%c22 = d11
    this%c23 = d13
    this%c24 = - d15
    this%c25 = d14
    this%c26 = - d16
    this%c33 = d33
    this%c34 = - d35
    this%c35 = d34
    this%c36 = - d36
    this%c44 = d55
    this%c45 = - d45
    this%c46 = d56
    this%c55 = d44
    this%c56 = - d46
    this%c66 = d66
  end subroutine glob2cart

  subroutine output_cijkl(this, c11, c12, c13, c14, c15, c16, &
                    c22, c23, c24, c25, c26, c33, &
                    c34, c35, c36, c44, c45, c46, &
                    c55, c56, c66)
    class(AnisoStruct), intent(inout) :: this
    real(kind=CUSTOM_REAL), intent(out) :: &
      c11,c12,c13,c14,c15,c16,&
      c22,c23,c24,c25,c26,&
      c33,c34,c35,c36,&
      c44,c45,c46,&
      c55,c56,c66

    ! Output the elastic tensor components
    c11 = this%c11
    c12 = this%c12
    c13 = this%c13
    c14 = this%c14
    c15 = this%c15
    c16 = this%c16
    c22 = this%c22
    c23 = this%c23
    c24 = this%c24
    c25 = this%c25
    c26 = this%c26
    c33 = this%c33
    c34 = this%c34
    c35 = this%c35
    c36 = this%c36
    c44 = this%c44
    c45 = this%c45
    c46 = this%c46
    c55 = this%c55
    c56 = this%c56
    c66 = this%c66

  end subroutine output_cijkl

end module aniso