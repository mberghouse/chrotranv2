module characteristic_curves_WIPP_invariant_module
#include "petsc/finclude/petscsys.h"

use petscsys ! Necessary for PETSC_TRUE / PETSC_FALSE
use option_module ! Necessary for characterisic curve base class arguments
use characteristic_curves_base_module ! Needed to define base type
implicit none

#define KRP_VG_range 1,8
#define KRP_BC_range 2,3,4,12

#define KRP_SE1_range 2,8,9,11
#define KRP_SE2_range 1,3,4,5

! **************************************************************************** !
! Common WIPP Saturation Function type
! **************************************************************************** !

type, public, extends(sat_func_base_type) :: sf_WIPP_type
  private
    procedure(set_k_type), public, pointer :: setK
    procedure(set_Swj_type), pointer :: setSwj

    procedure(calc_Pc_type), pointer :: KPCPc
    procedure(calc_Sw_type), pointer :: KPCSw

    procedure(calc_Pc_type), pointer :: KRPPc
    procedure(calc_Sw_type), pointer :: KRPSw

!   PFLOTRAN object parameter             BRAGFLO equivalent
!   PetscInt  :: KRP                      ! KRP
!   PetscInt  :: KPC                      ! KPC
!   function pointers replace these. No need for enumeration

!   Effective Saturation Parameters
    PetscReal :: Swr                      ! SWR or SOCZRO
    PetscReal :: Sgr                      ! SGR
    PetscReal :: Sgr_comp                 ! 1 - SGR
    PetscReal :: Se_span                  ! 1/SETERM or 1/SETERM2
    PetscReal :: dSe_dSw                  ! SETERM or SETERM2
    PetscReal :: Semin                    ! SOCEFFMIN for KRP12

!   Brooks-Corey Parameters
    PetscReal :: lambda                   ! 1/XLAM1
    PetscReal :: lambda_nrec              ! -XLAM1

!   Van Genuchten Parameters
    PetscReal :: m                        ! XLAM4
    PetscReal :: m_nrec                   ! XLAM6
    PetscReal :: m_comp                   ! XLAM7
    PetscReal :: n                        ! 1/XLAM7
    PetscReal :: n_rec                    ! -XLAM7
    PetscReal :: Pcm_Pct                  ! Precalculated ratio for KRP1

!   Unsaturated Extension Parameters
!   PetscReal :: Pcmax                    ! PCFIX          Defined in base
    PetscReal :: Swj, Pcj
    PetscReal :: dPcj_dSwj

!   Threshold Pressure Parameters
    PetscReal :: pct
    PetscReal :: pct_a
    PetscReal :: pct_exp
    PetscReal :: permeability

!   Coefficient for analytical derivatives
    PetscReal :: k_dSw_dSe                ! -m*n*Se_span or -lambda*Se_span
  contains
! Overridden function pointers from the PFLOTRAN base class
    procedure, public :: CapillaryPressure => SFWIPPCapillaryPressure
    procedure, public :: Saturation        => SFWIPPSaturation
end type

! **************************************************************************** !
! Function prototypes for function pointers
! **************************************************************************** !

abstract interface
  function set_Swj_type(this, Swj) result (error)
    import sf_WIPP_type
    class(sf_WIPP_type), intent(inout) :: this
    PetscReal, intent(in)  :: Swj
    PetscInt :: error
  end function
  pure subroutine calc_Pc_type(this, Sw, Pc, dPc_dSw)
    import sf_WIPP_type
    class(sf_WIPP_type), intent(in) :: this
    PetscReal, intent(in)  :: Sw
    PetscReal, intent(out) :: Pc, dPc_dSw
  end subroutine
  pure subroutine calc_Sw_type(this, Pc, Sw)
    import sf_WIPP_type
    class(sf_WIPP_type), intent(in) :: this
    PetscReal, intent(in)  :: Pc
    PetscReal, intent(out) :: Sw
  end subroutine
  subroutine set_k_type(this, k)
    import sf_WIPP_type
    class(sf_WIPP_type), intent(inout) :: this
    PetscReal, intent(in)  :: k
  end subroutine
end interface

! **************************************************************************** !
! Public/Private Procedure Declarations
! **************************************************************************** !

! WIPP constructor
public  :: SFWIPPctor

! Implemented WIPP KPC procedures
private :: SFWIPPKPC1Pc , &
           SFWIPPKPC1Sw , &
           SFWIPPKPC1Swj, &
           SFWIPPKPC2Pc , &
           SFWIPPKPC2Sw , &
           SFWIPPKPC2Swj, &
           SFWIPPKPC6Pc , &
           SFWIPPKPC6Sw , &
           SFWIPPKPC6Swj

! Implemented WIPP KRP Pc procedures
private :: SFWIPPKRP1Pc , &
           SFWIPPKRP1Sw , &
           SFWIPPKRP2Pc , &
           SFWIPPKRP2Sw , &
           SFWIPPKRP3Pc , &
           SFWIPPKRP3Sw , &
           SFWIPPKRP4Pc , &
           SFWIPPKRP4Sw , &
           SFWIPPKRP5Pc , &
           SFWIPPKRP5Sw , &
           SFWIPPKrp8Pc , &
           SFWIPPKRP8Sw , &
           SFWIPPKRP9Pc , &
           SFWIPPKRP9Sw , &
           SFWIPPKRP11Pc , &
           SFWIPPKRP11Sw , &
           SFWIPPKRP12Pc , &
           SFWIPPKRP12Sw

! Implemented WIPP Permeability derived Pct procedures
private :: SFWIPPSetK , &
           SFWIPPIgnoreK

contains 

! **************************************************************************** !
! WIPP Constructor
! **************************************************************************** !

function SFWIPPctor(KRP, KPC, Swr, Sgr, expon, Pct_ignore, Pct_alpha, &
                    Pct_expon, Pcmax, Swj, Semin) result (new)
  class(sf_WIPP_type), pointer :: new
  PetscInt, intent(inout)  :: KRP, KPC
  PetscReal, intent(in) :: Swr, Sgr, expon, Pct_alpha, Pct_expon, Swj, Pcmax
  PetscReal, intent(in) :: Semin
  PetscBool, intent(in) :: Pct_ignore
  PetscInt :: error

  ! Memory allocation
  allocate(new)
  if (.not. associated(new)) return ! Memory allocation failed, abort

  ! Data validation
                                        error = 0

! Derivatives are necessary to make smooth unsaturated extensions
  new%analytical_derivative_available = .TRUE.

!  For KRP 11, the cavity model, capillary pressure is always zero.
!  For all others except KRP 9, if PCT_A is zero, the same occurs.

  if (KRP == 11 .OR. (Pct_alpha == 0d0 .AND. KRP /= 9)) then
    new%setK => SFWIPPIgnoreK
    new%setSwj => SFWIPPKPC1Swj

    new%KRPPc  => SFWIPPKRP11Pc
    new%KRPSw  => SFWIPPKRP11Sw

    new%KPCPc  => new%KRPPc
    new%KPCSw  => new%KRPSW
    return
  end if

  ! Otherwise, if KRP is in branch table, set function pointers, else flag error
  select case(KRP)
  case (1)
    new%KRPPc => SFWIPPKRP1Pc
    new%KRPSw => SFWIPPKRP1Sw
  case (2)
    new%KRPPc => SFWIPPKRP2Pc
    new%KRPSw => SFWIPPKRP2Sw
  case (3)
    new%KRPPc => SFWIPPKRP3Pc
    new%KRPSw => SFWIPPKRP3Sw
  case (4)
    new%KRPPc => SFWIPPKRP4Pc
    new%KRPSw => SFWIPPKRP4Sw
  case (5)
    new%KRPPc => SFWIPPKRP5Pc
    new%KRPSw => SFWIPPKRP5Sw
  case (8)
    new%KRPPc => SFWIPPKRP8Pc
    new%KRPSw => SFWIPPKRP8Sw
  case (9)
    new%KRPPc => SFWIPPKRP9Pc
    new%KRPSw => SFWIPPKRP9Sw
  case (11) ! Cavity model, caught above
  case (12)
    new%KRPPc => SFWIPPKRP12Pc
    new%KRPSw => SFWIPPKRP12Sw
  case default
                                        error = error + 1
  end select

  ! If KPC is in branch table, set function pointer, else flag error
  select case(KPC)
  case (1) ! 0 at or below residual
    new%KPCPc  => SFWIPPKPC1Pc 
    new%KPCSw  => SFWIPPKPC1Sw
    new%setSwj => SFWIPPKPC1Swj
  case (2) ! Pcmax at or below residual
    new%KPCPc  => SFWIPPKPC2Pc
    new%KPCSw  => SFWIPPKPC2Sw
    new%setSwj => SFWIPPKPC2Swj
  case (6) ! Linear at or below junction
    new%KPCPc  => SFWIPPKPC6Pc
    new%KPCSw  => SFWIPPKPC6Sw
    new%setSwj => SFWIPPKPC6Swj
  case default
                                        error = error + 2
  end select

  ! Check the residual saturations, and their sum, are bound between 0 and 1
  if (Swr < 0d0 .or. Swr >= 1d0)        error = error + 4
  if (Sgr < 0d0 .or. Sgr >= 1d0)        error = error + 8
  if (Swr + Sgr > 1d0)                  error = error + 12

  ! Check the exponent parameter is valid for the KRP option chosen
  select case(KRP)
  case (KRP_BC_range) ! Brooks-Corey types
    if (expon <= 0d0)                   error = error + 16
  case (KRP_VG_range) ! Van Genuchten types
    if (expon <= 0d0 .or. expon >= 1d0) error = error + 16
  case default
    ! Other KRP functions do not use lambda or m
  end select

  ! Abort if errors caught in data validation
  if (error /= 0) then
    ! Potentally write which parameters were invalid to error stream
    deallocate(new)
    nullify(new)
    return
  end if

  ! Assign residual saturation parameters
  new%Swr = Swr
  new%Sgr = Sgr
  select case(KRP)
  case (KRP_SE1_range)
    new%Sgr_comp = 1d0
  case (KRP_SE2_range)
    new%Sgr_comp = 1d0 - new%Sgr
  case (12)
    new%Sgr_comp = 1d0
    new%Semin = Semin      ! Used for KRP12 only
    new%Swr = Swr - Semin  ! Swr is BRAGFLO SOCZRO = SOCMIN - SOCEFFMIN
  ! Note, BRAGFLO UM 6.02 indicates it should be +, but the code is -
  ! I.e. 1 - (Smin - Seffmin) /= 1 - Smin - Seffmin
  case default
  end select
  
  new%Se_span = new%Sgr_comp - new%Swr 
  new%dSe_dSw = 1d0 / new%Se_span

  ! Depending on KRP, expon is either VG m or BC lambda
  select case(KRP)
  case (KRP_BC_range) ! Brooks-Corey types
    new%lambda = expon
    new%lambda_nrec = -1d0/new%lambda
  case (KRP_VG_range) ! Van Genuchten types
! In BRAGFLO, VG functions use the closed form Mualem condition n = 1/(1-m)
    new%m = expon
    new%m_nrec = -1d0 / new%m
    new%m_comp =  1d0 - new%m
    new%n      =  1d0 / new%m_comp
    new%n_rec  =  1d0 / new%n

! KRP1 relies on an approximate conversion of BC to VG
    new%lambda = new%m/(1d0-new%m)
    new%lambda_nrec = -1d0/new%lambda
    new%pcm_pct = 0.5d0**new%lambda_nrec
  case default        ! Unused parameters
    new%lambda = 0d0
    new%m = 0d0
  end select

  ! Depending on pct_ignore, pct_alpha is either alpha or pct_a
  new%permeability = -1d0   ! Initialize
  if (pct_ignore) then
    new%setK  => SFWIPPIgnoreK
    new%pct_a   = 0d0       ! Not used
    new%pct_exp = 0d0       ! Not used
    new%pct     = pct_alpha ! Pct permanently set to alpha
  else
    new%setK  => SFWIPPSetK
    new%pct_a   = pct_alpha
    new%pct_exp = pct_expon
    new%pct     = 1E6       ! Initialize to 1 MPa to permit test function. 
  end if

  ! Initialize unsaturated extensions 
  error = new%setSwj(Swj)
  ! Abort if error in unsaturated extension (e.g. Swj < Swr)
  if (error /= 0) then
    ! Potentally write which parameters were invalid to error stream
    deallocate(new)
    nullify(new)
    return
  end if

end function

! **************************************************************************** !
! WIPP Saturation Function Roots
! **************************************************************************** !

subroutine SFWIPPCapillaryPressure(this, liquid_saturation, capillary_pressure,&
                                         dpc_dsatl, option)
  class(sf_WIPP_type)              :: this
  PetscReal, intent(in)            :: liquid_saturation
  PetscReal, intent(out)           :: capillary_pressure, dpc_dsatl
  type(option_type), intent(inout) :: option

  if (liquid_saturation <= this%Swj) then
    call this%KPCPc(liquid_saturation, capillary_pressure, dpc_dsatl)
  else
    call this%KRPPc(liquid_saturation, capillary_pressure, dpc_dsatl)
  end if
end subroutine

! **************************************************************************** !

subroutine SFWIPPSaturation(this, capillary_pressure, liquid_saturation, &
                                  dsat_dpres, option)
  class(sf_WIPP_type)              :: this
  PetscReal, intent(in)            :: capillary_pressure
  PetscReal, intent(out)           :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  if (capillary_pressure >= this%Pcj) then
    call this%KPCSw(capillary_pressure, liquid_saturation)
  else
    call this%KRPSw(capillary_pressure, liquid_saturation)
  end if
  dsat_dpres = 0d0 ! analytic derivatives not yet available for Richard's mode
end subroutine

! **************************************************************************** !
! WIPP PCT Subroutines
! **************************************************************************** !

subroutine SFWIPPSetK(this, permeability)
  class(sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: permeability
  PetscInt :: error
 
  ! Permeability changes with material region, but because this occurs in 
  ! blocks, frequently permeability has not changed from call to call.
  ! Because exponentiaton is expensive, avoid this if possible.
  if (permeability /= this%permeability) then
    this%permeability = permeability
    this%pct = this%pct_a * permeability ** this%pct_exp
    ! Update unsaturated extension to reflect new Pct
    error = this%setSwj(this%Swj)
  end if
end subroutine

! **************************************************************************** !

subroutine SFWIPPIgnoreK(this, permeability)
  class(sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: permeability
end subroutine

! **************************************************************************** !
! WIPP KPC Subroutines
! **************************************************************************** !

function SFWIPPKPC1Swj(this,Swj) result (error)
  class (sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Swj
  PetscInt :: error

  ! Ignore input, set Swj to Swr
  this%Swj = this%Swr
  this%Pcj = huge(this%Pcj)
  error = 0
end function

! **************************************************************************** !

pure subroutine SFWIPPKPC1Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  Pc = 0d0
  dPc_dSw = 0d0
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKPC1Sw(this, Pc, Sw)
 class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw

  Sw = this%Swr
end subroutine

! **************************************************************************** !

function SFWIPPKPC2Swj(this,Swj) result (error)
  class (sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Swj
  PetscInt :: error

  ! Ignore input, calculate Swj based on Pcmax
  call this%KRPSw(this%Pcmax, this%Swj)
  this%Pcj = this%Pcmax

  error = 0
end function

! **************************************************************************** !

pure subroutine SFWIPPKPC2Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  Pc = this%Pcmax
  dPc_dSw = 0d0
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKPC2Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw

  Sw = this%Swj
end subroutine

! **************************************************************************** !

function SFWIPPKPC6Swj(this,Swj) result (error)
  class(sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Swj
  PetscInt :: error

  if (Swj > this%Swr) then ! Linearly extrapolate from the valid Swj
    this%Swj = Swj
    call this%KRPPc(Swj, this%Pcj, this%dPcj_dSwj)
    this%Pcmax = this%Pcj - this%dPcj_dSwj * Swj
    error = 0
  else ! Invalid Swj
    error = 1
  end if
end function

! **************************************************************************** !

pure subroutine SFWIPPKPC6Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  if (Sw > 0d0) then                                  ! Linear interpolation
    Pc = this%Pcmax + this%dPcj_dSwj*Sw
  else                                                ! y-intercept
    Pc = this%Pcmax
  end if
  dPc_dSw = this%dPcj_dSwj
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKPC6Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw

  if (Pc < this%Pcmax) then                           ! Linear interpolation
    Sw = (Pc - this%Pcmax) / this%dPcj_dSwj
  else                                                ! y-intercept
    Sw = 0d0
  end if
end subroutine

! **************************************************************************** !
! BRAGFLO KRP Subroutines
! **************************************************************************** !

pure subroutine SFWIPPKRP1Pc(this, Sw, Pc, dPc_dSw)
! van Genuchten w/ Sgr > 0
! Author: Heeho Park; Modified by Jennifer Frederick
! Date: 11/17/16; Modified 04/26/2017 ; Refactored 10/12/2021
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se, Se_mrtrec, aPc_n

  if (Sw >= this%Sgr_comp) then ! Saturated limit
    Pc = 0d0
    dPc_dSw = -huge(dPc_dSw) ! TODO calculate finite difference limit
  else
    Se = (Sw-this%Swr) * this%dSe_dSw
    Se_mrtrec = Se**this%m_nrec
    aPc_n = Se_mrtrec - 1d0
    Pc = this%Pct*this%Pcm_Pct * aPc_n**this%m_comp
    dPc_dSw = (this%dSe_dSw*Pc/this%m/this%n) / (Se/Se_mrtrec-Se)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP1Sw(this, Pc, Sw)
!------- Modified van Genuchten/Parker model (SGR>=0)
! Author: Heeho Park; Modified by Jennifer Frederick
! Date: 11/17/16; Modified 04/26/2017
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se

  if (Pc <= 0d0) then
    Sw = 1d0
  else
    Se = ( Pc**this%n / (this%Pct * this%Pcm_Pct) + 1d0)**(-this%m)
    Sw = this%Swr + this%Se_span * Se
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP2Pc(this, Sw, Pc, dPc_dSw)
! Brooks-Corey w/ Sgr = 0
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se

  Se      = (Sw - this%Swr) * this%dSe_dSw
  Pc      = this%Pct * Se**this%lambda_nrec
  dPc_dSw = (-this%dSe_dSw/this%lambda) * Pc / Se
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP2Sw(this, Pc, Sw)
! Brooks-Corey w/ Sgr = 0
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se

  if (Pc <= this%Pct) then
    Sw = 1d0
  else
    Se = (this%Pct/Pc)**this%lambda
    Sw = this%Swr + this%Se_span*Se
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP3Pc(this, Sw, Pc, dPc_dSw)
! Brooks-Corey w/ Sgr > 0
! Flat above Sgr
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se

  if (Sw >= this%Sgr_comp) then
    Pc = this%Pct
    dPc_dSw = 0d0
  else
    Se = (Sw - this%Swr) * this%dSe_dSw
    Pc = this%Pct * Se**this%lambda_nrec
    dPc_dSw = -this%dSe_dSw*Pc / (this%lambda*Se)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP3Sw(this, Pc, Sw)
! Brooks-Corey w/ Sgr > 0
! Flat above Sgr
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se

  if (Pc <= this%Pct) then ! Note, this function is degenerate
    Sw = 1d0
  else
    Se = (this%Pct/Pc)**this%lambda
    Sw = this%Swr + this%Se_span*Se
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP4Pc(this, Sw, Pc, dPc_dSw)
! Brooks-Corey w/ Sgr > 0
! Curve continues above Sgr
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se

  Se = (Sw - this%Swr) * this%dSe_dSw
  Pc = this%Pct * Se**this%lambda_nrec
  dPc_dSw = -this%dSe_dSw * Pc / (this%lambda*Se)
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP4Sw(this, Pc, Sw)
! Brooks-Corey w/ Sgr > 0
! Curve continues above Sgr
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se

  if (Pc <= this%Pct) then
    Sw = 1d0
  else
    Se = (this%Pct/Pc)**this%lambda
    Sw = this%Swr + this%Se_span*Se
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP5Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
!------- Linear model (A)
  if (Sw <= this%Swr) then
    Pc = this%Pcmax
    dPc_dSw = 0d0
  else if (Sw <= this%Sgr_comp) then
    Pc = this%Pct
    dPc_dSw = 0d0
  else
    dPc_dSw = (this%Pct-this%Pcmax)*this%dSe_dSw
    Pc = dPc_dSw*(Sw-this%Swr) + this%Pcmax
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP5Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se

  if (Pc <= 0d0) then
   Sw = 1d0
  else
    Se = (Pc-this%Pcmax)/(this%Pct-this%Pcmax)
    Sw = this%Swr + this%Se_span*Se
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP8Pc(this, Sw, Pc, dPc_dSw)
! vG w/ Sgr = 0
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se, Se_mrtrec

  if (Sw >= 1d0) then
    Pc = 0d0
    dPc_dSw = -huge(dPc_dSw) ! TODO finite difference limit
  else
    Se = (Sw - this%Swr) * this%dSe_dSw
    Se_mrtrec = Se**this%m_nrec
    Pc = this%Pct * (Se_mrtrec - 1d0)**this%n_rec
    dPc_dSw = (this%dSe_dSw*this%m_nrec*this%n_rec*Pc) / (Se-Se/Se_mrtrec)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP8Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se

  if (Pc <= 0d0) then
   Sw = 1d0
  else
    Se = ((Pc**this%n)/this%Pct + 1d0)**this%m_nrec
    Sw = this%Swr + this%Se_span*Se
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP9Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  ! Computes the capillary_pressure as a function of saturation
  ! based on experimental measurements and analyses done by Vauclin et al.
  ! as discussed by Moridis and Pruess, and the BRAGFLO V6.02 Requirements
  ! Document and Verification and Validation Plan, Sandia National Laboratories,
  ! Carlsbad, NM. ERMS #558659.  
  ! Moridis, G. J., and K. Pruess.  1992.  TOUGH Simulations of 
  ! Updegraff's Set of Fluid and Heat Flow Problems.  LBL-32611, ERMS# 138458.
  ! Berkeley, CA:  Lawrence Berkeley Laboratory.
  !
  ! Author: Heeho Park
  ! Date: 03/26/15
  PetscReal :: Sg
  PetscReal :: Se9

  PetscReal, parameter :: a = 3783.0145d0
  PetscReal, parameter :: b = 2.9d0
  PetscReal, parameter :: b_rec = 1d0/b

  if (Sw <= this%Swr) then
    Pc = 0d0
    dPc_dSw = 0d0
  else
    Sg = 1d0 - Sw
    Se9 = Sg/Sw
    Pc = a*Se9**b_rec
    dPc_dSw = -Pc/(b*Sw*Sg)
  endif
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP9Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se9

  PetscReal, parameter :: a = 3783.0145d0
  PetscReal, parameter :: b = 2.9d0
  PetscReal, parameter :: b_rec = 1d0/b

  if (Pc <= 0d0) then
    Sw = 1d0
  else
    Se9 = (Pc/a)**b
    Sw = 1d0 / (Se9+1d0)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP11Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  Pc = 0d0
  dPc_dSw = 0d0
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP11Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  Sw = 1d0
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP12Pc(this, Sw, Pc, dPc_dSw)
! BC w/ flat extension
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se
 
  if (Sw >= 1d0) then
    Pc = this%pct
  else
    Se = (Sw - this%Swr) * this%dSe_dSw
    Se = max(Se,this%Semin)
    Pc = this%pct * Se**this%lambda_nrec
  end if

  dPc_dSw = -this%dSe_dSw*Pc / (this%lambda*Se)

end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP12Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se

  if (Pc < this%Pct) then
    Sw = 1d0
  else
    Se = (this%pct/Pc)**this%lambda
    Sw = this%Swr + this%Se_span*Se
  end if
end subroutine

end module
