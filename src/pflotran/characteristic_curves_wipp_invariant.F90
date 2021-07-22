module characteristic_curves_WIPP_invariant_module
#include "petsc/finclude/petscsys.h"

use petscsys ! Necessary for PETSC_TRUE / PETSC_FALSE
use option_module ! Necessary for global passing of PCT
use characteristic_curves_base_module ! Needed to define base type
implicit none

#define KRP_VG_range 1,8
#define KRP_BC_range 2,3,4,12

! **************************************************************************** !
! Common BRAGFLO Saturation Function type
! **************************************************************************** !

type, public, extends(sat_func_base_type) :: sf_WIPP_type
  private
    procedure(set_Pct_type), public, pointer :: setPct
    procedure(set_Swj_type), pointer :: setSwj
    procedure(calc_Pc_type), pointer :: KPCPc
    procedure(calc_Sw_type), pointer :: KPCSw

    procedure(calc_Pc_type), pointer :: KRPPc
    procedure(calc_Sw_type), pointer :: KRPSw

! KRP and KPC enums. could be stored locally, but it is so far unnecessary
! Branching statements should be evaluated at ctor time with function pointers
!
!   PetscInt  :: KRP                      ! KRP
!   PetscInt  :: KPC                      ! KPC

! Some memory could be saved with union/equivalence statements
! But, there are very few saturation_function objects
    PetscReal :: Swr                      ! SWR or SOCZRO
    PetscReal :: Sgr                      ! SGR
    PetscReal :: Sgr_comp
    PetscReal :: Se1_span                 ! 1/SETERM
    PetscReal :: dSe1_dSw                 ! SETERM
    PetscReal :: Se2_span                 ! 1/SETERM2
    PetscReal :: dSe2_dSw                 ! SETERM2
    PetscReal :: lambda                   ! 1/XLAM1
    PetscReal :: lambda_nrec              ! -XLAM1
    PetscReal :: m                        ! XLAM4
    PetscReal :: m_nrec                   ! XLAM6
    PetscReal :: m_comp                   ! XLAM7
    PetscReal :: n                        ! 1/XLAM7
    PetscReal :: n_rec                    ! -XLAM7
!   PetscReal :: Pcmax                    ! PCFIX
    PetscReal :: Semin                    ! SOCEFFMIN

! In BRAGFLO, VG functions use the closed form Mualem relationship n = 1/(1-m)
! XLAM2, XLAM3, and XLAM5 are only necessary for relative permeability

!   Precalculated ratio between midpoint BC Pc and Pt for KRP1
    PetscReal :: PcmPct

! Parameters for efficient Pct
    PetscReal :: pct
    PetscReal :: permeability
    PetscReal :: pct_a
    PetscReal :: pct_exp
    PetscReal :: alpha

! Derived parameters for unsaturated extensions
    PetscReal :: Swj, Pcj
    PetscReal :: dPcj_dSwj, dSwj_dPcj

  contains
    procedure, public :: CapillaryPressure => SF_BRAGFLO_CapillaryPressure
    procedure, public :: Saturation        => SF_BRAGFLO_Saturation
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
  subroutine set_pct_type(this, K)
    import sf_WIPP_type
    class(sf_WIPP_type), intent(inout) :: this
    PetscReal, intent(in)  :: K
  end subroutine
end interface

! **************************************************************************** !
! Public/Private Procedure Declarations
! **************************************************************************** !

! General BRAGFLO constructor
public  :: SFWIPPctor

! Implemented BRAGFLO KPC subroutines
private :: SFWIPPKPC1Swj, &
           SFWIPPKPC2Pc , &
           SFWIPPKPC2Sw , &
           SFWIPPKPC2Swj, &
           SFWIPPKPC6Pc , &
           SFWIPPKPC6Sw , &
           SFWIPPKPC6Swj

! Implemented BRAGFLO KRP Pc subroutines
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
           SFWIPPKRP8Sw , &
           SFWIPPKRP9Pc , &
           SFWIPPKRP9Sw , &
           SFWIPPKRP11Pc , &
           SFWIPPKRP11Sw , &
           SFWIPPKRP12Pc , &
           SFWIPPKRP12Sw

contains 

! **************************************************************************** !
! BRAGFLO Saturation Function Constructor
! **************************************************************************** !

function SFWIPPctor(KRP, KPC, Swr, Sgr, expon, Pct_ignore, Pct_alpha, &
                    Pct_expon, Pcmax, Swj, Smin, Semin) result (new)
  class(sf_WIPP_type), pointer :: new
  PetscInt, intent(inout)  :: KRP, KPC
  PetscReal, intent(in) :: Swr, Sgr, expon, Pct_alpha, Pct_expon, Swj, Pcmax
  PetscReal, intent(in) :: Smin, Semin
  PetscBool, intent(in) :: Pct_ignore
  PetscInt :: error

  ! Memory allocation
  allocate(new)
  if (.not. associated(new)) return ! Memory allocation failed

  ! TODO confirm PCT model.
  ! Except for KRP 9, if PCT_A, capilary pressure is always 0, as per KRP 11
  if (Pct_alpha == 0d0 .AND. KRP /= 9) then
    KRP = 11
  end if

  ! Data validation
                                        error = 0
  ! If KRP is in branch table, set function pointer, otherwise flag error
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
  case (11) ! Cavity model, Pc is always 0, nothing else matters
    new%KRPPc => SFWIPPKRP11Pc
    new%KRPSw => SFWIPPKRP11Sw
    new%KPCPc  => new%KRPPc
    new%KPCSw  => new%KRPSW
    new%setPct => SF_BRAGFLO_IgnorePermeability
    return
  case (12)
    new%KRPPc => SFWIPPKRP12Pc
    new%KRPSw => SFWIPPKRP12Sw
  case default
                                        error = error + 1
  end select

  ! If KPC is in branch table, set function pointer, otherwise flag error
  select case(KPC)
  case (1) ! KPC 1 points directly to KRP
    new%KPCPc  => new%KRPPc
    new%KPCSw  => new%KRPSW
    new%setSwj => SFWIPPKPC1Swj
  case (2)
    new%KPCPc  => SFWIPPKPC2Pc
    new%KPCSw  => SFWIPPKPC2Sw
    new%setSwj => SFWIPPKPC2Swj
  case (6)
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

  ! Assign canonical parameters
  new%Swr = Swr
  new%Sgr = Sgr

  ! Override Swr for KRP 12 only
  if (KRP == 12) then ! Note, BRAGFLO code does not match the BRAGFLO UM
    new%Swr = Smin - Semin ! Swr is BRAGFLO SOCZRO
    new%Semin = Semin
  end if
  
  ! Depending on the KRP option, canonical expon is either VG m or BC lambda
  select case(KRP)
  case (KRP_BC_range) ! Brooks-Corey types
    new%lambda = expon
    new%m = new%lambda/(1d0+new%lambda)
  case (KRP_VG_range) ! Van Genuchten types
    new%m = expon
    new%lambda = new%m/(1d0-new%m)
  case default        ! Unused parameters
    new%lambda = 0d0
    new%m = 0d0
  end select

  ! Depending on the PCT option, pct_alpha is either alpha or pct_a
  if (Pct_ignore) then
    new%setPct => SF_BRAGFLO_IgnorePermeability
    new%Pct = pct_alpha ! Pct permanently set to alpha
  else
    new%setPct => SF_BRAGFLO_SetPermeability
    new%pct_a = pct_alpha
    new%pct_exp = pct_expon
    new%Pct = 1E6 ! Initialize to 1 MPa for test function. 
  end if

  ! Derived 1 or 2 residual saturation parameters
  new%Sgr_comp = 1d0 - new%Sgr

  new%Se1_span = 1d0 - new%Swr
  new%dSe1_dSw = 1d0 / new%dSe1_dSw

  new%Se2_span = new%Sgr_comp - new%Swr
  new%dSe2_dSw = 1d0 / new%Se2_span

  ! Brooks-Corey and Van Genuchten-Mualem Derived Parameters
  new%lambda_nrec = -1d0/new%lambda
  new%PcmPct = 0.5d0**new%lambda_nrec
  new%m_nrec = -1d0 / new%m
  new%m_comp =  1d0 - new%m
  new%n      =  1d0 / new%m_comp
  new%n_rec  =  1d0 / new%n

  ! Initialize unsaturated extension
  error = new%setSwj(Swj)

end function

! **************************************************************************** !
! BRAGFLO Saturation Function Wrappers
! **************************************************************************** !

subroutine SF_BRAGFLO_CapillaryPressure(this, liquid_saturation, &
                     capillary_pressure, dpc_dsatl, option)
  class(sf_WIPP_type)              :: this
  PetscReal, intent(in)            :: liquid_saturation
  PetscReal, intent(out)           :: capillary_pressure, dpc_dsatl
  type(option_type), intent(inout) :: option

  ! Calls the function pointer KPCPc
  ! KPCPc in turn calls KRPPc as needed
  call this%KPCPc(liquid_saturation, capillary_pressure, dpc_dsatl)
end subroutine

! **************************************************************************** !

subroutine SF_BRAGFLO_Saturation(this, capillary_pressure, liquid_saturation, &
                                 dsat_dpres, option)
  class(sf_WIPP_type)              :: this
  PetscReal, intent(in)            :: capillary_pressure
  PetscReal, intent(out)           :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  ! Calls the function pointer KPCSw
  ! KPCSw in turn calls KRPSw as needed
  call this%KPCSw(capillary_pressure, liquid_saturation)
  dsat_dpres = 0d0 ! Or NaN, whatever
end subroutine

! **************************************************************************** !
! BRAGFLO PCT Subroutines
! **************************************************************************** !

subroutine SF_BRAGFLO_SetPermeability(this, permeability)
  class(sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: permeability
  PetscInt :: error
  
  if (permeability /= this%permeability) then
    this%permeability = permeability
    this%pct = this%pct_a * permeability ** this%pct_exp
    ! Refresh unsaturated extension
    error = this%setSwj(this%Swj)
  end if
end subroutine

! **************************************************************************** !

subroutine SF_BRAGFLO_IgnorePermeability(this, permeability)
  class(sf_WIPP_type), intent(inout)   :: this
  PetscReal, intent(in) :: permeability
end subroutine

! **************************************************************************** !
! BRAGFLO KPC Subroutines
! **************************************************************************** !

function SFWIPPKPC1Swj(this,Swj) result (error)
  class (sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Swj
  PetscInt :: error
  
  error = 0
end function

! **************************************************************************** !

pure subroutine SFWIPPKPC2Pc(this, Sw, Pc, dPc_dSw)
  ! Apply maximum Pc below Swj
  ! Apply KRP     Pc above
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  if (Sw <= this%Swj) then
    Pc = this%Pcmax
    dPc_dSw = 0d0
  else
    call this%KRPPc(Sw, Pc, dPc_dSw)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKPC2Sw(this, Pc, Sw)
  ! Apply maximum Pc below Swj
  ! Apply KRP     Pc abovecall this%KRPPc(Swj, this%Pcj, this%dPcj_dSwj)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw

  if (Pc >= this%Pcmax) then
    Sw = this%Swr
  else
    call this%KRPSw(Pc, Sw)
  end if
end subroutine

! **************************************************************************** !

function SFWIPPKPC2Swj(this,Swj) result (error)
  class (sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Swj
  PetscInt :: error

  call this%KRPSw(this%Pcmax,this%Swj)
  this%Pcj = this%Pcmax
  this%dPcj_dSwj = 0d0
  error = 0 
end function

! **************************************************************************** !

pure subroutine SFWIPPKPC6Pc(this, Sw, Pc, dPc_dSw)
  ! Apply linear Pc below Swr
  ! Apply KRP    Pc above
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  if (Sw <= this%Swj) then
    dPc_dSw = this%dPcj_dSwj
    if (Sw <= 0d0) then
      Pc = this%Pcmax
    else
      Pc = this%Pcmax + this%dPcj_dSwj*Sw
    end if
  else
    call this%KRPPc(Sw, Pc, dPc_dSw)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKPC6Sw(this, Pc, Sw)
  ! Apply linear Sw below Pcj
  ! Apply KRP    Sw above
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw

  if (Pc >= this%Pcj) then
    if (Pc >= this%Pcmax) then
      Sw = 0d0
    else
      Sw = (this%Pcmax - Pc) * this%dSwj_dPcj
    end if
  else
    call this%KRPSw(Pc, Sw)
  end if
end subroutine

! **************************************************************************** !

function SFWIPPKPC6Swj(this,Swj) result (error)
  class(sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Swj
  PetscInt :: error

  if (Swj > this%Swr) then
    error = 0
    this%Swj = Swj
    call this%KRPPc(Swj, this%Pcj, this%dPcj_dSwj)
    this%Pcmax = this%Pcj - this%dPcj_dSwj * Swj
    this%dSwj_dPcj = 1d0 / this%dPcj_dSwj
  else
    error = 1
  end if
end function

! **************************************************************************** !
! BRAGFLO KRP Subroutines
! **************************************************************************** !

pure subroutine SFWIPPKRP1Pc(this, Sw, Pc, dPc_dSw)
!------- Modified van Genuchten/Parker model (SGR>=0)
! Author: Heeho Park; Modified by Jennifer Frederick
! Date: 11/17/16; Modified 04/26/2017
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se2

  if (Sw > this%Swr) then
    if (Sw > this%Sgr_comp) then
      Pc = 0d0
    else
      Se2 = (Sw-this%Swr) * this%dSe2_dSw
      Pc = this%Pct*this%PcmPct*(Se2**this%m_nrec-1d0)**this%m_comp
    end if
  else
    Pc = 0d0
  end if

  ! TODO analytic derivatives
  dPc_dSw = 0d0
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP1Sw(this, Pc, Sw)
!------- Modified van Genuchten/Parker model (SGR>=0)
! Author: Heeho Park; Modified by Jennifer Frederick
! Date: 11/17/16; Modified 04/26/2017
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se2

  if (Pc <= 0d0) then
    Sw = 1d0
  else
    Se2 = ( Pc**this%n / (this%Pct * this%PcmPct) + 1d0)**(-this%m)
    Sw = this%Swr + this%Se2_span * Se2
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP2Pc(this, Sw, Pc, dPc_dSw)
!------- Original Brooks-Corey model
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se1

  if (Sw > this%Swr) then
    Se1     = (Sw-this%Swr)*this%dSe1_dSw
    Pc      = this%Pct * Se1**this%lambda_nrec
    dPc_dSw = -this%dSe1_dSw*Pc / (this%lambda*Se1)
  else
    Pc = 0d0
    dPc_dSw = 0d0
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP2Sw(this, Pc, Sw)
!------- Original Brooks-Corey model
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se1, Pct

  Pct = this%Pct

  if (Pc <= Pct) then
    Sw = 1d0
  else
    Se1     = (Pc/Pc)**(-this%lambda)
    Sw      = this%Swr + this%Se1_span*Se1
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP3Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se2

  if (Sw >= this%Sgr_comp) then
    Pc = this%Pct
    dPc_dSw = 0d0
  else if (Sw > this%Swr) then
    Se2 = (Sw - this%Swr)*this%dSe2_dSw
    Pc  = this%Pct*(Se2**this%lambda_nrec)
    dPc_dSw = -this%dSe2_dSw*Pc/(this%lambda*Se2)
  else
    Pc = 0d0
    dPc_dSw = 0d0
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP3Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se2, Pct

  Pct = this%Pct*(this%Se1_span/this%Se2_span)**this%lambda_nrec

  if (Pc <= Pct) then
    Sw = 1d0
  else
    Se2     = (Pct/Pc)**this%lambda
    Sw      = this%Swr + this%Se2_span*Se2
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP4Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se2
!------- Brooks-Corey model; modified non-wetting phase,
!-------   original wetting phase
! Note, Se2 can be greater than 1d0 in this model
  if (Sw > this%Swr) then
    Se2 = (Sw - this%Swr) * this%dSe2_dSw
    Pc = this%Pct * Se2**this%lambda_nrec
    dPc_dSw = -this%dSe2_dSw * Pc / (this%lambda*Se2)
  else
    Pc = 0d0
    dPc_dSw = 0d0
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP4Sw(this, Pc, Sw)
!------- Original Brooks-Corey model
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se2, Pct

  Pct = this%Pct*(this%Se1_span/this%Se2_span)**this%lambda_nrec

  if (Pc <= Pct) then
    Sw = 1d0
  else
    Se2     = (Pct/Pc)**this%lambda
    Sw      = this%Swr + this%Se2_span*Se2
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
    dPc_dSw = (this%Pct-this%Pcmax)*this%dSe2_dSw
    Pc = dPc_dSw*(Sw-this%Swr) + this%Pcmax
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP5Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se2

  if (Pc <= 0d0) then
   Sw = 1d0
  else
    Se2 = (Pc-this%Pcmax)/(this%Pct-this%Pcmax)
    Sw = this%Swr + this%Se2_span*Se2
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP8Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se1, Se1_mrtrec
!------- Original van Genuchten/Parker model (SGR=0)

  if (Sw > this%Swr .and. Sw < 1d0) then
    Se1 = (Sw - this%Swr) * this%dSe1_dSw
    Se1_mrtrec = Se1**this%m_nrec
    Pc = this%Pct * (Se1_mrtrec - 1d0)**this%n_rec
    dPc_dSw = (this%dSe1_dSw*this%m_nrec*this%n_rec*Pc) / (Se1-Se1/Se1_mrtrec)
  else
    Pc = 0d0
    dPc_dSw = 0d0
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP8Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se1

  if (Pc <= 0d0) then
   Sw = 1d0
  else
    Se1 = ((Pc**this%n)/this%Pct + 1d0)**this%m_nrec
    Sw = this%Swr + this%Se1_span*Se1
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
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se21
  ! Extra-modified Brooks-Corey
  ! The sign of Semin is inconsistent between BRAGFLO and the BRAGFLO UM 
  ! Szero is calculated in the constructor

  Se21 = (Sw - this%Swr) * this%dSe1_dSw
  Se21 = max(min(Se21,1d0),this%Semin)

  Pc = this%pct * Se21**this%lambda_nrec
  dPc_dSw = -this%dSe1_dSw*Pc / (this%lambda*Se21)
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP12Sw(this, Pc, Sw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw
  PetscReal :: Se21

  if (Pc < this%Pct) then
    Sw = 1d0
  else
    Se21 = (this%pct/Pc)**this%lambda
    Sw = this%Swr + this%Se1_span*Se21
  end if
end subroutine

end module
