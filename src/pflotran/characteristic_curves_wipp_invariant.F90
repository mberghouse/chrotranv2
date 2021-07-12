module characteristic_curves_WIPP_invariant_module
#include "petsc/finclude/petscsys.h"

use petscsys ! Necessary for PETSC_TRUE / PETSC_FALSE
use option_module ! Necessary for global passing of PCT
use characteristic_curves_base_module ! Needed to define base type
implicit none

#define KRP_VG_range 1,8
#define KRP_BC_range 2,3,4

! **************************************************************************** !
! Common BRAGFLO Saturation Function type
! **************************************************************************** !

type, public, extends(sat_func_base_type) :: sf_WIPP_type
  private
    procedure(calc_Pc_type), pointer :: KPCPc
    procedure(calc_Sw_type), pointer :: KPCSw

    procedure(calc_Pc_type), pointer :: KRPPc
    procedure(calc_Sw_type), pointer :: KRPSw

    procedure(get_Sw_type) , pointer  :: KRPSw_inflection

    PetscInt  :: KRP                      ! KRP
    PetscInt  :: KPC                      ! KPC
    PetscReal :: Swr                      ! SWR
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

! In BRAGFLO, VG functions use the closed form Mualem relationship n = 1/(1-m)
! XLAM2, XLAM3, and XLAM5 are only necessary for relative permeability

!   Precalculated ratio between midpoint BC Pc and Pt for KRP1
    PetscReal :: PcmPct

! SF_WIPP_Parameters to calculate Pct
    PetscReal :: permeability
    PetscReal :: pct_a
    PetscReal :: pct_exp
    PetscReal :: pct
    PetscReal :: alpha
    PetscBool :: ignore_permeability 

! Derived parameters for linear unsaturated extension
    PetscReal :: Swj, Pcj
    PetscReal :: dPcj_dSwj, dSwj_dPcj

  contains
    procedure, public :: CapillaryPressure => SF_BRAGFLO_CapillaryPressure
    procedure, public :: Saturation        => SF_BRAGFLO_Saturation
    procedure, public :: SetPcmax          => SF_BRAGFLO_SetPcmax
    procedure, public :: SetPermeability   => SF_BRAGFLO_SetPermeability
end type

! **************************************************************************** !
! Function prototypes for function pointers
! **************************************************************************** !

abstract interface
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
  pure function get_Sw_type(this) result (Sw)
    import sf_WIPP_type
    class(sf_WIPP_type), intent(in) :: this
    PetscReal :: Sw
  end function
end interface

! **************************************************************************** !
! Public/Private Procedure Declarations
! **************************************************************************** !

! General BRAGFLO constructor
public  :: SFWIPPctor

! Implemented BRAGFLO KPC subroutines
private :: SFWIPPKPC1Pc , &
           SFWIPPKPC1Sw , &
           SFWIPPKPC2Pc , &
           SFWIPPKPC2Sw , &
           SFWIPPKPC6Pc , &
           SFWIPPKPC6Sw

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
           SFWIPPKRP8Pc , &
           SFWIPPKRP8Sw , &
           SFWIPPKRP9Pc , &
           SFWIPPKRP9Sw , &
           SFWIPPKRP11Pc , &
           SFWIPPKRP11Sw
!           SFWIPPKRP12Pc , &
!           SFWIPPKRP12Sw

contains 

! **************************************************************************** !
! BRAGFLO Saturation Function Constructor
! **************************************************************************** !

function SFWIPPctor(KRP, KPC, Swr, Sgr, alpha, expon, Pcmax, Pct_a, Pct_exp, &
                    ignore_pct) result (new)
  class(sf_WIPP_type), pointer :: new
  PetscInt, intent(in)  :: KRP, KPC
  PetscReal, intent(in) :: Swr, Sgr, alpha, expon, Pcmax, Pct_a, Pct_exp
  PetscBool, intent(in) :: ignore_pct
  PetscInt :: error

  ! Memory allocation
  allocate(new)
  if (.not. associated(new)) return ! Memory allocation failed

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
  case (11)
    new%KRPPc => SFWIPPKRP11Pc
    new%KRPSw => SFWIPPKRP11Sw
!  case (12)
!    new%KRPPc => SFWIPPKRP12Pc
!    new%KRPSw => SFWIPPKRP12Sw
  case default
                                        error = error + 1
  end select

  ! If KPC is in branch table, set function pointer, otherwise flag error
  select case(KPC)
  case (1)
    new%KPCPc => SFWIPPKPC1Pc
    new%KPCSw => SFWIPPKPC1Sw
  case (2)
    new%KPCPc => SFWIPPKPC2Pc
    new%KPCSw => SFWIPPKPC2Sw
  case (6)
    new%KPCPc => SFWIPPKPC6Pc
    new%KPCSw => SFWIPPKPC6Sw
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

  ! Pcmax must more particularly be above Pct
  if (Pcmax < 0d0)                      error = error + 32

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
  new%Pcmax = Pcmax
  
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

  ! Calculate beta, Swj, and Pcj for linear unsaturated extension to Pcmax
!  error = new%setPcmax(Pcmax)

  ! Copy Pct parameters
  ! TODO
  ! If ignore permeability, alpha must be set
  ! If .NOT. ignore_permeability, pct_a and pct_exp must be set
  ! Unknown if it ever flips back and forth.
  ! In theory, 

  new%alpha = alpha
  new%pct_a = pct_a
  new%pct_exp = pct_exp
  new%ignore_permeability = ignore_pct

  ! Setup default Pct
  if (alpha == 0d0) then 
    new%Pct = 1E6 ! Default 1 GPa for test function. 
  else
    new%Pct = 1/alpha
  end if
end function

! **************************************************************************** !

function SF_BRAGFLO_setPcmax(this,Pcmax) result (error)
  class(sf_WIPP_type), intent(inout) :: this
  PetscReal, intent(in) :: Pcmax
  PetscInt :: error

  ! Use the bisection method to satisfy C1 continuity for linear extensions

  PetscReal :: Swa, Swb, Pce, dPcj_dSwj

  ! Set lower saturation bracket set where P(Sa) = Pcmax
  call this%KRPSw(Pcmax, Swa)

  ! Set upper saturation limit (Sb) to be the the inflection point
  Swb = this%KRPSw_inflection()

  ! Confirm Pcmax is above minimum extrapolating from inflection point
  call this%KRPPc(Swb, Pce, dPcj_dSwj)
  if (Pcmax > Pce - dPcj_dSwj*Swb) then
    error = 0
    this%Pcmax = Pcmax

    do while (Swb-Swa > epsilon(this%Swj)) ! Tolerance interval epsilon
      this%Swj = (Swa+Swb)/2d0             ! Bisect bracket
      call this%KRPPc(this%Swj, this%Pcj, this%dPcj_dSwj)

      ! Residual error = Pcmax + dPj_dSj * Sj - Pf(Sj)
      Pce = Pcmax + this%dPcj_dSwj*this%Swj - this%Pcj
      if (Pce < 0d0) then ! Sde a Error is negative below inflection point
        Swa = this%Swj
      else
        Swb = this%Swj
      end if
    end do
    this%dSwj_dPcj = 1d0 / this%dPcj_dSwj
  else
    error = 1
  end if
end function

! **************************************************************************** !
! BRAGFLO Saturation Function Wrappers
! **************************************************************************** !

subroutine SF_BRAGFLO_SetPermeability(this, permeability)
  class(sf_WIPP_type)              :: this
  PetscReal, intent(in)            :: permeability

  ! TODO could do data validation, i.e. permeability > 0d0
  if (permeability /= this%permeability) then
    this%permeability = permeability
    this%pct = this%pct_a * permeability ** this%pct_exp
    ! TODO update the extension. Fast for Sj options, slow for Pcmax options
    ! call this%set_Sj(this%Sj)
    ! Skip the option%flow%pct_updated. This is bad practice as it is process/thread unsafe.
    ! It would be appropriate to pass permeability to CapillaryPressure
  end if
end subroutine

subroutine SF_BRAGFLO_CapillaryPressure(this, liquid_saturation, &
                     capillary_pressure, dpc_dsatl, option)
  class(sf_WIPP_type)           :: this
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
  class(sf_WIPP_type)           :: this
  PetscReal, intent(in)            :: capillary_pressure
  PetscReal, intent(out)           :: liquid_saturation, dsat_dpres
  type(option_type), intent(inout) :: option

  if (this%ignore_permeability) then
    this%pct = 1d0/this%alpha
  else
    if (.not. option%flow%pct_updated) then
      option%io_buffer = '!! this%pct has not been updated: &
                         &sf_WIPP_type. STOPPING.'
      call PrintErrMsg(option)
    end if
    option%flow%pct_updated = PETSC_FALSE
  end if

  ! Calls the function pointer KPCSw
  ! KPCSw in turn calls KRPSw as needed
  call this%KPCSw(capillary_pressure, liquid_saturation)
  dsat_dpres = 0d0 ! Or NaN, whatever
end subroutine

! **************************************************************************** !
! BRAGFLO KPC Subroutines
! **************************************************************************** !

pure subroutine SFWIPPKPC1Pc(this, Sw, Pc, dPc_dSw)
  ! Apply KRP everywhere
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  call this%KRPPc(Sw, Pc, dPc_dSw)
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKPC1Sw(this, Pc, Sw)
  ! Apply KRP everywhere
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Pc
  PetscReal, intent(out) :: Sw

  call this%KRPSw(Pc, Sw)
end subroutine

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
  ! Apply KRP     Pc above
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

pure subroutine SFWIPPKPC5Pc(this, Sw, Pc, dPc_dSw)
  ! Apply maximum Pc below Swr
  ! Apply linear  Pc below Sgr_comp
  ! Apply KRP     Pc above
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  if (Sw < this%Swr) then 
    Pc = this%Pcmax
    dPc_dSw = 0d0
  else if (Sw < this%Sgr_comp) then
    dPc_dSw = (this%Pct - this%Pcmax) * this%dSe2_dSw
    Pc = this%Pcmax + dPc_dSw*(Sw - this%Swr)
  else
    call this%KRPPc(Sw, Pc, dPc_dSw)
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKPC6Pc(this, Sw, Pc, dPc_dSw)
  ! Apply linear Pc below Sjr
  ! Apply KRP    Pc above
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw

  if (Sw <= this%Swj) then
    if (Sw <= 0d0) then
      Pc = this%Pcmax
    else
      Pc = this%Pcmax - this%dPcj_dSwj*Sw
      dPc_dSw = this%dPcj_dSwj
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
! BRAGFLO KRP Subroutines
! **************************************************************************** !

pure subroutine SFWIPPKRP1Pc(this, Sw, Pc, dPc_dSw)
!------- Modified van Genuchten/Parker model (SGR>=0)
! Author: Heeho Park; Modified by Jennifer Frederick
! Date: 11/17/16; Modified 04/26/2017
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se2, Pct

  Pct = this%Pct

  if (Sw >= this%Swr) then
    Se2 = (Sw-this%Swr) * this%dSe2_dSw
    Se2 = MIN(Se2,1d0)
    Pc = Pct*this%PcmPct*(Se2**this%m_nrec-1d0)**this%m_comp
  else
    Pc = 0d0
    dPc_dSw = 0d0
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
  PetscReal :: Se2
  PetscReal :: Po

  if (Pc <= 0d0) then
    Sw = 1d0
  else
    Po = this%Pct * this%PcmPct
    Se2 = (((Pc**this%n)/Po) + 1d0)**(-this%m)
    Sw = this%Swr + Se2 * this%Se2_span
  end if
end subroutine

! **************************************************************************** !

pure subroutine SFWIPPKRP2Pc(this, Sw, Pc, dPc_dSw)
!------- Original Brooks-Corey model
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se1

  if (Sw <= this%Swr) then
    Pc = 0d0
    dPc_dSw = 0d0
  else
    Se1     = (Sw-this%Swr)*this%dSe1_dSw
    Pc      = this%Pct*(Se1**this%lambda_nrec)
    dPc_dSw = -this%dSe1_dSw*Pc/(this%lambda*Se1)
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
!------- Modified Brooks-Corey model
  PetscReal :: Se2

  if (Sw <= this%Sgr_comp) then
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

pure subroutine SFWIPPKRP4Pc(this, Sw, Pc, dPc_dSw)
  class(sf_WIPP_type), intent(in) :: this
  PetscReal, intent(in)  :: Sw
  PetscReal, intent(out) :: Pc, dPc_dSw
  PetscReal :: Se2
!------- Brooks-Corey model; modified non-wetting phase,
!-------   original wetting phase
! Note, Se2 can be greater than 1d0 in this model
  if (Sw > this%Swr) then
    Se2 = (Sw - this%Swr)*this%dSe2_dSw
    Pc = this%Pct*(Se2**this%lambda_nrec)
    dPc_dSw = -this%dSe2_dSw*Pc/(this%lambda*Se2)
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

end module
