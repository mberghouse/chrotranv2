module Reaction_Sandbox_Kinetics_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module

  implicit none
  
  private
  
  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_kinetics_type
    ! Aqueous species
    PetscInt :: species_Aaq_id
    PetscInt :: species_Baq_id
    PetscInt :: species_Caq_id
    PetscInt :: species_Daq_id
    PetscInt :: species_Eaq_id
    PetscInt :: species_Faq_id

    PetscInt :: species_na_id
    PetscInt :: species_cl_id
    PetscInt :: species_fe2_id
    PetscInt :: species_fe3_id
    PetscInt :: species_cr6_id
    PetscInt :: species_cr3_id
    PetscInt :: species_h_id
    PetscInt :: species_o2_id

! Immobile species (e.g. biomass)
    PetscInt :: species_Xim_id
    PetscInt :: species_Yim_id
    character(len=MAXWORDLENGTH) :: species_name
    PetscInt :: species_id
    PetscReal :: rate_constant
    character(len=MAXWORDLENGTH) :: model
  contains
    procedure, public :: ReadInput => KineticsRead
    procedure, public :: Setup => KineticsSetup
    procedure, public :: Evaluate => KineticsEvaluate
    procedure, public :: Destroy => KineticsDestroy
  end type reaction_sandbox_kinetics_type

  public :: KineticsCreate

contains

! ************************************************************************** !

function KineticsCreate()
  ! 
  ! Allocates kinetics reaction object.
  ! 
  ! Author: Peter Lichtner
  ! Date: 10/24/2021

  implicit none
  
  class(reaction_sandbox_kinetics_type), pointer :: KineticsCreate

  allocate(KineticsCreate)
  KineticsCreate%species_name = ''
  KineticsCreate%species_id = 0
  KineticsCreate%rate_constant = 0.d0
  KineticsCreate%model = ''
  nullify(KineticsCreate%next)
      
end function KineticsCreate

! ************************************************************************** !

subroutine KineticsRead(this,input,option)
  !
  ! Reads input deck for kinetics reaction parameters (if any)
  !
  ! Author: Peter Lichtner
  ! Date: 10/24/2021
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_kinetics_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units
  
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,KINETICS')
    call StringToUpper(word)

    select case(trim(word))

      ! Example Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !   : begin user-defined input
      !     KINETICS
      !       MODEL Cr
      !       MODEL Cr-Fe
      !     END
      !   : end user defined input
      !   END
      !   ...
      ! END

! 5. Add case statement for reading variables.
      case('SPECIES_NAME')
! 6. Read the variable
        ! Read the character string indicating which of the primary species
        ! is being decayed.
        call InputReadWord(input,option,this%species_name,PETSC_TRUE)
! 7. Inform the user of any errors if not read correctly.
        call InputErrorMsg(input,option,'species_name', &
                           'CHEMISTRY,REACTION_SANDBOX,KINETICS')
! 8. Repeat for other variables
      case('RATE_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'rate_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,KINETICS')
        ! Read the units
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) then
          ! If units do not exist, assume default units of 1/s which are the
          ! standard internal PFLOTRAN units for this rate constant.
          input%err_buf = 'REACTION_SANDBOX,KINETICS,RATE CONSTANT UNITS'
          call InputDefaultMsg(input,option)
        else
          ! If units exist, convert to internal units of 1/s
          internal_units = 'unitless/sec'
          this%rate_constant = this%rate_constant * &
            UnitsConvertToInternal(word,internal_units,option)
        endif

      case('MODEL')
        call InputReadWord(input,option,this%model,PETSC_TRUE)
        call InputErrorMsg(input,option,'model', &
                          'CHEMISTRY,REACTION_SANDBOX,KINETICS')
      case default
        call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,KINETICS',option)
    end select
  enddo
  call InputPopBlock(input,option)
  
end subroutine KineticsRead

! ************************************************************************** !

subroutine KineticsSetup(this,reaction,option)
  ! 
  ! Sets up the kinetics reaction with hardwired parameters
  ! 
  ! Author: Peter Lichtner
  ! Date: 10/24/2021

  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_kinetics_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word

  ! Aqueous species
  word = 'Aaq'
  this%species_Aaq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Baq'
  this%species_Baq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Caq'
  this%species_Caq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Daq'
  this%species_Daq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Eaq'
  this%species_Eaq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Faq'
  this%species_Faq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)

  word = 'Na+'
  this%species_na_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Cl-'
  this%species_cl_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Fe++'
  this%species_fe2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Fe+++'
  this%species_fe3_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'CrO4--'
  this%species_cr6_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Cr+++'
  this%species_cr3_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'H+'
  this%species_h_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'O2(aq)'
  this%species_o2_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)

  ! Immobile species
  word = 'Xim'
  this%species_Xim_id = &
    GetImmobileSpeciesIDFromName(word,reaction%immobile,option)
  word = 'Yim'
  this%species_Yim_id = &
    GetImmobileSpeciesIDFromName(word,reaction%immobile,option)

end subroutine KineticsSetup

! ************************************************************************** !

subroutine KineticsEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar,reaction, &
                          option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Peter Lichtner
  ! Date: 10/24/2021
  ! 
  use Option_module

  use String_module
  use Input_Aux_module

  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(reaction_sandbox_kinetics_type) :: this

  type(input_type) :: input
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative

  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: volume                 ! m^3 bulk volume
  PetscReal :: porosity               ! m^3 pore / m^3 bulk volume
  PetscReal :: liquid_saturation      ! m^3 water / m^3 pore space
  PetscReal :: L_water                ! L water
  
  PetscReal :: Aaq, Baq, Caq, Daq, Eaq, Faq  ! mol/L water
  PetscReal :: Xim, Yim  ! mol/m^3 bulk volume
  PetscReal :: Rate, Rate1, Rate2, RateAB
  PetscReal :: RateA, RateB, RateC, RateD, RateE, RateF, RateX, RateY  ! mol/sec
  PetscReal :: stoichA, stoichB, stoichC, stoichD, stoichE, stoichF
  PetscReal :: stoichX, stoichY
  PetscReal :: k0, k, kr, k1, k2, kTST, kf, kb  ! units are problem specific
  PetscReal :: K_Aaq, K_Baq ! [mol/L water]
  PetscReal :: kfe, kcr0, kcr, n, pO2, oh, ph, RateFe, RateCr, RateFeCr
  PetscReal :: na,cl,fe2,fe3,cro4,cr,h
  PetscReal :: m_na, m_cl, m_fe2, m_fe3, m_cr6, m_cr3, m_h, m_o2
  PetscReal :: Totalfe2
  PetscReal :: gh ! activity coefficients
  PetscReal :: stoichfe2, stoichfe3, stoichh, stoicho2, stoichcr, stoichcro4
  PetscReal :: Koh, Ko2, Keq
  PetscReal :: Ratefe2, Ratefe3, Rateh, Rateo2, Ratecr6, Ratecr3
  PetscReal :: RateABs, dRAA, dRAB, dRBA, dRBB, dRCC
  character(len=MAXWORDLENGTH) :: word

  porosity = material_auxvar%porosity               ! [m^3 pore/m^3 bulk volume]
  liquid_saturation = global_auxvar%sat(iphase)     ! [m^3 water/m^3 pore]
  volume = material_auxvar%volume                   ! [m^3 bulk volume]
            ! multiply by 1.d3 to convert m^3 water -> L water
  L_water = porosity*liquid_saturation*volume*1.d3

! print *, 'params: ',L_water,porosity,liquid_saturation,volume
  
  Aaq = rt_auxvar%pri_molal(this%species_Aaq_id) ! [mol/L water]
  Baq = rt_auxvar%pri_molal(this%species_Baq_id) ! [mol/L water]
  Caq = rt_auxvar%pri_molal(this%species_Caq_id) ! [mol/L water]
  Daq = rt_auxvar%pri_molal(this%species_Daq_id) ! [mol/L water]
  Eaq = rt_auxvar%pri_molal(this%species_Eaq_id) ! [mol/L water]
  Faq = rt_auxvar%pri_molal(this%species_Faq_id) ! [mol/L water]

  m_na = rt_auxvar%pri_molal(this%species_na_id) ! [mol/L water]
  m_cl = rt_auxvar%pri_molal(this%species_cl_id) ! [mol/L water]
  m_fe2 = rt_auxvar%pri_molal(this%species_fe2_id) ! [mol/L water]
  m_fe3 = rt_auxvar%pri_molal(this%species_fe3_id) ! [mol/L water]
  m_cr6 = rt_auxvar%pri_molal(this%species_cr6_id) ! [mol/L water]
  m_cr3 = rt_auxvar%pri_molal(this%species_cr3_id) ! [mol/L water]
  m_h = rt_auxvar%pri_molal(this%species_h_id) ! [mol/L water]
  m_o2 = rt_auxvar%pri_molal(this%species_o2_id) ! [mol/L water]

  Totalfe2 = rt_auxvar%total(this%species_fe2_id,iphase) ! [mol/L water]
  
  Xim = rt_auxvar%immobile(this%species_Xim_id)     ! [mol/m^3 bulk volume]
  Yim = rt_auxvar%immobile(this%species_Yim_id)     ! [mol/m^3 bulk volume]

  ! initialize all rates to zero
  Rate = 0.d0
  RateA = 0.d0
  RateB = 0.d0
  RateC = 0.d0
  RateD = 0.d0
  RateE = 0.d0
  RateF = 0.d0
  RateX = 0.d0
  RateY = 0.d0

  Ratefe2 = 0.d0
  Ratefe3 = 0.d0
  Rateh = 0.d0
  Rateo2 = 0.d0

  ! stoichiometries
  ! reactants have negative stoichiometry
  ! products have positive stoichiometry
  stoichA = 0.d0
  stoichB = 0.d0
  stoichC = 0.d0
  stoichD = 0.d0
  stoichE = 0.d0
  stoichF = 0.d0
  stoichX = 0.d0
  stoichY = 0.d0
  
  ! kinetic rate constants
  k = 0.d0
  kr = 0.d0

  ! Monod half-saturation constants
  K_Aaq = 0.d0
  K_Baq = 0.d0

  word = this%model
! print *,'word = ',word

  gh = rt_auxvar%pri_act_coef(this%species_h_id) ! [mol/L water]

  select case(trim(word))
! models: Cr, Fe, Cr-Fe, Cu, Cu-Fe, Cr-Cu-Fe, Fe-CO2, Fe-Cr-CO2, ...
!   TST, A+B=AB(aq)
    case('Cr')

      print *,'Model Cr not implemented'
      stop

    case('Fe')
 
  !++++++++++++++++++++++++++++++++++++++++++
  ! Oxidation of Fe(II) Phreeqc ex9
  ! Fe2+ + H+ + 1/4 O2 → Fe3+ + 1/2 H2O
  !++++++++++++++++++++++++++++++++++++++++++

      k0 = 2.91d-9
      k1 = 1.33d12
      pO2 = 0.2
      Ko2 = 10.d0**(2.8983d0)
      Koh = 10.d0**(-13.9951d0)
 
      ph = -log10(gh * m_h)

!     oh = Koh / (gh * m_h)
      oh = Koh / m_h

      kfe = k0 + k1 * oh**2 * pO2

      pO2 = Ko2 * m_o2
!     RateFe = (k0 + k1 * oh**2 * pO2) * Totalfe2 * 365.d0 ! * L_water
      RateFe = (k0 + k1 * oh**2 * pO2) * m_fe2 * 365.d0 ! * L_water

      stoichfe2 = -1.d0
      stoichfe3 =  1.d0
      stoichh = -1.d0
      stoicho2 = -0.25d0

      Ratefe2 = stoichfe2 * RateFe
      Ratefe3 = stoichfe3 * RateFe
      Rateh = stoichh * RateFe
      Rateo2 = stoicho2 * RateFe

!     print *,'kinetics1: ',ph,oh,RateFe,m_fe2,m_fe3,m_cro4,kfe

  
  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second

  !Residual(this%species_na_id) = Residual(this%species_na_id)
  !Residual(this%species_cl_id) = Residual(this%species_cl_id)
  Residual(this%species_fe2_id) = Residual(this%species_fe2_id) - Ratefe2
  Residual(this%species_fe3_id) = Residual(this%species_fe3_id) - Ratefe3
  !Residual(this%species_cr6_id) = Residual(this%species_cr6_id) - Ratecr6
  !Residual(this%species_cr3_id) = Residual(this%species_cr3_id) - Ratecr3
  Residual(this%species_h_id) = Residual(this%species_h_id) - Rateh
  Residual(this%species_o2_id) = Residual(this%species_o2_id) - Rateo2


    case('Cr-Fe')
  
  !++++++++++++++++++++++++++++++++++++++++++
  ! Cr(VI) reduction by Fe(II): Fendorf & Li (1996)
  ! 3 Fe2+ + CrO2− + 8 H+ → Cr3+ + 3 Fe3+ + 4 H2O
  !++++++++++++++++++++++++++++++++++++++++++

      ph = -log10(gh * m_h)
      oh = 10.d0**(ph-13.9951d0)
      pO2 = 0.2d0

!     kfe = 8.e13 * oh**2 * pO2 / 60.d0
      kfe = (k0 + k1 * oh**2 * pO2) / 60.d0

      n = 0.6
      kcr = 56.3d0 * 10.d0**(3.d0*n) / 60.d0

      RateFe = kfe * m_fe2 * L_water             ! [mol/sec]
      RateCr = kcr * m_fe2**n * m_cr6 * L_water  ! [mol/sec]

!     following stoichiometries are not used
      stoichfe2 = -1.d0
      stoichfe3 = 1.d0
      stoichcro4 = -1.d0
      stoichcr = 1.d0
      stoichh = -1.d0
      stoicho2 = -0.25d0

      Ratefe2 = -1.d0 * RateFe
      Ratefe3 = 1.d0 * RateFe
      Ratecr6 = -1.d0 * RateCr
      Ratecr3 = 1.d0 * RateCr
      Rateh = -1.d0 * RateFe - 5.d0 * RateCr
      Rateo2 = -0.25d0 * RateFe + 0.75d0 * RateCr


  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second

  !Residual(this%species_na_id) = Residual(this%species_na_id)
  !Residual(this%species_cl_id) = Residual(this%species_cl_id)
  Residual(this%species_fe2_id) = Residual(this%species_fe2_id) - Ratefe2
  Residual(this%species_fe3_id) = Residual(this%species_fe3_id) - Ratefe3
  Residual(this%species_cr6_id) = Residual(this%species_cr6_id) - Ratecr6
  Residual(this%species_cr3_id) = Residual(this%species_cr3_id) - Ratecr3
  !Residual(this%species_h_id) = Residual(this%species_h_id) - Rateh
  !Residual(this%species_o2_id) = Residual(this%species_o2_id) - Rateo2

    case('Cr-Fe-H')

!++++++++++++++++++++++++++++++++++++++++++
! Cr(VI) reduction by Fe(II): Fendorf & Li (1996)
! 3 Fe2+ + CrO42− + 8 H+ → Cr3+ + 3 Fe3+ + 4 H2O
!++++++++++++++++++++++++++++++++++++++++++

  ph = -log10(gh * m_h)
  oh = 10.d0**(ph-13.9951d0)
  pO2 = 0.2d0

! kfe = 8.e13 * oh**2 * pO2 / 60.d0
! kfe = (k0 + k1 * oh**2 * pO2) / 60.d0

  n = 0.6d0

  kcr0 = 56.3
  kcr  = kcr0 * (10.d0**(3.d0*n))/60.d0

! print *,'kcr= ',kcr,kcr0 ! times out at first step when print statement turned on

! print *,'kcr= ',kcr0

! kcr = 59.204830823724784 ! calculated using mathematica

! RateFe = kfe * m_fe2 * L_water             ! [mol/sec]
  RateFeCr = kcr * m_fe2**n * m_cr6 * L_water  ! [mol/sec]

! print *, 'params2: ',L_water,porosity,liquid_saturation,volume,kcr,kcr0,RateFeCr

! following stoichiometries are not used
  stoichfe2 = -3.d0
  stoichfe3 = 3.d0
  stoichcro4 = -1.d0
  stoichcr = 1.d0
  stoichh = -8.d0
  stoicho2 = 0.d0

  Ratefe2 = -3.d0 * RateFeCr
  Ratefe3 = 3.d0 * RateFeCr
  Ratecr6 = -1.d0 * RateFeCr
  Ratecr3 = 1.d0 * RateFeCr
  Rateh = -8.d0 * RateFeCr
  Rateo2 = 0.d0 * RateFeCr


! NOTES
! 1. Always subtract contribution from residual
! 2. Units of residual are moles/second

!Residual(this%species_na_id) = Residual(this%species_na_id)
!Residual(this%species_cl_id) = Residual(this%species_cl_id)
Residual(this%species_fe2_id) = Residual(this%species_fe2_id) - Ratefe2
Residual(this%species_fe3_id) = Residual(this%species_fe3_id) - Ratefe3
Residual(this%species_cr6_id) = Residual(this%species_cr6_id) - Ratecr6
Residual(this%species_cr3_id) = Residual(this%species_cr3_id) - Ratecr3
Residual(this%species_h_id) = Residual(this%species_h_id)     - Rateh
!Residual(this%species_o2_id) = Residual(this%species_o2_id) - Rateo2

    case('Cu')

      print *,'Model Cu not implemented'
      stop

    case('Cu-Fe')

      print *,'Model Cu-Fe not implemented'
      stop

    case('Cr-Cu-Fe')

      print *,'Model Cr-Cu-Fe not implemented'
      stop

    case('TST')

!     Dissolution reaction for a single component using TST rate law
!     Keq = 1.d2
      Keq = 1.d0
      kTST = 5.e-5
!     Rate = -kTST * (Baq - 1.d0/Keq) * volume * L_water
      Rate =  kTST / Keq * (1.d0 - Keq * Baq) * L_water

    case('Elem')

!     A + B = C(aq), A + B = C(s)
      kf = 1.e-3
      kb = 1.e-3
      RateAB = (kf*Aaq*Baq - kb*Caq) * L_water

      Keq = 1.d0
      kTST = 5.d-2
      RateABs = -kTST * (1-Keq*Aaq*Baq)

!     print *,'rate: ',RateAB,RateABs,L_water,Aaq,Baq,Caq

    if (compute_derivative) then
      dRAA = kf * Baq * L_water
      dRAB = kf * Aaq * L_water
      dRBA = kf * Baq * L_water
      dRBB = kf * Aaq * L_water
      dRCC = kb * L_water
    endif

    case default
      print *,'Kinetics error msg.: ', this%model, ' Model not recognized. Stop!'
      stop
  end select


  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second  
  Residual(this%species_Aaq_id) = Residual(this%species_Aaq_id) - (-RateAB - RateABs)
  Residual(this%species_Baq_id) = Residual(this%species_Baq_id) - (-RateAB - RateABs)
  Residual(this%species_Caq_id) = Residual(this%species_Caq_id) - RateAB
  !Residual(this%species_Daq_id) = Residual(this%species_Daq_id) - RateD
  !Residual(this%species_Eaq_id) = Residual(this%species_Eaq_id) - RateE
  !Residual(this%species_Faq_id) = Residual(this%species_Faq_id) - RateF

#if 0
  Residual(this%species_Xim_id + reaction%offset_immobile) = &
    Residual(this%species_Xim_id + reaction%offset_immobile) - RateX
  Residual(this%species_Yim_id + reaction%offset_immobile) = &
    Residual(this%species_Yim_id + reaction%offset_immobile) - RateY
#endif

#if 0
  if (compute_derivative) then
    option%io_buffer = 'Reaction Sandbox Kinetics does not support analytical &
                       &derivatives. Stop!'
    call PrintErrMsg(option)
  endif
#endif

  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for Jacobian evaluation

    ! always add contribution to Jacobian
    ! units = (mol/sec)*(kg water/mol) = kg water/sec
    Jacobian(this%species_Aaq_id,this%species_Aaq_id) = &
    Jacobian(this%species_Aaq_id,this%species_Aaq_id) - &
      (-dRAA) * L_water * &
      rt_auxvar%aqueous%dtotal(this%species_Aaq_id,this%species_Aaq_id,iphase)

    Jacobian(this%species_Aaq_id,this%species_Baq_id) = &
    Jacobian(this%species_Aaq_id,this%species_Baq_id) - &
      (-dRAB) * L_water * &
       rt_auxvar%aqueous%dtotal(this%species_Aaq_id,this%species_Baq_id,iphase)

    Jacobian(this%species_Baq_id,this%species_Aaq_id) = &
    Jacobian(this%species_Baq_id,this%species_Aaq_id) - &
      (-dRBA) * L_water * &
      rt_auxvar%aqueous%dtotal(this%species_Baq_id,this%species_Aaq_id,iphase)

    Jacobian(this%species_Baq_id,this%species_Baq_id) = &
    Jacobian(this%species_Baq_id,this%species_Baq_id) - &
      (-dRBB) * L_water * &
      rt_auxvar%aqueous%dtotal(this%species_Baq_id,this%species_Baq_id,iphase)

    Jacobian(this%species_Caq_id,this%species_Caq_id) = &
    Jacobian(this%species_Caq_id,this%species_Caq_id) - &
      (dRCC) * L_water * &
      rt_auxvar%aqueous%dtotal(this%species_Caq_id,this%species_Caq_id,iphase)

  endif
  
end subroutine KineticsEvaluate

! ************************************************************************** !

subroutine KineticsDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Peter Lichtner
  ! Date: 10/11/2021
  !
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_kinetics_type) :: this

! deallocate(this%species_na_id)
! call DeallocateArray(this%irow)
! call DeallocateArray(this%icol)
! call DeallocateArray(this%stoich_row)

end subroutine KineticsDestroy

end module Reaction_Sandbox_Kinetics_class
