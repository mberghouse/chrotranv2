module Reaction_Sandbox_ArcticGHG_class

#include "petsc/finclude/petscsys.h"
  use petscsys

! 1. Change all references to "Example" as desired to rename the module and
!    and subroutines within the module.

  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

! 2. Add module variables here.  Note that one must use the PETSc data types
!    PetscInt, PetscReal, PetscBool to declare variables of type integer
!    float/real*8, and logical respectively.  E.g.,
!
! PetscReal, parameter :: formula_weight_of_water = 18.01534d0

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_arcticGHG_type
! 3. Add variables/arrays associated with new reaction.
    ! Aqueous species
    PetscInt :: species_CH2Oaq_id
    PetscInt :: species_O2aq_id
    PetscInt :: species_CO2aq_id
    PetscInt :: species_Daq_id
    ! Immobile species (e.g. biomass)
    !PetscInt :: species_Xim_id
    character(len=MAXWORDLENGTH) :: species_name
    PetscInt :: species_id
    PetscReal :: rate_constant
    PetscReal :: k_max
    PetscReal :: K_CH2Oaq_n
    PetscReal :: K_O2aq
    PetscReal :: I_CO2aq
    PetscReal :: yield
    PetscReal :: k_decay
    PetscReal :: n
    PetscBool :: molarity_units
    PetscReal, pointer :: stoich(:)
  contains
    procedure, public :: ReadInput => ArcticGHGRead
    procedure, public :: Setup => ArcticGHGSetup
    procedure, public :: Evaluate => ArcticGHGEvaluate
    procedure, public :: Destroy => ArcticGHGDestroy
  end type reaction_sandbox_arcticGHG_type

  public :: ArcticGHGCreate

contains

! ************************************************************************** !

function ArcticGHGCreate()
  !
  ! Allocates example reaction object.
  !
  ! Author: John Doe (replace in all subroutine headers with name of developer)
  ! Date: 00/00/00 (replace in all subroutine headers with current date)
  !

  implicit none

  class(reaction_sandbox_arcticGHG_type), pointer :: ArcticGHGCreate

! 4. Add code to allocate the object, initialize all variables to zero and
!    nullify all pointers. E.g.,
  allocate(ArcticGHGCreate)
  ArcticGHGCreate%species_name = ''
  ArcticGHGCreate%species_id = 0
  ArcticGHGCreate%rate_constant = 0.d0
  ArcticGHGCreate%k_max = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%K_CH2Oaq_n = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%K_O2aq = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%I_CO2aq = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%yield = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%k_decay = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%n = 1.d0
  ArcticGHGCreate%molarity_units = PETSC_TRUE
  nullify(ArcticGHGCreate%next)
  !DF: call MicrobialCreate

end function ArcticGHGCreate

! ************************************************************************** !

subroutine ArcticGHGRead(this,input,option)
  !
  ! Reads input deck for monod parameters
  !
  ! Author: David Fukuyama
  ! Date: 07/17/23
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_arcticGHG_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units, error_string
  PetscReal :: K_CH2Oaq

  K_CH2Oaq = UNINITIALIZED_DOUBLE
  error_string = 'CHEMISTRY,REACTION_SANDBOX,FLEXIBLE_BIODEGRADATION'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,ARCTICGHG')
    call StringToUpper(word)

    select case(trim(word))
      ! case('REACTION')
      !   ! remainder of string should be the reaction equation
      !   this%reaction = trim(adjustl(input%buf))
      !   ! set flag for error message
      !   if (len_trim(this%reaction) < 2) input%ierr = 1
      !   call InputErrorMsg(input,option,'reaction string', &
      !                       'CHEMISTRY,MICROBIAL_REACTION,REACTION')

      ! Temperature list
      ! Reaction rate list
      case('AQUEOUS_CONCENTRATION_UNITS')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case(word)
          case('MOLARITY')
            this%molarity_units = PETSC_TRUE
          case('MOLALITY')
            this%molarity_units = PETSC_FALSE
          case default
            call InputKeywordUnrecognized(input,word, &
                         trim(error_string)//&
                         'AQUEOUS_CONCENTRATION_UNITS',option)
        end select
      case('MAX_SPECIFIC_UTILIZATION_RATE')
        call InputReadDouble(input,option,this%k_max)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%k_max,'1/sec|mol/mol-sec', &
                                      trim(error_string)//','//word,option)
      case('CH2OAQ_HALF_SATURATION_CONSTANT')
        call InputReadDouble(input,option,K_CH2Oaq)
        call InputErrorMsg(input,option,word,error_string)
      case('O2AQ_HALF_SATURATION_CONSTANT')
        call InputReadDouble(input,option,this%K_O2aq)
        call InputErrorMsg(input,option,word,error_string)
      case('CO2AQ_MONOD_INHIBITION_CONSTANT')
        call InputReadDouble(input,option,this%I_CO2aq)
        call InputErrorMsg(input,option,word,error_string)
      case('YIELD')
        call InputReadDouble(input,option,this%yield)
        call InputErrorMsg(input,option,word,error_string)
      case('BIOMASS_DECAY_RATE_CONSTANT')
        call InputReadDouble(input,option,this%k_decay)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%k_decay,'1/sec', &
                                      trim(error_string)//','//word,option)
      case('HILL_EXPONENT')
        call InputReadDouble(input,option,this%n)
        call InputErrorMsg(input,option,word,error_string)

! ! 5. Add a case statement for reading variables.
!       case('SPECIES_NAME')
! ! 6. Read the variable.
!         ! Read the character string indicating which of the primary species
!         ! is being decayed.
!         call InputReadWord(input,option,this%species_name,PETSC_TRUE)
! ! 7. Inform the user of any errors if not read correctly.
!         call InputErrorMsg(input,option,'SPECIES_NAME', &
!                            'CHEMISTRY,REACTION_SANDBOX,ARCTICGHG')
! ! 8. Repeat for other variables.
!       case('RATE_CONSTANT')
!         ! Read the double precision rate constant.
!         call InputReadDouble(input,option,this%rate_constant)
!         ! Note the use of character variable 'word' instead of 'RATE_CONSTANT'
!         ! in the error message, as they are identical.
!         call InputErrorMsg(input,option,word, &
!                            'CHEMISTRY,REACTION_SANDBOX,ARCTICGHG')
!         ! Read the optional units and convert to internal
!         ! units of 1/s.
!         internal_units = 'unitless/sec'
!         call InputReadAndConvertUnits(input,this%rate_constant, &
!                                 internal_units,'CHEMISTRY,REACTION_SANDBOX,&
!                                 &ARCTICGHG,RATE_CONSTANT',option)
      case default
        call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,ARCTICGHG',option)
    end select
  enddo
  call InputPopBlock(input,option)

  error_string = ''
  if (Uninitialized(this%k_max)) then
    error_string = 'MAX_SPECIFIC_UTILIZATION_RATE,'
  endif
  if (Uninitialized(this%k_decay)) then
    error_string = trim(error_string) // 'BIOMASS_DECAY_RATE_CONSTANT,'
  endif
  if (Uninitialized(K_CH2Oaq)) then
    error_string = trim(error_string) // 'CH2OAQ_HALF_SATURATION_CONSTANT,'
  endif
  if (Uninitialized(this%K_O2aq)) then
    error_string = trim(error_string) // 'O2AQ_HALF_SATURATION_CONSTANT,'
  endif
  if (Uninitialized(this%I_CO2aq)) then
    error_string = trim(error_string) // 'CO2AQ_MONOD_INHIBITION_CONSTANT,'
  endif
  if (Uninitialized(this%yield)) then
    error_string = trim(error_string) // 'YIELD,'
  endif

  if (len_trim(error_string) > 0) then
    option%io_buffer = 'Reaction Sandbox FLEXIBLE_BIODEGRADATION has &
      &uninitialized parameters: ' &
      // error_string(1:len_trim(error_string)-1)
    call PrintErrMsg(option)
  endif

  this%K_CH2Oaq_n = K_CH2Oaq**this%n

end subroutine ArcticGHGRead

! ************************************************************************** !

subroutine ArcticGHGSetup(this,reaction,option)
  !
  ! Sets up the example reaction with parameters either read from the
  ! input deck or hardwired.
  !
  ! Author: David Fukuyama
  ! Date: 07/17/23
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_arcticGHG_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word

! 9. Add code to initialize.
  ! this%species_id = &
  !   GetPrimarySpeciesIDFromName(this%species_name,reaction,option)

  ! Aqueous species
  word = 'CH2OAq'
  this%species_CH2Oaq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'O2Aq'
  this%species_O2aq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'CO2Aq'
  this%species_CO2aq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  word = 'Daq'
  this%species_Daq_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)

  !Immobile species
  ! word = 'Xim'
  ! this%species_Xim_id = &
  !   GetImmobileSpeciesIDFromName(word,reaction%immobile,option)

  allocate(this%stoich(reaction%ncomp))
  this%stoich = 0.d0
  this%stoich(this%species_CH2Oaq_id) = -1.d0
  this%stoich(this%species_O2aq_id) = -1.d0
  this%stoich(this%species_CO2aq_id) = 1.d0
  this%stoich(this%species_Daq_id) = 0.d0
  !this%stoich(this%species_Xim_id+reaction%offset_immobile) = this%yield

end subroutine ArcticGHGSetup

! ************************************************************************** !

subroutine ArcticGHGEvaluate(this,Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_arcticGHG_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: volume               ! [m^3 bulk volume]
  PetscReal :: molality_to_molarity ! [kg water / L water]
  PetscReal :: CH2Oaq, O2aq, CO2aq, Daq   ! [mole / L water] or [mole / kg water]
!  PetscReal :: Xim                  ! [mole biomass / m^3 bulk volume]
  PetscReal :: I_r                  ! [mole rxn / m^3 bulk volume-sec]
  PetscReal :: I                    ! [mole rxn / sec]
  PetscReal :: dIdx(reaction%ncomp) ! mixed units
  PetscInt :: icomp, jcomp, Xim_offset

!  Xim_offset = this%species_Xim_id + reaction%offset_immobile
  volume = material_auxvar%volume
  if (this%molarity_units) then
    molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
  else
    molality_to_molarity = 1.d0
  endif

  CH2Oaq = rt_auxvar%pri_molal(this%species_CH2Oaq_id)*molality_to_molarity
  O2aq = rt_auxvar%pri_molal(this%species_O2aq_id)*molality_to_molarity
  CO2aq = rt_auxvar%pri_molal(this%species_CO2aq_id)*molality_to_molarity
  Daq = rt_auxvar%pri_molal(this%species_Daq_id)*molality_to_molarity
  !Xim = rt_auxvar%immobile(this%species_Xim_id)

  ! Anaerobic decomposition of C
  ! CH2O -> CO2 + CH4
  ! Aerobic decomposition of C
  ! CH2O + O2 -> CO2
  ! Methanogenesis of CO2
  ! CO2 + H2O -> CH4
  ! 

  I_r = this%k_max * CH2Oaq**this%n / (this%K_CH2Oaq_n + CH2Oaq**this%n) * &
                           O2aq / (this%K_O2aq + O2aq) * &
                           this%I_CO2aq / (this%I_CO2aq + CO2aq)
  I = I_r * volume

  ! Note: Always subtract contribution from residual
  do icomp = 1, reaction%ncomp
    Residual(icomp) = Residual(icomp) - this%stoich(icomp) * I
  enddo
  ! Add due to negative stoichiometry for decay
  !Residual(Xim_offset) = Residual(Xim_offset) + this%k_decay * Xim * volume

  if (compute_derivative) then
    ! Note: Always subtract contribution from Jacobian
    dIdx = 0.d0
    ! [mole rxn - kg water / sec] for CH2Oaq, O2aq, CO2aq
    dIdx(this%species_CH2Oaq_id) = &
      this%n * I / rt_auxvar%pri_molal(this%species_CH2Oaq_id) - &
      I / (this%K_CH2Oaq_n + CH2Oaq**this%n) * &
        this%n * CH2Oaq**this%n / rt_auxvar%pri_molal(this%species_CH2Oaq_id)
    dIdx(this%species_O2aq_id) = &
      I / rt_auxvar%pri_molal(this%species_O2aq_id) - &
      I / (this%K_O2aq + O2aq) * &
        O2aq / rt_auxvar%pri_molal(this%species_O2aq_id)
    dIdx(this%species_CO2aq_id) = &
      -1.d0 * I / (this%I_CO2aq + CO2aq) * &
        CO2aq / rt_auxvar%pri_molal(this%species_CO2aq_id)
    ! [mole rxn - m^3 bulk / sec] for Xim
    ! dIdx(Xim_offset) = I / Xim
    ! do icomp = 1, reaction%ncomp
    !   do jcomp = 1, reaction%ncomp
    !     Jacobian(icomp,jcomp) = Jacobian(icomp,jcomp) - &
    !       this%stoich(icomp) * dIdx(jcomp)
    !   enddo
    ! enddo
    ! Jacobian(Xim_offset,Xim_offset) = &
    !   Jacobian(Xim_offset,Xim_offset) + this%k_decay * volume
  endif

  ! Description of subroutine arguments:

  ! Residual - 1D array storing Residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entries in units kg water/sec
  !
  !  Jacobian [kg water/sec] * dc [mol/kg water] = -Res [mol/sec]
  !
  ! compute_derivative - Flag indicating whether analytical derivatives will
  !   be calculated.  The user must provide either the analytical derivatives
  !   or a numerical approximation unless always running with
  !   NUMERICAL_JACOBIAN defined within the NUMERICAL_METHODS TRANSPORT,
  !   NEWTON_SOLVER block of the input deck.  If the use of
  !   NUMERICAL_JACOBIAN is assumed, the user should provide an error
  !   message when compute_derivative is true.  E.g.,
  !
  !   if (compute_derivative) then
  !     option%io_buffer = 'NUMERICAL_JACOBIAN must be specified within &
  !       &the NEWTON_SOLVER block of NUMERICAL_METHODS TRANSPORT due to &
  !       &assumptions made in ExampleEvaluate.'
  !     call PrintErrMsg(option)
  !   endif
  !
  ! rt_auxvar - Object holding chemistry information (e.g., concentrations,
  !   activity coefficients, mineral volume fractions, etc.).  See
  !   reactive_transport_aux.F90.
  !
  !   Useful variables:
  !     rt_auxvar%total(:,iphase) - total component concentrations
  !                                 [mol/L water] for phase
  !     rt_auxvar%pri_molal(:) - free ion concentrations [mol/kg water]
  !     rt_auxvar%pri_act_coef(:) - activity coefficients for primary species
  !     rt_auxvar%aqueous%dtotal(:,iphase) - derivative of total component
  !                 concentration with respect to free ion [kg water/L water]
  !
  ! global_auxvar - Object holding information on flow (e.g., saturation,
  !   density, viscosity, temperature, etc)
  !
  !   Useful variables:
  !     global_auxvar%den(iphase) - liquid density [mol/m^3]
  !     global_auxvar%den_kg(iphase) - liquid density [kg/m^3]
  !     global_auxvar%sat(iphase) - liquid saturation [m^3 water/m^3 pore]
  !     global_auxvar%temp - temperature [C]
  !
  ! porosity - effective porosity of grid cell [m^3 pore/m^3 bulk]
  ! volume - volume of grid cell [m^3]
  ! reaction - Provides access to variable describing chemistry.  E.g.,
  !   reaction%ncomp - # chemical degrees of freedom (mobile and immobile)
  !   reaction%naqcomp - # chemical degrees of freedom on water
  !   reaction%primary_species_names(:) - names of primary species
  !
  ! option - Provides handle for controlling simulation, catching and
  !          reporting errors.
  
! 10. Add code for the Residual evaluation.

  ! Units of the Residual must be in moles/second.
  ! 1.d3 converts m^3 water -> L water
!   L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
!             material_auxvar%volume*1.d3
!   ! Always "subtract" the contribution from the Residual.
!   Residual(this%species_id) = Residual(this%species_id) - &
!     (-1.d0) * & ! negative stoichiometry
!     this%rate_constant * &  ! 1/sec
!     L_water * & ! L water
!     rt_auxvar%total(this%species_id,iphase) ! mol/L water

!   if (compute_derivative) then

! ! 11. If using an analytical Jacobian, add code for the Jacobian evaluation.

!     ! Always "add" the contribution to the Jacobian.
!     ! Units = (mol/sec)*(kg water/mol) = kg water/sec
!     Jacobian(this%species_id,this%species_id) = &
!     Jacobian(this%species_id,this%species_id) - &
!       (-1.d0) * & ! negative stoichiometry
!       this%rate_constant * & ! 1/sec
!       L_water * & ! L water
!       ! rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) =
!       !   derivative of total component concentration with respect to the
!       !   free ion concentration of the same species.
!       ! kg water/L water
!       rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase)

!   endif

end subroutine ArcticGHGEvaluate

! ************************************************************************** !

subroutine ArcticGHGDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !
  implicit none

  class(reaction_sandbox_arcticGHG_type) :: this

! 12. Add code to deallocate dynamic members of reaction_sandbox_arcticGHG_type.

end subroutine ArcticGHGDestroy

end module Reaction_Sandbox_ArcticGHG_class
