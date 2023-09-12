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
    PetscInt :: species_CH2O_id, species_O2_id, species_CO2_id
    PetscInt :: species_NO3_id, species_N2_id, species_HCO3_id
    PetscInt :: species_MNO2_id, species_MN_id
    PetscInt :: species_FEOH3_id, species_FE_id
    PetscInt :: species_SO4_id, species_H2S_id
    PetscInt :: species_CH4_id
    PetscInt :: species_D_id

    PetscInt :: respiration_id, methanogenesis_id
    PetscInt :: methane_oxidation_id
    ! Immobile species (e.g. biomass)
    !PetscInt :: species_Xim_id
    character(len=MAXWORDLENGTH) :: species_name
    PetscInt :: species_id
    PetscReal :: rate_constant
    PetscReal :: k_max
    PetscInt :: nrxn
    PetscReal :: K_CH2Oaq
    PetscReal :: K_CH2Oaq_methano
    PetscReal :: K_O2aq
    PetscReal :: I_CO2aq
    PetscReal :: yield
    PetscReal :: k_decay
    PetscReal :: k_CH2O
    PetscReal :: n
    PetscBool :: molarity_units
    PetscBool :: respiration
    PetscReal :: k_o2
    PetscReal :: k_c
    PetscReal :: r_o2_max
    PetscReal :: ae_o2
    PetscReal :: slope_o2
    PetscReal :: intercept_o2
    PetscReal :: alpha_o2
    PetscReal :: k_in_o2
    PetscBool :: methanogenesis
    PetscReal :: k_ch4
    PetscReal :: k_ch4ox
    PetscReal :: k_ch4ox_o2
    PetscReal :: r_ch4_max
    PetscReal :: r_ch4ox_max
    PetscReal :: ae_ch4
    PetscReal :: ae_ch4ox
    PetscReal :: slope_ch4
    PetscReal :: intercept_ch4
    PetscReal :: alpha_ch4
    PetscReal :: slope_ch4ox
    PetscReal :: intercept_ch4ox
    PetscReal :: alpha_ch4ox
    PetscReal :: k_in_ch4
    PetscBool :: methane_oxidation

    PetscReal, pointer :: stoich(:,:)
    PetscReal, pointer :: I_r(:)
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
  ArcticGHGCreate%nrxn = 0.d0
  ArcticGHGCreate%k_max = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%K_CH2Oaq = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%K_O2aq = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%I_CO2aq = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%yield = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%k_decay = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%n = 1.d0
  ArcticGHGCreate%molarity_units = PETSC_TRUE
  
  ArcticGHGCreate%respiration = PETSC_FALSE
  ArcticGHGCreate%r_o2_max = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%k_o2 = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%k_c = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%ae_o2 = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%slope_o2 = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%intercept_o2 = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%alpha_o2 = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%k_in_o2 = UNINITIALIZED_DOUBLE

  ArcticGHGCreate%methanogenesis = PETSC_FALSE
  ArcticGHGCreate%r_ch4_max = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%K_CH2Oaq_methano = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%ae_ch4 = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%slope_ch4 = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%intercept_ch4 = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%alpha_ch4 = UNINITIALIZED_DOUBLE
  ArcticGHGCreate%k_in_ch4 = UNINITIALIZED_DOUBLE
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

  character(len=MAXWORDLENGTH) :: word, internal_units, error_string, keyword, subkeyword
  PetscReal :: K_CH2Oaq
  PetscInt :: nrxn

  K_CH2Oaq = UNINITIALIZED_DOUBLE
  error_string = 'ARCTIC_GHG'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
           'CHEMISTRY,REACTION_SANDBOX,ARCTIC_GHG')
    call StringToUpper(word)

    select case(trim(word))
      case('REACTION')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit

          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'CONSTANT', &
                 'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_AQUEOUS')
          call StringToUpper(word)

          select case(trim(word))
            case('RESPIRATION','AEROBIC_RESPIRATION')
              this%respiration = PETSC_TRUE
              this%nrxn = this%nrxn + 1
              ! call InputReadWord(input,option,this%name_aqueous,PETSC_TRUE)
              ! call InputErrorMsg(input,option,'PARTICLE_NAME_AQ', &
              !        'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
              call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'CONSTANT', &
                       'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_AQUEOUS')
                call StringToUpper(word)
                select case(trim(word))
                  case('HALF_SATURATION_CONSTANT_O2','MONOD_CONSTANT_O2')
                    call InputReadDouble(input,option,this%k_o2)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_o2, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                 &RESPIRATION,HALF_SATURATION_CONSTANT_O2',option)
                  case('HALF_SATURATION_CONSTANT_C','MONOD_CONSTANT_C')
                    call InputReadDouble(input,option,this%k_CH2Oaq)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_CH2Oaq, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &RESPIRATION,HALF_SATURATION_CONSTANT_C',option)
                  case('RMAX')
                    call InputReadDouble(input,option,this%r_o2_max)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%r_o2_max, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &RESPIRATION,ALPHA',option)
                  case('ACTIVATION_ENERGY')
                    call InputReadDouble(input,option,this%ae_o2)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%ae_o2, &
                                       'J/mol','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &RESPIRATION,ACTIVATION_ENERGY',option)
                  case('MONOD_SLOPE')
                    call InputReadDouble(input,option,this%slope_o2)
                    call InputErrorMsg(input,option,word,error_string)
                  case('MONOD_INTERCEPT')
                    call InputReadDouble(input,option,this%intercept_o2)
                    call InputErrorMsg(input,option,word,error_string)
                  case('ALPHA')
                    call InputReadDouble(input,option,this%alpha_o2)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%alpha_o2, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &RESPIRATION,ALPHA',option)
                  case default
                    call InputKeywordUnrecognized(input,word, &
                           'CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,RESPIRATION',option)
                end select
              enddo
            case('METHANOGENESIS')
              this%nrxn = this%nrxn + 1
              this%methanogenesis = PETSC_TRUE
              call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'CONSTANT', &
                       'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_AQUEOUS')
                call StringToUpper(word)
                select case(trim(word))
                  case('HALF_SATURATION_CONSTANT_C','MONOD_CONSTANT_C')
                    call InputReadDouble(input,option,this%k_CH2Oaq_methano)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_CH2Oaq_methano, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                 &METHANOGENESIS,HALF_SATURATION_CONSTANT_C',option)
                  case('HALF_SATURATION_CONSTANT','MONOD_CONSTANT')
                    call InputReadDouble(input,option,this%k_ch4)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_ch4, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &METHANOGENESIS,HALF_SATURATION_CONSTANT',option)
                  case('RMAX')
                    call InputReadDouble(input,option,this%r_ch4_max)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%r_ch4_max, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &METHANOGENESIS,ALPHA',option)
                  case('ACTIVATION_ENERGY')
                    call InputReadDouble(input,option,this%ae_ch4)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%ae_ch4, &
                                       'J/mol','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &METHANOGENESIS,ACTIVATION_ENERGY',option)
                  case('MONOD_SLOPE')
                    call InputReadDouble(input,option,this%slope_ch4)
                    call InputErrorMsg(input,option,word,error_string)
                  case('MONOD_INTERCEPT')
                    call InputReadDouble(input,option,this%intercept_ch4)
                    call InputErrorMsg(input,option,word,error_string)
                  case('ALPHA')
                    call InputReadDouble(input,option,this%alpha_ch4)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%alpha_ch4, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &METHANOGENESIS,ALPHA',option)
                  case('INHIBITION_CONSTANT')
                    call InputReadDouble(input,option,this%k_in_o2)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_in_o2, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &METHANOGENESIS,INHIBITION_CONSTANT',option)
                  case default
                    call InputKeywordUnrecognized(input,word, &
                           'CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,METHANOGENESIS',option)
                end select
              enddo
            case('METHANE_OXIDATION')
              this%nrxn = this%nrxn + 1
              this%methane_oxidation = PETSC_TRUE
              call InputPushBlock(input,option)
              do
                call InputReadPflotranString(input,option)
                if (InputError(input)) exit
                if (InputCheckExit(input,option)) exit

                call InputReadCard(input,option,word)
                call InputErrorMsg(input,option,'CONSTANT', &
                       'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_AQUEOUS')
                call StringToUpper(word)
                select case(trim(word))
                  case('HALF_SATURATION_CONSTANT_CH4','MONOD_CONSTANT_CH4')
                    call InputReadDouble(input,option,this%k_ch4ox)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_ch4ox, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                 &METHANE_OXIDATION,MONOD_CONSTANT_CH4',option)
                  case('HALF_SATURATION_CONSTANT_O2','MONOD_CONSTANT_O2')
                    call InputReadDouble(input,option,this%k_ch4ox_o2)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%k_ch4ox_o2, &
                                       'mol/L','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                 &METHANE_OXIDATION,MONOD_CONSTANT_CH4',option)
                  case('RMAX')
                    call InputReadDouble(input,option,this%r_ch4ox_max)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%r_ch4ox_max, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &METHANE_OXIDATION,RMAX',option)
                  case('ACTIVATION_ENERGY')
                    call InputReadDouble(input,option,this%ae_ch4ox)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%ae_ch4ox, &
                                       'J/mol','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &METHANE_OXIDATION,ACTIVATION_ENERGY',option)
                  case('MONOD_SLOPE')
                    call InputReadDouble(input,option,this%slope_ch4ox)
                    call InputErrorMsg(input,option,word,error_string)
                  case('MONOD_INTERCEPT')
                    call InputReadDouble(input,option,this%intercept_ch4ox)
                    call InputErrorMsg(input,option,word,error_string)
                  case('ALPHA')
                    call InputReadDouble(input,option,this%alpha_ch4ox)
                    call InputErrorMsg(input,option,word,error_string)
                    call InputReadAndConvertUnits(input,this%alpha_ch4ox, &
                                       'mol/L-s','CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,REACTION,&
                                                &METHANE_OXIDATION,ALPHA',option)
                  case default
                    call InputKeywordUnrecognized(input,word, &
                           'CHEMISTRY,REACTION_SANDBOX,ARCTICGHG,METHANE_OXIDATION',option)
                end select
              enddo
          end select
              call InputPopBlock(input,option)
        end do
        call InputPopBlock(input,option)

      case default
        call InputKeywordUnrecognized(input,word, &
               'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
               &DECAY_AQUEOUS,CONSTANT',option)
        end select
  end do
  call InputPopBlock(input,option)



  ! error_string = ''
  ! if (Uninitialized(this%k_max)) then
  !   error_string = 'MAX_SPECIFIC_UTILIZATION_RATE,'
  ! endif
  ! if (Uninitialized(this%k_decay)) then
  !   error_string = trim(error_string) // 'BIOMASS_DECAY_RATE_CONSTANT,'
  ! endif
  ! if (Uninitialized(K_CH2Oaq)) then
  !   error_string = trim(error_string) // 'CH2OAQ_HALF_SATURATION_CONSTANT,'
  ! endif
  ! if (Uninitialized(this%K_O2aq)) then
  !   error_string = trim(error_string) // 'O2AQ_HALF_SATURATION_CONSTANT,'
  ! endif
  ! if (Uninitialized(this%I_CO2aq)) then
  !   error_string = trim(error_string) // 'CO2AQ_MONOD_INHIBITION_CONSTANT,'
  ! endif
  ! if (Uninitialized(this%yield)) then
  !   error_string = trim(error_string) // 'YIELD,'
  ! endif

  ! if (len_trim(error_string) > 0) then
  !   option%io_buffer = 'Reaction Sandbox FLEXIBLE_BIODEGRADATION has &
  !     &uninitialized parameters: ' &
  !     // error_string(1:len_trim(error_string)-1)
  !   call PrintErrMsg(option)
  ! endif

  ! this%K_CH2Oaq_n = K_CH2Oaq**this%n

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
  this%respiration_id = 1
  this%methanogenesis_id = 2
  this%methane_oxidation_id = 3

  allocate(this%stoich(3,reaction%ncomp))

  ! word = 'CH2OAq'
  ! this%species_CH2Oaq_id = &
  !   GetPrimarySpeciesIDFromName(word,reaction,option)
  ! word = 'O2Aq'
  ! this%species_O2aq_id = &
  !   GetPrimarySpeciesIDFromName(word,reaction,option)
  ! word = 'CO2Aq'
  ! this%species_CO2aq_id = &
  !   GetPrimarySpeciesIDFromName(word,reaction,option)
  ! word = 'Daq'
  ! this%species_Daq_id = &
  !   GetPrimarySpeciesIDFromName(word,reaction,option)

  this%stoich = 0.d0

  word = 'CH2O'
  this%species_CH2O_id = &
    GetPrimarySpeciesIDFromName(word,reaction,option)
  if (this%respiration) then
    word = 'O2(aq)'
    this%species_O2_id = &
      GetPrimarySpeciesIDFromName(word,reaction,option)
    word = 'CO2'
    this%species_CO2_id = &
      GetPrimarySpeciesIDFromName(word,reaction,option)
    this%stoich(this%respiration_id,this%species_CH2O_id) = -1.d0
    this%stoich(this%respiration_id,this%species_O2_id) = -1.d0
    this%stoich(this%respiration_id,this%species_CO2_id) = 1.d0
  endif

  if (this%methanogenesis) then
    word = 'CH4'
    this%species_CH4_id = &
      GetPrimarySpeciesIDFromName(word,reaction,option)
    if (.not. this%respiration) then
      word = 'CO2'
      this%species_CO2_id = &
        GetPrimarySpeciesIDFromName(word,reaction,option)
    endif
    this%stoich(this%methanogenesis_id,this%species_CH2O_id) = -1.d0
    this%stoich(this%methanogenesis_id,this%species_CH4_id) = 0.5d0
    this%stoich(this%methanogenesis_id,this%species_CO2_id) = 0.5d0
  endif

  if (this%methane_oxidation) then
    if (.not. this%methanogenesis) then
      word = 'CH4'
      this%species_CH4_id = &
        GetPrimarySpeciesIDFromName(word,reaction,option)
    endif
    if (.not. this%respiration) then
      word = 'CO2'
      this%species_CO2_id = &
        GetPrimarySpeciesIDFromName(word,reaction,option)
    endif
    this%stoich(this%methane_oxidation_id,this%species_O2_id) = -1.d0
    this%stoich(this%methane_oxidation_id,this%species_CH4_id) = -1.d0
    this%stoich(this%methane_oxidation_id,this%species_CO2_id) = 1.d0
  endif

  !Immobile species
  ! word = 'Xim'
  ! this%species_Xim_id = &
  !   GetImmobileSpeciesIDFromName(word,reaction%immobile,option)


  ! this%stoich = 0.d0
  ! this%stoich(nrxn,this%species_CH2Oaq_id) = -1.d0
  ! this%stoich(nrxn,this%species_O2aq_id) = -1.d0
  ! this%stoich(nrxn,this%species_CO2aq_id) = 1.d0
  ! this%stoich(nrxn,this%species_Daq_id) = 0.d0
  !this%stoich(this%species_Xim_id+reaction%offset_immobile) = this%yield

end subroutine ArcticGHGSetup

! ************************************************************************** !

subroutine ArcticGHGEvaluate(this,Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar, &
                           reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: David Fukuyama
  ! Date: 07/17/23
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
  PetscReal :: volume_L_water
  PetscReal :: molality_to_molarity ! [kg water / L water]
  PetscReal :: CH2O,O2,CO2,NO3,MN,FE,SO4,CH4  ! [mole / L water] or [mole / kg water]
!  PetscReal :: Xim                  ! [mole biomass / m^3 bulk volume]
!  PetscReal :: I_r(nrxn)                  ! [mole rxn / m^3 bulk volume-sec]
  PetscReal :: I                    ! [mole rxn / sec]
  PetscReal :: I_r2, I2
  PetscReal :: dIdx(reaction%ncomp) ! mixed units
  PetscInt :: icomp, jcomp, Xim_offset, irxn!, nrxn
  PetscReal :: T, R

  allocate(this%I_r(3))

  T = global_auxvar%temp + 273.15d0 ! K
  R = IDEAL_GAS_CONSTANT            ! J/mol-K

!  Xim_offset = this%species_Xim_id + reaction%offset_immobile
  volume = material_auxvar%volume   ! m^3 bulk volume
  volume_L_water = volume * material_auxvar%porosity * 1.d3 * global_auxvar%sat(1)! L pore water
  if (this%molarity_units) then
    molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
 else
    
    molality_to_molarity = 1.d0
  endif

  CH2O = rt_auxvar%pri_molal(this%species_CH2O_id) * molality_to_molarity ! mol/L

  ! Respiration
  if (this%respiration) then
    O2   = rt_auxvar%pri_molal(this%species_O2_id) * molality_to_molarity
    CO2  = rt_auxvar%pri_molal(this%species_CO2_id) * molality_to_molarity
    if (Initialized(this%alpha_o2)) then
      this%r_o2_max = this%alpha_o2 * exp(-1.d0*this%ae_o2/(R*T)) 
      ! alpha_o2: [mol/L-s]
      ! ae_o2   : [J/mol]
      ! R       : [J/mol-K]
      ! T       : [K]
    endif
    if (Initialized(this%slope_o2)) then
      this%k_o2 = this%intercept_o2 + (T-273.15d0) * this%slope_o2
    endif
    ! I_r units: [mol/L-s] * [mol/L]/([mol/L]+[mol/L]) * [mol/L]/([mol/L]+[mol/L])
    this%I_r(this%respiration_id) = this%r_o2_max * CH2O / (this%k_CH2O + CH2O) * O2 / (this%k_o2 + O2)
    ! Residual units: [mol/s]
    Residual(this%species_O2_id) = Residual(this%species_O2_id) - this%I_r(this%respiration_id) * &
                                   volume_L_water * this%stoich(this%respiration_id,this%species_O2_id)
    Residual(this%species_CO2_id) = Residual(this%species_CO2_id) - this%I_r(this%respiration_id) * &
                                   volume_L_water * this%stoich(this%respiration_id,this%species_CO2_id)
    Residual(this%species_CH2O_id) = Residual(this%species_CH2O_id) - this%I_r(this%respiration_id) * &
                                   volume_L_water * this%stoich(this%respiration_id,this%species_CH2O_id)
  endif

  ! Methanogenesis
  if (this%methanogenesis) then
    CH4 = rt_auxvar%pri_molal(this%species_CH4_id) * molality_to_molarity
    if (.not. this%respiration) then
      O2   = rt_auxvar%pri_molal(this%species_O2_id) * molality_to_molarity
    endif
    if (Initialized(this%alpha_ch4)) then
      this%r_ch4_max = this%alpha_ch4 * exp(-1.d0*this%ae_ch4/(R*T))
    endif
    if (Initialized(this%slope_o2)) then
      this%k_ch4 = this%intercept_ch4 + (T-273.15d0) * this%slope_ch4
    endif
    this%I_r(this%methanogenesis_id) = this%r_ch4_max * CH2O / (this%k_CH2Oaq_methano + CH2O) * &
                                       this%k_in_o2 / (this%k_in_o2 + O2)
    Residual(this%species_CH4_id) = Residual(this%species_CH4_id) - this%I_r(this%methanogenesis_id) * &
                                   volume_L_water * this%stoich(this%methanogenesis_id,this%species_CH4_id)
    Residual(this%species_CO2_id) = Residual(this%species_CO2_id) - this%I_r(this%methanogenesis_id) * &
                                   volume_L_water * this%stoich(this%methanogenesis_id,this%species_CO2_id)
    Residual(this%species_CH2O_id) = Residual(this%species_CH2O_id) - this%I_r(this%methanogenesis_id) * &
                                   volume_L_water * this%stoich(this%methanogenesis_id,this%species_CH2O_id)
  endif
  !print *, 'ch4 reaction rate: ',this%r_ch4_max
  ! Methane oxidation
  if (this%methane_oxidation) then
    if (.not. this%methanogenesis) then
      CH4 = rt_auxvar%pri_molal(this%species_CH4_id) * molality_to_molarity
    endif
    if (.not. this%respiration) then
      O2  = rt_auxvar%pri_molal(this%species_O2_id) * molality_to_molarity
      CO2 = rt_auxvar%pri_molal(this%species_CO2_id) * molality_to_molarity
    endif
    if (Initialized(this%alpha_ch4ox)) then
      this%r_ch4_max = this%alpha_ch4ox * exp(-1.d0*this%ae_ch4ox/(R*T))
    endif
    if (Initialized(this%slope_o2)) then
      this%k_ch4 = this%intercept_ch4ox + (T-273.15d0) * this%slope_ch4ox
    endif
    this%I_r(this%methane_oxidation_id) = this%r_ch4ox_max * CH4 / (this%k_ch4ox + CH4) * &
                                          O2 / (this%k_ch4ox_o2 + O2)
    Residual(this%species_CH4_id) = Residual(this%species_CH4_id) - this%I_r(this%methane_oxidation_id) * &
                                   volume_L_water * this%stoich(this%methane_oxidation_id,this%species_CH4_id)
    Residual(this%species_O2_id) = Residual(this%species_O2_id) - this%I_r(this%methane_oxidation_id) * &
                                   volume_L_water * this%stoich(this%methane_oxidation_id,this%species_O2_id)
    Residual(this%species_CO2_id) = Residual(this%species_CO2_id) - this%I_r(this%methane_oxidation_id) * &
                                   volume_L_water * this%stoich(this%methane_oxidation_id,this%species_CO2_id)
  endif

  ! do irxn = 1, this%nrxn
  !   do icomp = 1, reaction%ncomp
  !      Residual(icomp) = Residual(icomp) - this%I_r(irxn) * volume * this%stoich(irxn,icomp)
  !   end do
  ! end do

  ! CH2O = rt_auxvar%pri_molal(this%species_CH2O_id)*molality_to_molarity
  ! O2 = rt_auxvar%pri_molal(this%species_O2_id)*molality_to_molarity
  ! CO2 = rt_auxvar%pri_molal(this%species_CO2_id)*molality_to_molarity
  !Xim = rt_auxvar%immobile(this%species_Xim_id)

  ! Anaerobic decomposition of C
  ! CH2O -> CO2 + CH4
  ! Aerobic decomposition of C
  ! CH2O + O2 -> CO2
  ! Methanogenesis of CO2
  ! CO2 + H2O -> CH4
  ! 

!   I_r

!   I_r = this%k_max * CH2O**this%n / (this%K_CH2O_n + CH2O**this%n) * &
!                            O2 / (this%K_O2 + O2aq)! * &
! !                           this%I_CO2aq / (this%I_CO2aq + CO2aq)
!   I = I_r * volume + 1.0 - 1.0

!   I_r2 = this%k_max/2.3 * CH2O**this%n / (this%K_CH2O_n + CH2O**this%n) * &
!                            this%I_CO2aq / (this%I_CO2aq + O2aq)
!   I2 = I_r2 * volume

!   ! Note: Always subtract contribution from residual
!   do icomp = 1, reaction%ncomp
!     Residual(icomp) = Residual(icomp) - this%stoich(icomp) * I
!   enddo
 
!   Residual(this%species_CH2O_id) = Residual(this%species_CH2O_id) + 1.d0 * I2
!   Residual(this%species_CO2aq_id) = Residual(this%species_CO2aq_id) - 0.4d0 * I2
!   Residual(this%species_Daq_id) = Residual(this%species_Daq_id) -0.6d0 * I2

!   ! Add due to negative stoichiometry for decay
!   !Residual(Xim_offset) = Residual(Xim_offset) + this%k_decay * Xim * volume

!   if (compute_derivative) then
!     ! Note: Always subtract contribution from Jacobian
!     dIdx = 0.d0
!     ! [mole rxn - kg water / sec] for CH2O, O2aq, CO2aq
!     dIdx(this%species_CH2O_id) = &
!       this%n * I / rt_auxvar%pri_molal(this%species_CH2O_id) - &
!       I / (this%K_CH2O_n + CH2O**this%n) * &
!         this%n * CH2O**this%n / rt_auxvar%pri_molal(this%species_CH2O_id)
!     dIdx(this%species_O2aq_id) = &
!       I / rt_auxvar%pri_molal(this%species_O2aq_id) - &
!       I / (this%K_O2aq + O2aq) * &
!         O2aq / rt_auxvar%pri_molal(this%species_O2aq_id)
!     dIdx(this%species_CO2aq_id) = &
!       -1.d0 * I / (this%I_CO2aq + CO2aq) * &
!         CO2aq / rt_auxvar%pri_molal(this%species_CO2aq_id)

!     ! dIdx(this%species(CH2O_id) = &
!     !   this%n * I / rt_auxvar%pri_molal(this%species_CH2O_id) - &
!     !   I / (this%K_CH2O_n + CH2O**this%n) * &


!     ! dIdx(this%
!     ! [mole rxn - m^3 bulk / sec] for Xim
!     ! dIdx(Xim_offset) = I / Xim
!     ! do icomp = 1, reaction%ncomp
!     !   do jcomp = 1, reaction%ncomp
!     !     Jacobian(icomp,jcomp) = Jacobian(icomp,jcomp) - &
!     !       this%stoich(icomp) * dIdx(jcomp)
!     !   enddo
!     ! enddo
!     ! Jacobian(Xim_offset,Xim_offset) = &
!     !   Jacobian(Xim_offset,Xim_offset) + this%k_decay * volume
!   endif

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
