module Reaction_Sandbox_BioTH_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_bioth_type

    !ID number of bioparticle species
    PetscInt :: species_Vaq_id ! Aqueous species
    PetscInt :: species_Vim_id ! Immobile species

    !Name of bioparticle species
    character(len=MAXWORDLENGTH) :: name_aqueous
    character(len=MAXWORDLENGTH) :: name_immobile

    !Decay rates (Temperature model)
    PetscReal :: logDref_aqueous
    PetscReal :: Tref_aqueous
    PetscReal :: zT_aqueous
    PetscReal :: nAq_aqueous

    PetscReal :: logDref_adsorbed
    PetscReal :: Tref_adsorbed
    PetscReal :: zT_adsorbed
    PetscReal :: nAq_adsorbed

    !Decay rates (Constant)
    PetscReal :: decay_aqueous
    PetscReal :: decay_adsorbed

    !Attachment rate (Filtration Model)
    PetscReal :: diam_collector
    PetscReal :: diam_particle
    PetscReal :: hamaker_constant
    PetscReal :: BOLTZMANN_CONSTANT = 1.380649d-23 !J/K
                 !(Not found on pflotran_constants.f90)
    PetscReal :: density_particle
    PetscReal :: alpha_efficiency

    !Attachment/Detachment rates (Constant)
    PetscReal :: rate_attachment
    PetscReal :: rate_detachment

    !Debug?
    PetscBool :: debug_option
	
	character(len=MAXWORDLENGTH) :: name_D_mobile   ! This the mobile food in the reaction
    character(len=MAXWORDLENGTH) :: name_D_immobile ! This the immobile food in the reaction
    character(len=MAXWORDLENGTH) :: name_C     ! This is for Cr(VI) in the reaction
    character(len=MAXWORDLENGTH) :: name_B ! This is for the biomass in the reaction
    character(len=MAXWORDLENGTH) :: name_I ! This is for the alcohol in the reaction
    character(len=MAXWORDLENGTH) :: name_X ! This is for the biocide in the reaction
    character(len=MAXWORDLENGTH) :: name_biomineral ! This is for the dummny bio mineral
	character(len=MAXWORDLENGTH) :: name_O2 ! This is for the dummny bio mineral
	character(len=MAXWORDLENGTH) :: name_CO2 ! This is for the dummny bio mineral

    PetscInt :: B_id
    PetscInt :: C_id
    PetscInt :: D_immobile_id  ! Immobile food
    PetscInt :: D_mobile_id    ! Mobile food
    PetscInt :: I_id
    PetscInt :: X_id
    PetscInt :: biomineral_id
    PetscInt :: O2_id
	PetscInt :: CO2_id

    ! Decay and inhibition parameters in our sophisticated model
    PetscReal :: background_concentration_B    ! Minimum background concentration of the biomass
    PetscReal :: rate_B_1
    PetscReal :: rate_B_2
    PetscReal :: rate_C
    PetscReal :: rate_D
    PetscReal :: rate_D_i
    PetscReal :: rate_D_m
    PetscReal :: inhibition_B
    PetscReal :: inhibition_C
    PetscReal :: monod_D
    PetscReal :: inhibition_I
    PetscReal :: mass_action_B
    PetscReal :: mass_action_CD
    PetscReal :: mass_action_X
    PetscReal :: stoichiometric_C
    PetscReal :: stoichiometric_D_1
    PetscReal :: stoichiometric_D_2
    PetscReal :: exponent_B
    PetscReal :: density_B
    PetscReal :: beta
    PetscReal :: alpha
    PetscReal :: beta_vel
    PetscReal :: alpha_vel

    PetscReal :: k    ! Maximum respiration rate
    PetscReal :: K_O  ! Half-saturation constant for oxygen

  contains
    procedure, public :: ReadInput => BioTH_Read
    procedure, public :: Setup => BioTH_Setup
    procedure, public :: Evaluate => BioTH_React
	procedure, public :: UpdateKineticState => BioTH_KineticState
    procedure, public :: Destroy => BioTH_Destroy

  end type reaction_sandbox_bioth_type

  public :: BioTH_Create

contains

! ************************************************************************** !

function BioTH_Create()
  !
  ! Allocates particle transport variables.
  !
  ! Author: Edwin Saavedra C.
  ! Date: 10/01/2020
  !

  implicit none

  class(reaction_sandbox_bioth_type), pointer :: BioTH_Create

  allocate(BioTH_Create)

  !ID number of bioparticle species
  BioTH_Create%species_Vaq_id = 0
  BioTH_Create%species_Vim_id = 0

  !Name of bioparticle species
  BioTH_Create%name_aqueous = ''
  BioTH_Create%name_immobile = ''

  !Decay rates (Temperature model)
  BioTH_Create%logDref_aqueous = 0.d0
  BioTH_Create%Tref_aqueous = 0.d0
  BioTH_Create%zT_aqueous = 0.d0
  BioTH_Create%nAq_aqueous = 0.d0

  BioTH_Create%logDref_adsorbed = 0.d0
  BioTH_Create%Tref_adsorbed = 0.d0
  BioTH_Create%zT_adsorbed = 0.d0
  BioTH_Create%nAq_adsorbed = 0.d0

  !Decay rates (Constant)
  BioTH_Create%decay_aqueous = -1.d0
  BioTH_Create%decay_adsorbed = -1.d0

  !Filtration Model
  BioTH_Create%diam_collector = 0.d0
  BioTH_Create%diam_particle = 0.d0
  BioTH_Create%hamaker_constant = 0.d0
  BioTH_Create%density_particle = 0.d0
  BioTH_Create%alpha_efficiency = 1.d0

  !Attachment rates
  BioTH_Create%rate_attachment = -1.d0
  BioTH_Create%rate_detachment = 0.d0

  !Attachment rates
  BioTH_Create%debug_option = .False.
  BioTH_Create%name_D_mobile = ''
  BioTH_Create%name_D_immobile = ''
  BioTH_Create%name_C = ''
  BioTH_Create%name_B = ''
  BioTH_Create%name_I = ''
  BioTH_Create%name_X = ''
  BioTH_Create%name_biomineral = ''

  BioTH_Create%D_mobile_id = 0
  BioTH_Create%D_immobile_id = 0
  BioTH_Create%C_id = 0
  BioTH_Create%B_id = 0
  BioTH_Create%I_id = 0
  BioTH_Create%X_id = 0
  BioTH_Create%biomineral_id = 0

  BioTH_Create%stoichiometric_D_1 = 0.d0
  BioTH_Create%rate_D = 0.d0
  BioTH_Create%rate_C = 0.d0
  BioTH_Create%inhibition_C = 0.d0
  BioTH_Create%rate_B_2 = 0.d0
  BioTH_Create%rate_B_1 = 0.d0
  BioTH_Create%monod_D = 0.d0
  BioTH_Create%inhibition_B = 0.d0
  BioTH_Create%background_concentration_B = 0.d0
  BioTH_Create%mass_action_CD = 0.d0
  BioTH_Create%stoichiometric_C = 1.d0
  BioTH_Create%stoichiometric_D_2 = 1.d0
  BioTH_Create%rate_D_i = 0.d0
  BioTH_Create%rate_D_m = 0.d0
  BioTH_Create%exponent_B = 0.d0
  BioTH_Create%inhibition_I = 0.d0
  BioTH_Create%mass_action_B = 0.d0
  BioTH_Create%mass_action_X = 0.d0
  BioTH_Create%density_B = 0.d0

  nullify(BioTH_Create%next)

end function BioTH_Create

! ************************************************************************** !
subroutine BioTH_Read(this,input,option)
  !
  ! Reads input deck for reaction sandbox parameters
  !
  ! Author: Edwin Saavedra C.
  ! Date: 10/01/2020 - created
  !       05/18/2020 - fixed formatting
  !

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_bioth_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
    call StringToUpper(word)

    select case(trim(word))
	  case('NAME_B')
        call InputReadWord(input,option,this%name_B,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_B', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('NAME_C')
        call InputReadWord(input,option,this%name_C,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_C', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('NAME_D_MOBILE')
        call InputReadWord(input,option,this%name_D_mobile,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_D_mobile', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('NAME_D_IMMOBILE')
        call InputReadWord(input,option,this%name_D_immobile,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_D_immobile', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('NAME_I')
        call InputReadWord(input,option,this%name_I,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_I', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('NAME_X')
        call InputReadWord(input,option,this%name_X,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_X', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('NAME_BIOMINERAL')
        call InputReadWord(input,option,this%name_biomineral,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_biomineral', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('BETA') 
        call InputReadDouble(input,option,this%beta) 
        call InputErrorMsg(input,option,'beta', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')     
      case('ALPHA_VEL') 
        call InputReadDouble(input,option,this%alpha_vel) 
        call InputErrorMsg(input,option,'alpha_vel', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')    

      case('BETA_VEL') 
        call InputReadDouble(input,option,this%beta_vel) 
        call InputErrorMsg(input,option,'beta_vel', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')						   
      
      case('K')
        call InputReadDouble(input, option, this%k)
        call InputErrorMsg(input, option, 'k', 'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('K_O')
        call InputReadDouble(input, option, this%K_O)
        call InputErrorMsg(input, option, 'K_O', 'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('NAME_O2')
        call InputReadWord(input,option,this%name_O2,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_O2', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')  
	  case('NAME_CO2')
        call InputReadWord(input,option,this%name_CO2,PETSC_TRUE)
        call InputErrorMsg(input,option,'name_CO2', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')  
	  
      case('EXPONENT_B')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%exponent_B)
        call InputErrorMsg(input,option,'exponent_B', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('BACKGROUND_CONC_B')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%background_concentration_B )
        call InputErrorMsg(input,option,'background_concentration_B', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('MASS_ACTION_B')
        call InputReadDouble(input,option,this%mass_action_B)
        call InputErrorMsg(input,option,'mass_action_B', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('MASS_ACTION_CD')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%mass_action_CD)
        call InputErrorMsg(input,option,'mass_action_CD', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('MASS_ACTION_X')
        call InputReadDouble(input,option,this%mass_action_X)
        call InputErrorMsg(input,option,'mass_action_X', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('INHIBITION_B')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%inhibition_B)
        call InputErrorMsg(input,option,'inhibition_B', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('INHIBITION_C')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%inhibition_C)
        call InputErrorMsg(input,option,'inhibition_C', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('MONOD_D')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%monod_D)
        call InputErrorMsg(input,option,'monod_D', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('INHIBITION_I')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%inhibition_I)
        call InputErrorMsg(input,option,'inhibition_I', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('RATE_B_1')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_B_1)
        call InputErrorMsg(input,option,'rate_B_1', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('RATE_B_2')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_B_2)
        call InputErrorMsg(input,option,'rate_B_2', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('RATE_C')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_C)
        call InputErrorMsg(input,option,'rate_C', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('RATE_D')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_D)
        call InputErrorMsg(input,option,'rate_D', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('RATE_D_IMMOB')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%rate_D_i)
        call InputErrorMsg(input,option,'rate_D_i', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('RATE_D_MOBIL')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%rate_D_m)
        call InputErrorMsg(input,option,'rate_D_m', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('DENSITY_B')
        call InputReadDouble(input,option,this%density_B)
        call InputErrorMsg(input,option,'density_B', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('STOICHIOMETRIC_C')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%stoichiometric_C )
        call InputErrorMsg(input,option,'stoichiometric_C', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('STOICHIOMETRIC_D_1')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%stoichiometric_D_1)
        call InputErrorMsg(input,option,'stoichiometric_D_1', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
      case('STOICHIOMETRIC_D_2')
        ! Read the double precision background concentration in kg/m^3
        call InputReadDouble(input,option,this%stoichiometric_D_2)
        call InputErrorMsg(input,option,'stoichiometric_D_2', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('PARTICLE_NAME_AQ')
        ! Bioparticle name while in suspension
        call InputReadWord(input,option,this%name_aqueous,PETSC_TRUE)
        call InputErrorMsg(input,option,'PARTICLE_NAME_AQ', &
               'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('PARTICLE_NAME_IM')
        ! Bioparticle name while immobilized
        call InputReadWord(input,option,this%name_immobile,PETSC_TRUE)
        call InputErrorMsg(input,option,'PARTICLE_NAME_IM', &
               'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')

      case('DECAY_AQUEOUS')
        ! Decay rate while in the aqueous phase
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'DECAY_AQUEOUS', &
               'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
        select case(trim(word))
          case('CONSTANT')
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
                case('VALUE')
                ! Read the double precision rate constant
                  call InputReadDouble(input,option,this%decay_aqueous)
                  call InputErrorMsg(input,option,'VALUE', &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_AQUEOUS,CONSTANT')
                  ! Read the units
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (InputError(input)) then
                  ! If units do not exist, assume default units of 1/s which are the
                  ! standard internal PFLOTRAN units for this rate constant.
                    input%err_buf = 'REACTION_SANDBOX,BIOPARTICLE,&
                                     &DECAY_AQUEOUS,RATE CONSTANT UNITS'
                    call InputDefaultMsg(input,option)
                  else
                    ! If units exist, convert to internal units of 1/s
                    internal_units = 'unitless/sec'
                    this%decay_aqueous = this%decay_aqueous * &
                      UnitsConvertToInternal(word,internal_units, &
                                  'REACTION_SANDBOX,BIOPARTICLE,&
                                     &DECAY_AQUEOUS,RATE CONSTANT UNITS', &
                                   option)
                  endif
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_AQUEOUS,CONSTANT',option)
              end select
            end do
            call InputPopBlock(input,option)

          case('TEMPERATURE_MODEL')
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'TEMPERATURE_MODEL', &
                     'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_AQUEOUS')
              call StringToUpper(word)

              select case(trim(word))
                case('TREF')
                ! Reference temperature (Probably 4°C)
                  call InputReadDouble(input,option,this%Tref_aqueous)
                  call InputErrorMsg(input,option,'TREF',&
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_AQUEOUS,TEMPERATURE_MODEL')
                case('ZT')
                ! Model parameter zT
                  call InputReadDouble(input,option,this%zT_aqueous)
                  call InputErrorMsg(input,option,'ZT', &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_AQUEOUS,TEMPERATURE_MODEL')
                case('N')
                  ! Model parameter n (Probably 1.0 or 2.0 )
                  call InputReadDouble(input,option,this%nAq_aqueous)
                  call InputErrorMsg(input,option,'N', &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_AQUEOUS,TEMPERATURE_MODEL')
                case('LOGDREF')
                ! D reference value (Probably 2.3)
                  call InputReadDouble(input,option,this%logDref_aqueous)
                  call InputErrorMsg(input,option,'LOGDREF', &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_AQUEOUS,TEMPERATURE_MODEL')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_AQUEOUS,TEMPERATURE_MODEL',option)
              end select
            end do
            call InputPopBlock(input,option)

          case default
            call InputKeywordUnrecognized(input,word, &
                    'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_AQUEOUS', &
                    option)
        end select

      case('DECAY_ADSORBED')
        ! Decay rate while in the immobile phase
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'DECAY_ADSORBED', &
               'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
        select case(trim(word))
          case('CONSTANT')
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'CONSTANT', &
                     'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_ADSORBED')
              call StringToUpper(word)

              select case(trim(word))
                case('VALUE')
                ! Read the double precision rate constant
                  call InputReadDouble(input,option,this%decay_adsorbed)
                  call InputErrorMsg(input,option,'VALUE', &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_ADSORBED,CONSTANT')
                  ! Read the units
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (InputError(input)) then
                  ! If units do not exist, assume default units of 1/s which are the
                  ! standard internal PFLOTRAN units for this rate constant.
                    input%err_buf = 'REACTION_SANDBOX,BIOPARTICLE,&
                                    &DECAY_ADSORBED,CONSTANT,VALUE,&
                                    &RATE CONSTANT UNITS'
                    call InputDefaultMsg(input,option)
                  else
                    ! If units exist, convert to internal units of 1/s
                    internal_units = 'unitless/sec'
                    this%decay_adsorbed = this%decay_adsorbed * &
                      UnitsConvertToInternal(word,internal_units, &
                            'REACTION_SANDBOX,BIOPARTICLE,DECAY_ADSORBED,&
                            &CONSTANT,VALUE,RATE CONSTANT UNITS',option)
                  endif
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_ADSORBED,CONSTANT',option)
              end select
            end do
            call InputPopBlock(input,option)

          case('TEMPERATURE_MODEL')
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'TEMPERATURE_MODEL', &
                     'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_ADSORBED')
              call StringToUpper(word)

              select case(trim(word))
                case('TREF')
                  call InputReadDouble(input,option,this%Tref_adsorbed)
                  call InputErrorMsg(input,option,'TREF',&
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_ADSORBED,TEMPERATURE_MODEL')
                case('ZT')
                  call InputReadDouble(input,option,this%zT_adsorbed)
                  call InputErrorMsg(input,option,'ZT', &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_ADSORBED,TEMPERATURE_MODEL')
                case('N')
                  call InputReadDouble(input,option,this%nAq_adsorbed)
                  call InputErrorMsg(input,option,'N', &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_ADSORBED,TEMPERATURE_MODEL')

                case('LOGDREF')
                  call InputReadDouble(input,option,this%logDref_adsorbed)
                  call InputErrorMsg(input,option,'LOGDREF', &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_ADSORBED,TEMPERATURE_MODEL')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_ADSORBED,TEMPERATURE_MODEL',option)
              end select
            end do
            call InputPopBlock(input,option)
          case default
            call InputKeywordUnrecognized(input,word, &
                   'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_ADSORBED', &
                   option)
        end select

      case('RATE_ATTACHMENT')
        ! Decay rate while in the aqueous phase
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'RATE_ATTACHMENT', &
               'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
        select case(trim(word))
          case('CONSTANT')
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'CONSTANT', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE&
                           &,RATE_ATTACHMENT')
              call StringToUpper(word)

              select case(trim(word))
                case('VALUE')
                ! Read the double precision rate constant
                  call InputReadDouble(input,option,this%rate_attachment)
                  call InputErrorMsg(input,option,'VALUE', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE&
                           &,RATE_ATTACHMENT,CONSTANT')
                  ! Read the units
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (InputError(input)) then
                  ! If units do not exist, assume default units of 1/s which are the
                  ! standard internal PFLOTRAN units for this rate constant.
                    input%err_buf = 'RATE CONSTANT UNITS assumed as 1/s'
                    call InputDefaultMsg(input,option)
                  else
                    ! If units exist, convert to internal units of 1/s
                    internal_units = 'unitless/sec'
                    this%rate_attachment = this%rate_attachment * &
                      UnitsConvertToInternal(word,internal_units, &
                                    'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE&
                                    &,RATE_ATTACHMENT,CONSTANT',option)
                  endif
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &RATE_ATTACHMENT',option)
              end select
            end do
            call InputPopBlock(input,option)

          case('FILTRATION_MODEL')
            call InputPushBlock(input,option)
            ! Gotta turn this on so Darcy velocity is stored in global
            ! and can be used in the reaction sandbox
            option%flow%store_darcy_vel = PETSC_TRUE
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'FILTRATION_MODEL', &
                     'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
              call StringToUpper(word)

              select case(trim(word))
                case('DIAMETER_COLLECTOR')
                ! Diameter of the collector, i.e., soil grain size [m]
                  call InputReadDouble(input,option,this%diam_collector)
                  call InputErrorMsg(input,option,'DIAMETER_COLLECTOR', &
                        'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE, &
                         &FILTRATION_MODEL')

                case('DIAMETER_PARTICLE')
                ! Diameter of the bioparticle [m]
                  call InputReadDouble(input,option,this%diam_particle)
                  call InputErrorMsg(input,option,'DIAMETER_PARTICLE', &
                        'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE, &
                         &FILTRATION_MODEL')

                case('HAMAKER_CONSTANT')
                ! Hamaker constant particle-soil pair [Joules]
                  call InputReadDouble(input,option,this%hamaker_constant)
                  call InputErrorMsg(input,option,'HAMAKER_CONSTANT', &
                        'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE, &
                         &FILTRATION_MODEL')

                case('DENSITY_PARTICLE')
                ! Density of the bioparticles [kg/m3]
                  call InputReadDouble(input,option,this%density_particle)
                  call InputErrorMsg(input,option,'DENSITY_PARTICLE', &
                        'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE, &
                         &FILTRATION_MODEL')

                case('ALPHA_EFFICIENCY')
                ! Collision/attachment efficiency [-]
                  call InputReadDouble(input,option,this%alpha_efficiency)
                  call InputErrorMsg(input,option,'ALPHA_EFFICIENCY', &
                        'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE, &
                         &FILTRATION_MODEL')

                case('DEBUG')
                ! Print stuff on screen
                  print *, "Will debug -> print CFT stuff: "
                  this%debug_option = .True. ! Edwin debugging

                ! Something else
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &DECAY_ADSORBED,FILTRATION_MODEL',option)
              end select
            end do
            call InputPopBlock(input,option)

          case default
            call InputKeywordUnrecognized(input,word, &
                   'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,DECAY_ADSORBED', &
                   option)
        end select

      ! Detachment rate
      case('RATE_DETACHMENT')
        ! Decay rate while in the aqueous phase
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'RATE_DETACHMENT', &
               'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE')
        select case(trim(word))
          case('CONSTANT')
            call InputPushBlock(input,option)
            do
              call InputReadPflotranString(input,option)
              if (InputError(input)) exit
              if (InputCheckExit(input,option)) exit

              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'CONSTANT', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE&
                           &,RATE_DETACHMENT')
              call StringToUpper(word)

              select case(trim(word))
                case('VALUE')
                ! Read the double precision rate constant
                  call InputReadDouble(input,option,this%rate_detachment)
                  call InputErrorMsg(input,option,'VALUE', &
                           'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                            &RATE_DETACHMENT,CONSTANT')
                  ! Read the units
                  call InputReadWord(input,option,word,PETSC_TRUE)
                  if (InputError(input)) then
                  ! If units do not exist, assume default units of 1/s which are the
                  ! standard internal PFLOTRAN units for this rate constant.
                    input%err_buf = 'RATE CONSTANT UNITS assumed as 1/s'
                    call InputDefaultMsg(input,option)
                  else
                    ! If units exist, convert to internal units of 1/s
                    internal_units = 'unitless/sec'
                    this%rate_detachment = this%rate_detachment * &
                      UnitsConvertToInternal(word,internal_units, &
                                  'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                                  &RATE_DETACHMENT,CONSTANT',option)
                  endif
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE,&
                         &RATE_DETACHMENT',option)
              end select
            end do
            call InputPopBlock(input,option)

          case default
            call InputKeywordUnrecognized(input,word, &
                    'CHEMISTRY,REACTION_SANDBOX,BIOPARTICLE',option)
        end select
      case default
        call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,',option)
    end select
  end do

  call InputPopBlock(input,option)

end subroutine BioTH_Read

! ************************************************************************** !
subroutine BioTH_Setup(this,reaction,option)
  !
  ! Sets up the kinetic attachment/dettachment reactions
  !
  ! Author: Edwin Saavedra C.
  ! Date: 10/01/2020
  !

  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_bioth_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  this%species_Vaq_id = &
    GetPrimarySpeciesIDFromName(this%name_aqueous,reaction,option)

  this%species_Vim_id = &
    GetImmobileSpeciesIDFromName(this%name_immobile,reaction%immobile,option)
  this%D_mobile_id = &
    GetPrimarySpeciesIDFromName(this%name_D_mobile, &
                                reaction,option)

  this%C_id = &
    GetPrimarySpeciesIDFromName(this%name_C, &
                                reaction,option)

  this%I_id = &
    GetPrimarySpeciesIDFromName(this%name_I, &
                                reaction,option)

  this%X_id = &
    GetPrimarySpeciesIDFromName(this%name_X, &
                                reaction,option)
  this%O2_id = &
    GetPrimarySpeciesIDFromName(this%name_O2, &
                                reaction,option)
  this%CO2_id = &
    GetPrimarySpeciesIDFromName(this%name_CO2, &
                                reaction,option)                               
  this%D_immobile_id = &
    GetImmobileSpeciesIDFromName(this%name_D_immobile, &
                                 reaction%immobile,option)
  this%biomineral_id = &
    GetMineralIDFromName(this%name_biomineral, &
                         reaction%mineral,option)

end subroutine BioTH_Setup

! ************************************************************************** !
subroutine BioTH_React(this,Residual,Jacobian,compute_derivative, &
                        rt_auxvar,global_auxvar,material_auxvar, &
                        reaction, option)
  !
  ! Evaluates reaction
  !
  ! Author: Edwin
  ! Date: 04/09/2020 - created
  !       05/18/2021 - remove truncate concentrations
  !

  use Option_module
  use String_module
  use Reaction_Aux_module, only : reaction_rt_type
  use Reaction_Immobile_Aux_module
  use Material_Aux_module, only : material_auxvar_type

  implicit none

  class(reaction_sandbox_bioth_type) :: this
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
  PetscReal :: volume                 ! m^3 bulk
  PetscReal :: porosity               ! m^3 pore space / m^3 bulk
  PetscReal :: liquid_saturation      ! m^3 water / m^3 pore space
  PetscReal :: L_water                ! L water
  PetscReal :: temperature
  PetscReal :: rho_f
  PetscReal :: g
  PetscReal :: viscosity


  PetscReal :: Vaq  ! mol/L
  PetscReal :: Vim  ! mol/m^3
  PetscReal :: Rate
  PetscReal :: RateAtt, RateDet  ! mol/sec
  PetscReal :: RateDecayAq, RateDecayIm !Check units
  PetscReal :: stoichVaq
  PetscReal :: stoichVim

  ! Decay model parameters
  PetscReal :: decayAq, decayIm
  PetscReal :: logDrefAq, logDrefAd
  PetscReal :: TrefAq, TrefAd
  PetscReal :: zTAq, zTAd
  PetscReal :: nAq, nAd

  ! Filtration model parameters for attachment
  PetscReal :: katt
  PetscReal :: dc,dp,Hamaker,rho_p, alpha, qMag, diffusionCoeff, kB
  PetscReal :: Gm, Gm5, Happel
  PetscReal :: NR, NPe, NvdW, Ngr
  PetscReal :: Eta_D, Eta_I, Eta_G, Eta_0

  ! Detachment rate
  PetscReal :: kdet
  PetscReal :: mu_B_im, mu_B_im_residual
  PetscReal :: mu_B_mob, mu_B_mob_residual
  PetscReal :: sum_food
  PetscReal :: immobile_mole_fraction, mobile_mole_fraction

  ! Global stuff (Check global_aux.F90)
  volume = material_auxvar%volume
  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  L_water = porosity*liquid_saturation*volume*1.d3
  ! 1.d3 converts m^3 water -> L water

  viscosity = 0.0008891
  ! Ns/m2 (This should be a field but not found in global_auxvar)

  temperature = global_auxvar%temp
  rho_f = global_auxvar%den_kg(iphase)
  g = EARTH_GRAVITY


! Assign concentrations of Vaq and Vim
  Vaq = rt_auxvar%total(this%species_Vaq_id,iphase)
  Vim = rt_auxvar%immobile(this%species_Vim_id)

  ! initialize all rates to zero
  Rate = 0.d0
  RateAtt = 0.d0
  RateDet = 0.d0
  RateDecayAq = 0.d0
  RateDecayIm = 0.d0

  ! stoichiometries
  ! reactants have negative stoichiometry
  ! products have positive stoichiometry
  stoichVaq = -1.d0
  stoichVim = 1.d0

  ! kinetic rate constants
  katt = 0.d0
  kdet = 0.d0
  kdet = this%rate_detachment

  !!!!!!!!!!!!!!!!!!!
  ! Decay rate - Aqueous phase
  !
  !  Check Guillier et al. for model equation (2020)
  !  decay = ln(10)/D
  !  logD = logDref - [(T-Tref)/zT]^n
  !
  !!!!!!!!!!!!!!!!!!!
  logDrefAq = 0.0
  TrefAq = 0.0
  zTAq = 0.0
  nAq = 0.0

  IF (this%decay_aqueous < 0.d0) THEN
    logDrefAq = this%logDref_aqueous
    TrefAq = this%Tref_aqueous
    zTAq = this%zT_aqueous
    nAq = this%nAq_aqueous

    decayAq = (2.302585/(10.0 ** (logDrefAq - (((temperature - TrefAq)/zTAq)**nAq))))/3600  ! 1/s
  ELSE
    decayAq = this%decay_aqueous
  END IF


!!!!!!!!!!!!!!!!!!!
! Decay rate - Immobile phase
!!!!!!!!!!!!!!!!!!!

  logDrefAd = 0.0
  TrefAd = 0.0
  zTAd = 0.0
  nAd = 0.0

  IF (this%decay_adsorbed < 0.d0) THEN
    logDrefAd = this%logDref_adsorbed
    TrefAd = this%Tref_adsorbed
    zTAd = this%zT_adsorbed
    nAd = this%nAq_adsorbed
    decayIm = (2.302585/(10.0 ** (logDrefAd - (((temperature - TrefAd)/zTAd)**nAd))))/3600  ! 1/s

  ELSE
    decayIm = this%decay_adsorbed
  END IF

!!!!!!!!!!!!!!!!!!!
! Attachment rate
!!!!!!!!!!!!!!!!!!!
  dc = 0.0
  dp = 0.0
  Hamaker = 0.0
  rho_p = 0.0
  alpha = 0.0

  IF (this%rate_attachment < 0.d0) THEN
    kB = this%BOLTZMANN_CONSTANT
    dc = this%diam_collector
    dp = this%diam_particle
    rho_p = this%density_particle
    alpha = this%alpha_efficiency
    Hamaker = this%hamaker_constant

    qMag = MAX(global_auxvar%darcy_vel(iphase),1.0d-20)

    diffusionCoeff = this%BOLTZMANN_CONSTANT*(temperature+273.15) / &
                   (3.0 * PI * viscosity * dp)

    ! Non-dimensional parameters
    !! Happel parameter As
    Gm = (1.0 - porosity)**(1./3.)
    Gm5 = Gm*Gm*Gm*Gm*Gm
    Happel = (2.0 * (1.0 - Gm5)) / (2.0 - (3.0*Gm) + (3.0*Gm5) - (2.0*Gm*Gm5))

    !! Aspect ratio
    NR = dp/dc

    !! Péclet number
    NPe = (qMag * dc) / (diffusionCoeff)

    !! van der Waals number
    NvdW = Hamaker/(kB*(temperature + 273.15))

    !! Gravitational number
    NGr = PI/12.0 * (dp*dp*dp*dp) * (rho_p - rho_f) * g /&
          (kB*(temperature + 273.15))

    ! Collector efficiencies
    ! ( see Tufenkji & Elimelech 2014
    !   DOI : 10.1021/es034049r )
    !! Transport by diffusion
    Eta_D = 2.4 &
            * Happel**(1./3.) &
            * NR**(-0.081) &
            * NPe**(-0.715) &
            * NvdW**(0.052)

    !! Transport by interception
    Eta_I = 0.55 &
            * Happel &
            * NR**(1.55) &
            * NPe**(-0.125) &
            * NvdW**(0.125)

    !! Transport due to gravity
    Eta_G = 0.475 &
            * NR**(-1.35) &
            * NPe**(-1.11) &
            * NvdW**(0.053) &
            * NGr**(1.11)

    !! Single collector efficiency
    Eta_0 = Eta_D + Eta_I + Eta_G

    ! Rate of attachment according to CFT
    katt = 1.5 * (1 - porosity) * qMag * alpha * Eta_0 &
           / (dc * porosity)

    ! Edwin debugging
    if(this%debug_option) then
      print '(3x,"porosity = ", ES12.4)', porosity
      print '(3x,"temp C   = ", ES12.4)', Temperature
      print '(3x,"viscosit = ", ES12.4)', viscosity
      print '(3x,"densityF = ", ES12.4)', rho_f
      print '(3x,"densityP = ", ES12.4)', rho_p

      print '(3x,"DarcyqMa = ", ES12.4)', qMag
      print '(3x,"diffCoef = ", ES12.4)', diffusionCoeff
      print '(3x,"Gm       = ", ES12.4)', Gm
      print '(3x,"Gm5      = ", ES12.4)', Gm5
      print '(3x,"Happel   = ", ES12.4)', Happel
      print '(3x,"NR       = ", ES12.4)', NR
      print '(3x,"NPe      = ", ES12.4)', NPe
      print '(3x,"NvW      = ", ES12.4)', NvdW
      print '(3x,"NGr      = ", ES12.4)', NGr
      print '(3x,"EtaD     = ", ES12.4)', Eta_D
      print '(3x,"EtaI     = ", ES12.4)', Eta_I
      print '(3x,"EtaG     = ", ES12.4)', Eta_G
      print '(3x,"Eta0     = ", ES12.4)', Eta_0
      print '(3x,"katt     = ", ES12.4)', katt
      print *, "--------------------"
    endif

  else
    ! A constant rate of attachment
    katt = this%rate_attachment
  end if

!!!!!!!!!!!!!!!!!!!
! Detachment rate
!!!!!!!!!!!!!!!!!!!
  if (this%rate_detachment < 0.d0) then
    kdet = 0.0
  else
    ! A constant rate of attachment
    kdet = this%rate_detachment
  end if

  RateAtt = 0.0
  RateDet = 0.0
  RateDecayAq = 0.0
  RateDecayIm = 0.0


  ! Build here for attachment/detachment
  ! first-order forward - reverse (A <-> C)
  Rate = katt * Vaq * L_water - kdet * Vim * volume
  RateAtt = stoichVaq * Rate
  RateDet = stoichVim * Rate

  ! Build here for inactivation reactions
  ! first-order (A -> X)
  Rate = decayAq * Vaq * L_water
  RateDecayAq = - Rate

  Rate = decayIm * Vim * volume
  RateDecayIm = - Rate
  !immobile_to_water_vol = material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0
  sum_food = rt_auxvar%total(this%D_mobile_id,iphase)*L_water + &
             rt_auxvar%immobile(this%D_immobile_id)*volume               ! in mol

  mu_B_im = this%rate_B_1*Vim* &      ! mol/Ls
			(sum_food/(sum_food + this%monod_D))* &
			! !(rt_auxvar%total(idof_O2,iphase) / &        !oxygen 
			! !(this%K_O + rt_auxvar%total(idof_O2,iphase)))*&             ! limitation
            (this%inhibition_B/ (Vim + this%inhibition_B))**this%exponent_B 
  mu_B_mob = this%rate_B_1*.5*Vaq* &      ! mol/Ls
			(sum_food/(sum_food + this%monod_D))* &
			! !(rt_auxvar%total(idof_O2,iphase) / &        !oxygen 
			! !(this%K_O + rt_auxvar%total(idof_O2,iphase)))*&             ! limitation
            (this%inhibition_B/ (Vaq + this%inhibition_B))**this%exponent_B
			
  qMag = MAX(global_auxvar%darcy_vel(iphase),1.0d-20)
  
  mu_B_mob_residual = -1* mu_B_mob*L_water + & 
  ((this%alpha_vel*qMag)**this%beta_vel)* & 
  this%rate_B_2*(Vaq - this%background_concentration_B)* L_water  

  mu_B_im_residual = -1* mu_B_im*volume + & 
  ((this%alpha_vel*qMag)**this%beta_vel)* & 
  this%rate_B_2*(Vim - this%background_concentration_B)* volume  


  mobile_mole_fraction = rt_auxvar%total(this%D_mobile_id,iphase)*L_water/sum_food
  immobile_mole_fraction = 1 - mobile_mole_fraction


  Residual(this%D_mobile_id) = Residual(this%D_mobile_id) + &
                           ! ! Growth usage, mol/s
                           this%stoichiometric_D_1* &
                           mobile_mole_fraction* &                                ! dimensionless
                           mu_B_mob*L_water + &                        ! mol/m3 bulk/s * m3 bulk
                           ! ! Direct usage, mol/s
                           this%rate_D* &                          ! 1/s
                           mobile_mole_fraction* &                                ! dimensionless
                           Vaq*L_water + &                        ! mol/L/s * m3 bulk
                           ! ! immobilization, mol/s
                           this%rate_D_i* &                   ! 1/s
                           rt_auxvar%total(this%D_mobile_id,iphase)* &            ! mol/L
                           L_water - &   ! L water/m3 bulk
                           ! ! remobilization, mol/s
                           this%rate_D_m* &                   ! 1/s
                           rt_auxvar%immobile(this%D_immobile_id)* &           ! mol/m3 bulk
                           volume                                 ! m3 bulk


  Residual(this%D_immobile_id) = Residual(this%D_immobile_id) + &
                           ! ! Growth usage, mol/s
                           this%stoichiometric_D_1* &                          ! unitless
                           immobile_mole_fraction* &                              ! dimensionless
                           mu_B_im*volume + &                        ! mol/m3 bulk/s * m3 bulk
                           ! ! Direct usage, mol/s
                           this%rate_D* &                          ! 1/s
                           immobile_mole_fraction* &                              ! dimensionless
                           Vim* volume - &                              ! L water/m3 bulk
                           ! ! immobilization, mol/s
                           this%rate_D_i* &                   ! 1/s
                           rt_auxvar%total(this%D_mobile_id,iphase)* &            ! mol/L
                           L_water + & ! remobilization, mol/s
                           this%rate_D_m* &                   ! 1/s
                           rt_auxvar%immobile(this%D_immobile_id)* &           ! mol/m3 bulk
                           volume
! This awful block just tries to
! avoid concentrations below 1E-50
! (Is this avoided with TRUNCATE_CONCENTRATION ?) > sure it does
  ! if ( Vaq > 0.0 ) then
  !   if ( Vim > 0.0 ) then
  !     !Do nothing
  !   else if ( Vim <= 0.0 ) then
  !     Vim = 1.0d-50
  !     RateDet = 0.0
  !     RateDecayIm = 0.0
  !   end if
  ! else if ( Vaq <= 0.0 ) then
  !   if ( Vim > 0.0 ) then
  !     Vaq = 1.0d-50
  !     RateAtt = 0.0
  !     RateDecayAq = 0.0
  !   else if ( Vim <= 0.0 ) then
  !     Vim = 1.0d-50
  !     Vaq = 1.0d-50
  !     RateAtt = 0.0
  !     RateDet = 0.0
  !     RateDecayAq = 0.0
  !     RateDecayIm = 0.0
  !   end if
  ! end if

  ! The actual calculation:

  Residual(this%species_Vaq_id) = &
    Residual(this%species_Vaq_id) - RateAtt - RateDecayAq + RateDet + mu_B_mob_residual

  Residual(this%species_Vim_id + reaction%offset_immobile) = &
    Residual(this%species_Vim_id + reaction%offset_immobile) &
    - RateDet - RateDecayIm + mu_B_im_residual

  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second
  ! Residual(this%species_Vaq_id) = &
  !   Residual(this%species_Vaq_id) - RateAtt - RateDecayAq

  ! Residual(this%species_Vim_id + reaction%offset_immobile) = &
  !   Residual(this%species_Vim_id + reaction%offset_immobile) &
  !   - RateDet - RateDecayIm

end subroutine BioTH_React

subroutine BioTH_KineticState(this,rt_auxvar,global_auxvar, &
                                               material_auxvar, &
                                               reaction,option)
  !
  ! Updates the kinetic state for the sandbox
  !
  use Option_module
  use String_module
  use Reaction_Aux_module, only : reaction_rt_type
  use Reaction_Immobile_Aux_module
  use Material_Aux_module, only : material_auxvar_type

  implicit none

  class(reaction_sandbox_bioth_type) :: this
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
  PetscReal :: volume                 ! m^3 bulk
  PetscReal :: porosity               ! m^3 pore space / m^3 bulk
  PetscReal :: liquid_saturation      ! m^3 water / m^3 pore space
  PetscReal :: L_water                ! L water
  PetscReal :: temperature
  PetscReal :: rho_f
  PetscReal :: g
  PetscReal :: viscosity
  PetscReal :: delta_volfrac
  PetscReal :: Vaq  ! mol/L
  PetscReal :: Vim  ! mol/m^3
  PetscReal :: Rate
  PetscReal :: RateAtt, RateDet  ! mol/sec
  PetscReal :: RateDecayAq, RateDecayIm !Check units
  PetscReal :: stoichVaq
  PetscReal :: stoichVim

  ! Decay model parameters
  PetscReal :: decayAq, decayIm
  PetscReal :: logDrefAq, logDrefAd
  PetscReal :: TrefAq, TrefAd
  PetscReal :: zTAq, zTAd
  PetscReal :: nAq, nAd

  ! Filtration model parameters for attachment
  PetscReal :: katt
  PetscReal :: dc,dp,Hamaker,rho_p, alpha, qMag, diffusionCoeff, kB
  PetscReal :: Gm, Gm5, Happel
  PetscReal :: NR, NPe, NvdW, Ngr
  PetscReal :: Eta_D, Eta_I, Eta_G, Eta_0

  ! Detachment rate
  PetscReal :: kdet
  PetscReal :: mu_B_im, mu_B_im_residual
  PetscReal :: mu_B_mob, mu_B_mob_residual
  PetscReal :: sum_food
  PetscReal :: immobile_mole_fraction, mobile_mole_fraction

  ! Global stuff (Check global_aux.F90)
  volume = material_auxvar%volume
  porosity = material_auxvar%porosity
  liquid_saturation = global_auxvar%sat(iphase)
  L_water = porosity*liquid_saturation*volume*1.d3
  ! 1.d3 converts m^3 water -> L water

  viscosity = 0.0008891
  ! Ns/m2 (This should be a field but not found in global_auxvar)

  temperature = global_auxvar%temp
  rho_f = global_auxvar%den_kg(iphase)
  g = EARTH_GRAVITY


! Assign concentrations of Vaq and Vim
  Vaq = rt_auxvar%total(this%species_Vaq_id,iphase)
  Vim = rt_auxvar%immobile(this%species_Vim_id)

  ! initialize all rates to zero
  Rate = 0.d0
  RateAtt = 0.d0
  RateDet = 0.d0
  RateDecayAq = 0.d0
  RateDecayIm = 0.d0

  ! stoichiometries
  ! reactants have negative stoichiometry
  ! products have positive stoichiometry
  stoichVaq = -1.d0
  stoichVim = 1.d0

  ! kinetic rate constants
  katt = 0.d0
  kdet = 0.d0
  kdet = this%rate_detachment

  !!!!!!!!!!!!!!!!!!!
  ! Decay rate - Aqueous phase
  !
  !  Check Guillier et al. for model equation (2020)
  !  decay = ln(10)/D
  !  logD = logDref - [(T-Tref)/zT]^n
  !
  !!!!!!!!!!!!!!!!!!!
  logDrefAq = 0.0
  TrefAq = 0.0
  zTAq = 0.0
  nAq = 0.0

  IF (this%decay_aqueous < 0.d0) THEN
    logDrefAq = this%logDref_aqueous
    TrefAq = this%Tref_aqueous
    zTAq = this%zT_aqueous
    nAq = this%nAq_aqueous

    decayAq = (2.302585/(10.0 ** (logDrefAq - (((temperature - TrefAq)/zTAq)**nAq))))/3600  ! 1/s
  ELSE
    decayAq = this%decay_aqueous
  END IF


!!!!!!!!!!!!!!!!!!!
! Decay rate - Immobile phase
!!!!!!!!!!!!!!!!!!!

  logDrefAd = 0.0
  TrefAd = 0.0
  zTAd = 0.0
  nAd = 0.0

  IF (this%decay_adsorbed < 0.d0) THEN
    logDrefAd = this%logDref_adsorbed
    TrefAd = this%Tref_adsorbed
    zTAd = this%zT_adsorbed
    nAd = this%nAq_adsorbed
    decayIm = (2.302585/(10.0 ** (logDrefAd - (((temperature - TrefAd)/zTAd)**nAd))))/3600  ! 1/s

  ELSE
    decayIm = this%decay_adsorbed
  END IF

!!!!!!!!!!!!!!!!!!!
! Attachment rate
!!!!!!!!!!!!!!!!!!!
  dc = 0.0
  dp = 0.0
  Hamaker = 0.0
  rho_p = 0.0
  alpha = 0.0

  IF (this%rate_attachment < 0.d0) THEN
    kB = this%BOLTZMANN_CONSTANT
    dc = this%diam_collector
    dp = this%diam_particle
    rho_p = this%density_particle
    alpha = this%alpha_efficiency
    Hamaker = this%hamaker_constant

    qMag = MAX(global_auxvar%darcy_vel(iphase),1.0d-20)

    diffusionCoeff = this%BOLTZMANN_CONSTANT*(temperature+273.15) / &
                   (3.0 * PI * viscosity * dp)

    ! Non-dimensional parameters
    !! Happel parameter As
    Gm = (1.0 - porosity)**(1./3.)
    Gm5 = Gm*Gm*Gm*Gm*Gm
    Happel = (2.0 * (1.0 - Gm5)) / (2.0 - (3.0*Gm) + (3.0*Gm5) - (2.0*Gm*Gm5))

    !! Aspect ratio
    NR = dp/dc

    !! Péclet number
    NPe = (qMag * dc) / (diffusionCoeff)

    !! van der Waals number
    NvdW = Hamaker/(kB*(temperature + 273.15))

    !! Gravitational number
    NGr = PI/12.0 * (dp*dp*dp*dp) * (rho_p - rho_f) * g /&
          (kB*(temperature + 273.15))

    ! Collector efficiencies
    ! ( see Tufenkji & Elimelech 2014
    !   DOI : 10.1021/es034049r )
    !! Transport by diffusion
    Eta_D = 2.4 &
            * Happel**(1./3.) &
            * NR**(-0.081) &
            * NPe**(-0.715) &
            * NvdW**(0.052)

    !! Transport by interception
    Eta_I = 0.55 &
            * Happel &
            * NR**(1.55) &
            * NPe**(-0.125) &
            * NvdW**(0.125)

    !! Transport due to gravity
    Eta_G = 0.475 &
            * NR**(-1.35) &
            * NPe**(-1.11) &
            * NvdW**(0.053) &
            * NGr**(1.11)

    !! Single collector efficiency
    Eta_0 = Eta_D + Eta_I + Eta_G

    ! Rate of attachment according to CFT
    katt = 1.5 * (1 - porosity) * qMag * alpha * Eta_0 &
           / (dc * porosity)

    ! Edwin debugging
    if(this%debug_option) then
      print '(3x,"porosity = ", ES12.4)', porosity
      print '(3x,"temp C   = ", ES12.4)', Temperature
      print '(3x,"viscosit = ", ES12.4)', viscosity
      print '(3x,"densityF = ", ES12.4)', rho_f
      print '(3x,"densityP = ", ES12.4)', rho_p

      print '(3x,"DarcyqMa = ", ES12.4)', qMag
      print '(3x,"diffCoef = ", ES12.4)', diffusionCoeff
      print '(3x,"Gm       = ", ES12.4)', Gm
      print '(3x,"Gm5      = ", ES12.4)', Gm5
      print '(3x,"Happel   = ", ES12.4)', Happel
      print '(3x,"NR       = ", ES12.4)', NR
      print '(3x,"NPe      = ", ES12.4)', NPe
      print '(3x,"NvW      = ", ES12.4)', NvdW
      print '(3x,"NGr      = ", ES12.4)', NGr
      print '(3x,"EtaD     = ", ES12.4)', Eta_D
      print '(3x,"EtaI     = ", ES12.4)', Eta_I
      print '(3x,"EtaG     = ", ES12.4)', Eta_G
      print '(3x,"Eta0     = ", ES12.4)', Eta_0
      print '(3x,"katt     = ", ES12.4)', katt
      print *, "--------------------"
    endif

  else
    ! A constant rate of attachment
    katt = this%rate_attachment
  end if

!!!!!!!!!!!!!!!!!!!
! Detachment rate
!!!!!!!!!!!!!!!!!!!
  if (this%rate_detachment < 0.d0) then
    kdet = 0.0
  else
    ! A constant rate of attachment
    kdet = this%rate_detachment
  end if

  RateAtt = 0.0
  RateDet = 0.0
  RateDecayAq = 0.0
  RateDecayIm = 0.0


  ! Build here for attachment/detachment
  ! first-order forward - reverse (A <-> C)
  Rate = katt * Vaq * L_water - kdet * Vim * volume
  RateAtt = stoichVaq * Rate
  RateDet = stoichVim * Rate

  ! Build here for inactivation reactions
  ! first-order (A -> X)
  Rate = decayAq * Vaq * L_water
  RateDecayAq = - Rate

  Rate = decayIm * Vim * volume
  RateDecayIm = - Rate
  !immobile_to_water_vol = material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0
  ! sum_food = rt_auxvar%total(this%D_mobile_id,iphase)*L_water + &
            ! rt_auxvar%immobile(this%D_immobile_id)*volume                                                ! in mol

  ! mu_B_im = this%rate_B_1*Vim* &      ! mol/Ls
			! (sum_food/(sum_food + this%monod_D))* &
			! !(rt_auxvar%total(idof_O2,iphase) / &        !oxygen 
			! !(this%K_O + rt_auxvar%total(idof_O2,iphase)))*&             ! limitation
            ! (this%inhibition_B/ (Vim + this%inhibition_B))**this%exponent_B 
  ! mu_B_mob = this%rate_B_1*.25*Vaq* &      ! mol/Ls
			! (sum_food/(sum_food + this%monod_D))* &
			! !(rt_auxvar%total(idof_O2,iphase) / &        !oxygen 
			! !(this%K_O + rt_auxvar%total(idof_O2,iphase)))*&             ! limitation
            ! (this%inhibition_B/ (Vaq + this%inhibition_B))**this%exponent_B

  ! mu_B_mob_residual = -1* mu_B_mob*L_water + & 
  ! ((this%alpha_vel*global_auxvar%darcy_vel(iphase))**this%beta_vel)* & 
  ! this%rate_B_2*(Vaq - this%background_concentration_B)* L_water  

   ! mu_B_im_residual = -1* mu_B_im*volume + & 
  ! ((this%alpha_vel*global_auxvar%darcy_vel(iphase))**this%beta_vel)* & 
  ! this%rate_B_2*(Vim - this%background_concentration_B)* volume  


  ! mobile_mole_fraction = rt_auxvar%total(this%D_mobile_id,iphase)*L_water/sum_food
  ! immobile_mole_fraction = 1 - mobile_mole_fraction


  ! Residual(this%D_mobile_id) = Residual(this%D_mobile_id) + &
                           ! ! Growth usage, mol/s
                           ! this%stoichiometric_D_1* &
                           ! mobile_mole_fraction* &                                ! dimensionless
                           ! mu_B_mob*L_water + &                        ! mol/m3 bulk/s * m3 bulk
                           ! ! Direct usage, mol/s
                           ! this%rate_D* &                          ! 1/s
                           ! mobile_mole_fraction* &                                ! dimensionless
                           ! Vaq*L_water + &                        ! mol/L/s * m3 bulk
                           ! ! immobilization, mol/s
                           ! this%rate_D_i* &                   ! 1/s
                           ! rt_auxvar%total(this%D_mobile_id,iphase)* &            ! mol/L
                           ! L_water - &   ! L water/m3 bulk
                           ! ! remobilization, mol/s
                           ! this%rate_D_m* &                   ! 1/s
                           ! rt_auxvar%immobile(this%D_immobile_id)* &           ! mol/m3 bulk
                           ! volume                                 ! m3 bulk


  ! Residual(this%D_immobile_id) = Residual(this%D_immobile_id) + &
                           ! ! Growth usage, mol/s
                           ! this%stoichiometric_D_1* &                          ! unitless
                           ! immobile_mole_fraction* &                              ! dimensionless
                           ! mu_B_im*volume + &                        ! mol/m3 bulk/s * m3 bulk
                           ! ! Direct usage, mol/s
                           ! this%rate_D* &                          ! 1/s
                           ! immobile_mole_fraction* &                              ! dimensionless
                           ! Vim* volume - &                              ! L water/m3 bulk
                           ! ! immobilization, mol/s
                           ! this%rate_D_i* &                   ! 1/s
                           ! rt_auxvar%total(this%D_mobile_id,iphase)* &            ! mol/L
                           ! L_water + & ! remobilization, mol/s
                           ! this%rate_D_m* &                   ! 1/s
                           ! rt_auxvar%immobile(this%D_immobile_id)* &           ! mol/m3 bulk
                           ! volume
! This awful block just tries to
! avoid concentrations below 1E-50
! (Is this avoided with TRUNCATE_CONCENTRATION ?) > sure it does
  ! if ( Vaq > 0.0 ) then
  !   if ( Vim > 0.0 ) then
  !     !Do nothing
  !   else if ( Vim <= 0.0 ) then
  !     Vim = 1.0d-50
  !     RateDet = 0.0
  !     RateDecayIm = 0.0
  !   end if
  ! else if ( Vaq <= 0.0 ) then
  !   if ( Vim > 0.0 ) then
  !     Vaq = 1.0d-50
  !     RateAtt = 0.0
  !     RateDecayAq = 0.0
  !   else if ( Vim <= 0.0 ) then
  !     Vim = 1.0d-50
  !     Vaq = 1.0d-50
  !     RateAtt = 0.0
  !     RateDet = 0.0
  !     RateDecayAq = 0.0
  !     RateDecayIm = 0.0
  !   end if
  ! end if

  ! The actual calculation:

  Residual(this%species_Vaq_id) = &
    Residual(this%species_Vaq_id) - RateAtt - RateDecayAq + RateDet !+ mu_B_mob_residual

  Residual(this%species_Vim_id + reaction%offset_immobile) = &
    Residual(this%species_Vim_id + reaction%offset_immobile) &
    - RateDet - RateDecayIm !+ mu_B_im_residual




  delta_volfrac = &
            (Residual(this%species_Vaq_id)+ & 
			Residual(this%species_Vim_id + reaction%offset_immobile)) / &
            (this%density_B*1000.d0) / &                           ! mol/L * L/m3
            volume * &                             ! m3 bulk
            option%tran_dt                                         ! s

  rt_auxvar%mnrl_volfrac(this%biomineral_id) = rt_auxvar%mnrl_volfrac(this%biomineral_id) + &
                                      delta_volfrac

end subroutine BioTH_KineticState

! ************************************************************************** !
subroutine BioTH_Destroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Edwin S
  ! Date: 10/01/2020
  !

  implicit none

  class(reaction_sandbox_bioth_type) :: this

end subroutine BioTH_Destroy

end module Reaction_Sandbox_BioTH_class
