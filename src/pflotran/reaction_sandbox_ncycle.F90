module Reaction_Sandbox_Ncycle_class

!Nitrogen Biogeochemical Transformations
!After Jianqiu Zheng's White Paper
!Current Coeffs are written for rxn stoich for fs0 derived from glucose as organic C source

! Primary Species:
! (1) CH2O (DOC)
! (2) O2
! (3) NO3-
! (4) NO2-
! (5) NH4+
! (6) N2O
! (7) N2(aq)
! (8) CO2
! IMMOBILE
! (9) C5H7O2N (Biomass)


#include "petsc/finclude/petscsys.h"
  use petscsys
  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

! 2. Add module variables here.  Note that one must use the PETSc data types
!    PetscInt, PetscReal, PetscBool to declare variables of type integer
!    float/real*8, and logical respectively.  E.g.
!
! PetscReal, parameter :: formula_weight_of_water = 18.01534d0

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_Ncycle_type

! 3. Add variables/arrays associated with new reaction.
    !character(len=MAXWORDLENGTH) :: species_name
    PetscInt :: species_ch2o_id
    PetscInt :: species_o2_id
    PetscInt :: species_no3_id
    PetscInt :: species_no2_id
    PetscInt :: species_nh4_id
    PetscInt :: species_n2o_id
    PetscInt :: species_n2_id
    PetscInt :: species_co2_id
    PetscInt :: species_c5h7o2n_id

   ! PetscReal :: f1
   !PetscReal :: f2
   ! PetscReal :: f3

    PetscReal :: K1, K2, K3, K4, K5, K6, K7, K8
    PetscReal :: Kd1, Kd2, Kd3, Kd4, Kd5, Kd6, Kd7, Kd8
    PetscReal :: Ka1, Ka2, Ka3, Ka4, Ka5, Ka6, Ka7, Ka8

    PetscReal :: k_deg

    PetscReal :: stoich_r1_ch2o, stoich_r2_ch2o, stoich_r3_ch2o
    PetscReal :: stoich_r4_ch2o, stoich_r5_ch2o, stoich_r6_ch2o
    PetscReal :: stoich_r7_ch2o, stoich_r8_ch2o

    PetscReal :: stoich_r1_o2, stoich_r2_o2, stoich_r3_o2
    PetscReal :: stoich_r4_o2, stoich_r5_o2, stoich_r6_o2
    PetscReal :: stoich_r7_o2, stoich_r8_o2

    PetscReal :: stoich_r1_no3, stoich_r2_no3, stoich_r3_no3
    PetscReal :: stoich_r4_no3, stoich_r5_no3, stoich_r6_no3
    PetscReal :: stoich_r7_no3, stoich_r8_no3

    PetscReal :: stoich_r1_no2, stoich_r2_no2, stoich_r3_no2
    PetscReal :: stoich_r4_no2, stoich_r5_no2, stoich_r6_no2
    PetscReal :: stoich_r7_no2, stoich_r8_no2

    PetscReal :: stoich_r1_nh4, stoich_r2_nh4, stoich_r3_nh4
    PetscReal :: stoich_r4_nh4, stoich_r5_nh4, stoich_r6_nh4
    PetscReal :: stoich_r7_nh4, stoich_r8_nh4

    PetscReal :: stoich_r1_n2o, stoich_r2_n2o, stoich_r3_n2o
    PetscReal :: stoich_r4_n2o, stoich_r5_n2o, stoich_r6_n2o
    PetscReal :: stoich_r7_n2o, stoich_r8_n2o

    PetscReal :: stoich_r1_n2, stoich_r2_n2, stoich_r3_n2
    PetscReal :: stoich_r4_n2, stoich_r5_n2, stoich_r6_n2
    PetscReal :: stoich_r7_n2, stoich_r8_n2

    PetscReal :: stoich_r1_co2, stoich_r2_co2, stoich_r3_co2
    PetscReal :: stoich_r4_co2, stoich_r5_co2, stoich_r6_co2
    PetscReal :: stoich_r7_co2, stoich_r8_co2

    PetscReal :: stoich_r1_c5h7o2n, stoich_r2_c5h7o2n, stoich_r3_c5h7o2n
    PetscReal :: stoich_r4_c5h7o2n, stoich_r5_c5h7o2n, stoich_r6_c5h7o2n
    PetscReal :: stoich_r7_c5h7o2n, stoich_r8_c5h7o2n

    !PetscReal, dimension(9,8) :: stoich

    PetscReal :: ref_temp = 298.15d0

  contains
    procedure, public :: ReadInput => ExampleRead
    procedure, public :: Setup => ExampleSetup
    procedure, public :: Evaluate => ExampleEvaluate
    procedure, public :: Destroy => ExampleDestroy
  end type reaction_sandbox_Ncycle_type

  public :: NcycleCreate

contains

! ************************************************************************** !

function NcycleCreate()
  !
  ! Allocates example reaction object.
  !
  ! Author: John Doe (replace in all subroutine headers with name of developer)
  ! Date: 00/00/00 (replace in all subroutine headers with current date)
  !

  implicit none

  class(reaction_sandbox_Ncycle_type), pointer :: NcycleCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(NcycleCreate)
  !NcycleCreate%species_name = ''
  NcycleCreate%species_ch2o_id = UNINITIALIZED_INTEGER
  NcycleCreate%species_o2_id = UNINITIALIZED_INTEGER
  NcycleCreate%species_no3_id = UNINITIALIZED_INTEGER
  NcycleCreate%species_no2_id = UNINITIALIZED_INTEGER
  NcycleCreate%species_nh4_id = UNINITIALIZED_INTEGER
  NcycleCreate%species_n2o_id = UNINITIALIZED_INTEGER
  NcycleCreate%species_n2_id = UNINITIALIZED_INTEGER
  NcycleCreate%species_co2_id = UNINITIALIZED_INTEGER

  NcycleCreate%species_c5h7o2n_id = UNINITIALIZED_INTEGER

  !NcycleCreate%f1 = UNINITIALIZED_DOUBLE
  !NcycleCreate%f2 = UNINITIALIZED_DOUBLE
 ! NcycleCreate%f3 = UNINITIALIZED_DOUBLE

  NcycleCreate%K1 = UNINITIALIZED_DOUBLE
  NcycleCreate%K2 = UNINITIALIZED_DOUBLE
  NcycleCreate%K3 = UNINITIALIZED_DOUBLE
  NcycleCreate%K4 = UNINITIALIZED_DOUBLE
  NcycleCreate%K5 = UNINITIALIZED_DOUBLE
  NcycleCreate%K6 = UNINITIALIZED_DOUBLE
  NcycleCreate%K7 = UNINITIALIZED_DOUBLE
  NcycleCreate%K8 = UNINITIALIZED_DOUBLE

  NcycleCreate%Kd1 = UNINITIALIZED_DOUBLE
  NcycleCreate%Kd2 = UNINITIALIZED_DOUBLE
  NcycleCreate%Kd3 = UNINITIALIZED_DOUBLE
  NcycleCreate%Kd4 = UNINITIALIZED_DOUBLE
  NcycleCreate%Kd5 = UNINITIALIZED_DOUBLE
  NcycleCreate%Kd6 = UNINITIALIZED_DOUBLE
  NcycleCreate%Kd7 = UNINITIALIZED_DOUBLE
  NcycleCreate%Kd8 = UNINITIALIZED_DOUBLE

  NcycleCreate%Ka1 = UNINITIALIZED_DOUBLE
  NcycleCreate%Ka2 = UNINITIALIZED_DOUBLE
  NcycleCreate%Ka3 = UNINITIALIZED_DOUBLE
  NcycleCreate%Ka4 = UNINITIALIZED_DOUBLE
  NcycleCreate%Ka5 = UNINITIALIZED_DOUBLE
  NcycleCreate%Ka6 = UNINITIALIZED_DOUBLE
  NcycleCreate%Ka7 = UNINITIALIZED_DOUBLE
  NcycleCreate%Ka8 = UNINITIALIZED_DOUBLE

  NcycleCreate%k_deg = UNINITIALIZED_DOUBLE

  NcycleCreate%stoich_r1_ch2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r1_o2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r1_no3 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r1_no2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r1_nh4 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r1_n2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r1_n2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r1_co2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r1_c5h7o2n = UNINITIALIZED_DOUBLE

  NcycleCreate%stoich_r2_ch2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r2_o2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r2_no3 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r2_no2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r2_nh4 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r2_n2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r2_n2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r2_co2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r2_c5h7o2n = UNINITIALIZED_DOUBLE

  NcycleCreate%stoich_r3_ch2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r3_o2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r3_no3 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r3_no2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r3_nh4 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r3_n2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r3_n2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r3_co2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r3_c5h7o2n = UNINITIALIZED_DOUBLE

  NcycleCreate%stoich_r4_ch2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r4_o2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r4_no3 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r4_no2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r4_nh4 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r4_n2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r4_n2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r4_co2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r4_c5h7o2n = UNINITIALIZED_DOUBLE

  NcycleCreate%stoich_r5_ch2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r5_o2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r5_no3 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r5_no2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r5_nh4 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r5_n2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r5_n2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r5_co2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r5_c5h7o2n = UNINITIALIZED_DOUBLE

  NcycleCreate%stoich_r6_ch2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r6_o2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r6_no3 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r6_no2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r6_nh4 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r6_n2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r6_n2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r6_co2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r6_c5h7o2n = UNINITIALIZED_DOUBLE

  NcycleCreate%stoich_r7_ch2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r7_o2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r7_no3 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r7_no2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r7_nh4 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r7_n2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r7_n2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r7_co2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r7_c5h7o2n = UNINITIALIZED_DOUBLE

  NcycleCreate%stoich_r8_ch2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r8_o2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r8_no3 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r8_no2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r8_nh4 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r8_n2o = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r8_n2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r8_co2 = UNINITIALIZED_DOUBLE
  NcycleCreate%stoich_r8_c5h7o2n = UNINITIALIZED_DOUBLE

  nullify(NcycleCreate%next)

end function NcycleCreate

! ************************************************************************** !

subroutine ExampleRead(this,input,option)
  !!
  !! Reads input deck for example reaction parameters (if any)
  !!
  !! Author: John Doe
  !! Date: 00/00/00
  !!
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_Ncycle_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units, units
  character(len=MAXWORDLENGTH) :: error_string

  error_string = 'CHEMSITRY,REACTION_SANDBOX,NCYCLE'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    select case(trim(word))
#if 0
      case('F1')
        call InputReadDouble(input,option,this%f1)
        call InputErrorMsg(input,option,'f1',error_string)

      case('F2')
        call InputReadDouble(input,option,this%f2)
        call InputErrorMsg(input,option,'f2',error_string)

       case('F3')
        call InputReadDouble(input,option,this%f3)
        call InputErrorMsg(input,option,'f3',error_string)
#endif
       case('K1')
        call InputReadDouble(input,option,this%K1)
        call InputErrorMsg(input,option,'K1',error_string)
        call InputReadAndConvertUnits(input,this%K1,'1/sec',&
	      trim(error_string)//',K1',option)

       case('K2')
        call InputReadDouble(input,option,this%K2)
        call InputErrorMsg(input,option,'K2',error_string)
        call InputReadAndConvertUnits(input,this%K2,'1/sec',&
	      trim(error_string)//',K2',option)

       case('K3')
        call InputReadDouble(input,option,this%K3)
        call InputErrorMsg(input,option,'K3',error_string)
        call InputReadAndConvertUnits(input,this%K3,'1/sec',&
	      trim(error_string)//',K3',option)

       case('K4')
        call InputReadDouble(input,option,this%K4)
        call InputErrorMsg(input,option,'K4',error_string)
        call InputReadAndConvertUnits(input,this%K4,'1/sec',&
	      trim(error_string)//',K4',option)

       case('K5')
        call InputReadDouble(input,option,this%K5)
        call InputErrorMsg(input,option,'K5',error_string)
        call InputReadAndConvertUnits(input,this%K5,'1/sec',&
	      trim(error_string)//',K5',option)

       case('K6')
        call InputReadDouble(input,option,this%K6)
        call InputErrorMsg(input,option,'K6',error_string)
        call InputReadAndConvertUnits(input,this%K6,'1/sec',&
	      trim(error_string)//',K6',option)

       case('K7')
        call InputReadDouble(input,option,this%K7)
        call InputErrorMsg(input,option,'K7',error_string)
        call InputReadAndConvertUnits(input,this%K7,'1/sec',&
	      trim(error_string)//',K7',option)

       case('K8')
        call InputReadDouble(input,option,this%K8)
        call InputErrorMsg(input,option,'K8',error_string)
        call InputReadAndConvertUnits(input,this%K8,'1/sec',&
	      trim(error_string)//',K8',option)

      case('KD1')
        call InputReadDouble(input,option,this%Kd1)
        call InputErrorMsg(input,option,'Kd1',error_string)
        call InputReadAndConvertUnits(input,this%Kd1,'M',&
	      trim(error_string)//',Kd1',option)

      case('KD2')
        call InputReadDouble(input,option,this%Kd2)
        call InputErrorMsg(input,option,'Kd2',error_string)
        call InputReadAndConvertUnits(input,this%Kd2,'M',&
	      trim(error_string)//',Kd2',option)

      case('KD3')
        call InputReadDouble(input,option,this%Kd3)
        call InputErrorMsg(input,option,'Kd3',error_string)
        call InputReadAndConvertUnits(input,this%Kd3,'M',&
	      trim(error_string)//',Kd3',option)

        case('KD4')
        call InputReadDouble(input,option,this%Kd4)
        call InputErrorMsg(input,option,'Kd4',error_string)
        call InputReadAndConvertUnits(input,this%Kd4,'M',&
	      trim(error_string)//',Kd4',option)

        case('KD5')
        call InputReadDouble(input,option,this%Kd5)
        call InputErrorMsg(input,option,'Kd5',error_string)
        call InputReadAndConvertUnits(input,this%Kd5,'M',&
	      trim(error_string)//',Kd5',option)

        case('KD6')
        call InputReadDouble(input,option,this%Kd6)
        call InputErrorMsg(input,option,'Kd6',error_string)
        call InputReadAndConvertUnits(input,this%Kd6,'M',&
	      trim(error_string)//',Kd6',option)

        case('KD7')
        call InputReadDouble(input,option,this%Kd7)
        call InputErrorMsg(input,option,'Kd7',error_string)
        call InputReadAndConvertUnits(input,this%Kd7,'M',&
	      trim(error_string)//',Kd7',option)

        case('KD8')
        call InputReadDouble(input,option,this%Kd8)
        call InputErrorMsg(input,option,'Kd8',error_string)
        call InputReadAndConvertUnits(input,this%Kd8,'M',&
	      trim(error_string)//',Kd8',option)


      case('KA1')
        call InputReadDouble(input,option,this%Ka1)
        call InputErrorMsg(input,option,'Ka1', error_string)
        call InputReadAndConvertUnits(input,this%Ka1,'M',&
	      trim(error_string)//',Ka1',option)

      case('KA2')
        call InputReadDouble(input,option,this%Ka2)
        call InputErrorMsg(input,option,'Ka2',error_string)
        call InputReadAndConvertUnits(input,this%Ka2,'M',&
	      trim(error_string)//',Ka2',option)

      case('KA3')
        call InputReadDouble(input,option,this%Ka3)
        call InputErrorMsg(input,option,'Ka3',error_string)
        call InputReadAndConvertUnits(input,this%Ka3,'M',&
	      trim(error_string)//',Ka3',option)

      case('KA4')
        call InputReadDouble(input,option,this%Ka4)
        call InputErrorMsg(input,option,'Ka4',error_string)
        call InputReadAndConvertUnits(input,this%Ka4,'M',&
	      trim(error_string)//',Ka4',option)

      case('KA5')
        call InputReadDouble(input,option,this%Ka5)
        call InputErrorMsg(input,option,'Ka5',error_string)
        call InputReadAndConvertUnits(input,this%Ka5,'M',&
	      trim(error_string)//',Ka5',option)

      case('KA6')
        call InputReadDouble(input,option,this%Ka6)
        call InputErrorMsg(input,option,'Ka6',error_string)
        call InputReadAndConvertUnits(input,this%Ka6,'M',&
	      trim(error_string)//',Ka6',option)

      case('KA7')
        call InputReadDouble(input,option,this%Ka7)
        call InputErrorMsg(input,option,'Ka7',error_string)
        call InputReadAndConvertUnits(input,this%Ka7,'M',&
	      trim(error_string)//',Ka7',option)

      case('KA8')
        call InputReadDouble(input,option,this%Ka8)
        call InputErrorMsg(input,option,'Ka8',error_string)
        call InputReadAndConvertUnits(input,this%Ka8,'M',&
	      trim(error_string)//',Ka8',option)

      case('K_DEG')
        call InputReadDouble(input,option,this%k_deg)
        call InputErrorMsg(input,option,'k_deg',error_string)
        call InputReadAndConvertUnits(input,this%k_deg,'1/sec',&
			trim(error_string)//',k_deg',option)

      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine ExampleRead

! ************************************************************************** !

subroutine ExampleSetup(this,reaction,option)
  !
  ! Sets up the example reaction either with parameters either
  ! read from the input deck or hardwired.
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !

  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_Ncycle_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize
 character(len=MAXWORDLENGTH) :: word

 word = 'CH2O(aq)'
 this%species_ch2o_id = &
   GetPrimarySpeciesIDFromName(word,reaction,option)

 word = 'O2(aq)'
 this%species_o2_id = &
   GetPrimarySpeciesIDFromName(word,reaction,option)

 word = 'NO3-'
 this%species_no3_id = &
   GetPrimarySpeciesIDFromName(word,reaction,option)

 word = 'NO2-'
 this%species_no2_id = &
   GetPrimarySpeciesIDFromName(word,reaction,option)

 word = 'NH4+'
 this%species_nh4_id = &
   GetPrimarySpeciesIDFromName(word,reaction,option)

 word = 'N2O'
 this%species_n2O_id = &
   GetPrimarySpeciesIDFromName(word,reaction,option)

 word = 'N2(aq)'
 this%species_n2_id = &
   GetPrimarySpeciesIDFromName(word,reaction,option)

 word = 'CO2(aq)'
 this%species_co2_id = &
   GetPrimarySpeciesIDFromName(word,reaction,option)

 word = 'C5H7O2N'
 this%species_c5h7o2n_id = &
   GetImmobileSpeciesIDFromName(word,reaction%immobile,option)  !! Check this

 this%stoich_r1_ch2o = -1.39d0
 this%stoich_r2_ch2o = 1.d-12
 this%stoich_r3_ch2o = 1.d-12
 this%stoich_r4_ch2o = -1.37d0
 this%stoich_r5_ch2o = -1.3d0
 this%stoich_r6_ch2o = -1.27d0
 this%stoich_r7_ch2o = -1.88d0
 this%stoich_r8_ch2o = 1.d-12

 this%stoich_r1_o2 = -0.39d0
 this%stoich_r2_o2 = -9.53d0
 this%stoich_r3_o2 = -14.16d0
 this%stoich_r4_o2 = 1.d-12
 this%stoich_r5_o2 = 1.d-12
 this%stoich_r6_o2 = 1.d-12
 this%stoich_r7_o2 = 1.d-12
 this%stoich_r8_o2 = 1.d-12

 this%stoich_r1_no3 = 1.d-20
 this%stoich_r2_no3 = 1.d-12
 this%stoich_r3_no3 = 30.91d0
 this%stoich_r4_no3 = -0.74d0
 this%stoich_r5_no3 = 1.d-12
 this%stoich_r6_no3 = 1.d-12
 this%stoich_r7_no3 = 1.d-12
 this%stoich_r8_no3 = 1.d-12

 this%stoich_r1_no2 = 1.d-20
 this%stoich_r2_no2 = 7.02d0
 this%stoich_r3_no2 = -31.11d0
 this%stoich_r4_no2 = 0.74d0
 this%stoich_r5_no2 = -0.6d0
 this%stoich_r6_no2 = 1.d-12
 this%stoich_r7_no2 = -0.59d0
 this%stoich_r8_no2 = -5.44d0

 this%stoich_r1_nh4 = -0.2d0
 this%stoich_r2_nh4 = -7.22d0
 this%stoich_r3_nh4 = 1.d-12
 this%stoich_r4_nh4 = -0.2d0
 this%stoich_r5_nh4 = -0.2d0
 this%stoich_r6_nh4 = -0.2d0
 this%stoich_r7_nh4 = 0.39d0
 this%stoich_r8_nh4 = -6.97d0

 this%stoich_r1_n2o = 1.d-20
 this%stoich_r2_n2o = 1.d-12
 this%stoich_r3_n2o = 1.d-12
 this%stoich_r4_n2o = 1.d-12
 this%stoich_r5_n2o = 0.3d0
 this%stoich_r6_n2o = -0.53d0
 this%stoich_r7_n2o = 1.d-12
 this%stoich_r8_n2o = 1.d-12

 this%stoich_r1_n2 = 1.d-20
 this%stoich_r2_n2 = 1.d-20
 this%stoich_r3_n2 = 1.d-12
 this%stoich_r4_n2 = 1.d-12
 this%stoich_r5_n2 = 1.d-12
 this%stoich_r6_n2 = 0.53d0
 this%stoich_r7_n2 = 1.d-12
 this%stoich_r8_n2 = 6.1d0

 this%stoich_r1_co2 = 0.39d0
 this%stoich_r2_co2 = -1.d0
 this%stoich_r3_co2 = -1.d0
 this%stoich_r4_co2 = 0.37d0
 this%stoich_r5_co2 = 0.3d0
 this%stoich_r6_co2 = 0.27d0
 this%stoich_r7_co2 = 0.88d0
 this%stoich_r8_co2 = -1.d0

 this%stoich_r1_c5h7o2n = 0.2d0
 this%stoich_r2_c5h7o2n = 0.2d0
 this%stoich_r3_c5h7o2n = 0.2d0
 this%stoich_r4_c5h7o2n = 0.2d0
 this%stoich_r5_c5h7o2n = 0.2d0
 this%stoich_r6_c5h7o2n = 0.2d0
 this%stoich_r7_c5h7o2n = 0.2d0
 this%stoich_r8_c5h7o2n = 0.2d0

end subroutine ExampleSetup

! ************************************************************************** !

subroutine ExampleEvaluate(this,Residual,Jacobian,compute_derivative, &
                           rt_auxvar,global_auxvar,material_auxvar,reaction, &
                           option)
  !
  ! Evaluates reaction storing residual and/or Jacobian
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class

  implicit none

  class(reaction_sandbox_Ncycle_type) :: this
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
  PetscReal :: L_water
  PetscReal :: kg_water
  PetscReal :: volume
  PetscReal :: porosity
  PetscReal :: liquid_saturation
  PetscReal :: molality_to_molarity

  PetscReal :: C_ch2o, C_o2, C_no3, C_no2, C_nh4, C_n2o, C_n2, C_co2
  PetscReal :: C_c5h7o2n
  PetscReal :: r1kin, e1rel, r1
  PetscReal :: r2kin, e2rel, r2
  PetscReal :: r3kin, e3rel, r3
  PetscReal :: r4kin, e4rel, r4
  PetscReal :: r5kin, e5rel, r5
  PetscReal :: r6kin, e6rel, r6
  PetscReal :: r7kin, e7rel, r7
  PetscReal :: r8kin, e8rel, r8
  PetscReal :: sumkin
  PetscReal :: Rate_ch2o, Rate_o2, Rate_no3, Rate_no2, Rate_nh4, Rate_n2o
  PetscReal :: Rate_n2, Rate_co2
  PetscReal :: Rate_c5h7o2n


  porosity = material_auxvar%porosity               ! [m^3 pore/m^3 bulk volume]
  liquid_saturation = global_auxvar%sat(iphase)     ! [m^3 water/m^3 pore]
  volume = material_auxvar%volume                   ! [m^3 bulk volume]
            ! multiplying by 1.d3 converts m^3 water -> L water

  L_water = porosity*liquid_saturation*volume*1.d3
  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3

  C_ch2o = rt_auxvar%pri_molal(this%species_ch2o_id)*molality_to_molarity
  C_o2 = rt_auxvar%pri_molal(this%species_o2_id)*molality_to_molarity
  C_no3 = rt_auxvar%pri_molal(this%species_no3_id)*molality_to_molarity
  C_no2 = rt_auxvar%pri_molal(this%species_no2_id)*molality_to_molarity
  C_nh4 = rt_auxvar%pri_molal(this%species_nh4_id)*molality_to_molarity
  C_n2o = rt_auxvar%pri_molal(this%species_n2o_id)*molality_to_molarity
  C_n2 = rt_auxvar%pri_molal(this%species_n2_id)*molality_to_molarity
  C_co2 = rt_auxvar%pri_molal(this%species_co2_id)*molality_to_molarity

  C_c5h7o2n = rt_auxvar%immobile(this%species_c5h7o2n_id) ![mol/m^3 bulk volume]

 ! Unregulated reaction rates based on Monod kinetics
 ! CHECK EXPRESSIONS--- No Biomass term in rxns in white paper ----
 ! UNITS:  K1 = u,max [=] umol/L-sec ; Ka,Kd [=] umol/L; C_a, C_d [=] umol/L
 !         rkin [=] umol/sec
 !CHECK: Does L_water mean units of L !!!!
  r1kin = 0.d0
  r2kin = 0.d0
  r3kin = 0.d0
  r4kin = 0.d0
  r5kin = 0.d0
  r6kin = 0.d0
  r7kin = 0.d0
  r8kin = 0.d0


  r1kin = this%K1 * (C_ch2o / (this%Kd1 + C_ch2o)) * &
               (C_o2 / (this%Ka1 + C_o2)) * L_water ![mol/sec]

  r2kin = this%K2 * (C_nh4 / (this%Kd2 + C_nh4)) * &
               (C_o2 / (this%Ka2 + C_o2)) * L_water ![mol/sec]

  r3kin = this%K3 *  (C_no2 / (this%Kd3 + C_no2)) &
             * (C_o2 / (this%Ka3 + C_o2)) * L_water ![mol/sec]

  r4kin = this%K4 * (C_ch2o / (this%Kd4 + C_ch2o)) &
             * (C_no3 / (this%Ka4 + C_no3)) * L_water ![mol/sec]

  r5kin = this%K5 * (C_ch2o / (this%Kd5 + C_ch2o)) &
             * (C_no2 / (this%Ka5 + C_no2)) * L_water ![mol/sec]

  r6kin = this%K6 * (C_ch2o / (this%Kd6 + C_ch2o)) &
             * (C_n2o / (this%Ka6 + C_n2o)) * L_water ![mol/sec]

  r7kin = this%K7 * (C_ch2o / (this%Kd7 + C_ch2o)) &
             * (C_no2 / (this%Ka7 + C_no2)) * L_water ![mol/sec]

  r8kin = this%K8 *  (C_nh4 / (this%Kd8 + C_nh4)) &
             * (C_no2 / (this%Ka8 + C_no2)) * L_water ![mol/sec]
#if 0
#endif

! Relative contribution (unitless)
  sumkin = r1kin + r2kin + r3kin + r4kin + r5kin + r6kin + r7kin + r8kin
  e1rel = r1kin / sumkin
  e2rel = r2kin / sumkin
  e3rel = r3kin / sumkin
  e4rel = r4kin / sumkin
  e5rel = r5kin / sumkin
  e6rel = r6kin / sumkin
  e7rel = r7kin / sumkin
  e8rel = r8kin / sumkin

!print *, 'sumkin: ', sumkin
!stop
! Relative rates with enzyme levels considered
  r1 = r1kin * e1rel
  r2 = r2kin * e2rel
  r3 = r3kin * e3rel
  r4 = r4kin * e4rel
  r5 = r5kin * e5rel
  r6 = r6kin * e6rel
  r7 = r7kin * e7rel
  r8 = r8kin * e8rel

! Rates
! Reactions are modulated by BM concentration (i.e., C5H7O2N)
  ! If BM mobile, then units are: [mol/mol-BM-s * mol-BM/L-water *L-water] = mol/sec
  ! If BM immobile, then units are: [mol/mol-BM-s * mol-BM/m^3-bulk *Volume] = mol/sec

 ! DOC (i.e., CH2O) includes +kdeg on biomass (growth) (based on expressions in White Paper)
 Rate_ch2o = (this%stoich_r1_ch2o * r1 + this%stoich_r2_ch2o  * r2 &
    + this%stoich_r3_ch2o * r3 + this%stoich_r4_ch2o * r4 &
    + this%stoich_r5_ch2o * r5 + this%stoich_r6_ch2o * r6 &
    + this%stoich_r7_ch2o * r7 + this%stoich_r8_ch2o * r8 &
    + this%k_deg) * C_c5h7o2n * L_water

 Rate_o2 = (this%stoich_r1_o2 * r1 + this%stoich_r2_o2* r2 &
    + this%stoich_r3_o2 * r3  + this%stoich_r4_o2 * r4 &
    + this%stoich_r5_o2 * r5  + this%stoich_r6_o2 * r6 &
    + this%stoich_r7_o2 * r7  + this%stoich_r8_o2 * r8) &
    * C_c5h7o2n * L_water

 Rate_no3 = (this%stoich_r1_no3  * r1 + this%stoich_r2_no3 * r2 &
   + this%stoich_r3_no3 * r3 + this%stoich_r4_no3 * r4 &
   + this%stoich_r5_no3 * r5 + this%stoich_r6_no3 * r6 &
   + this%stoich_r7_no3 * r7 + this%stoich_r8_no3 * r8) &
   * C_c5h7o2n * L_water

 Rate_no2 = (this%stoich_r1_no2 * r1 + this%stoich_r2_no2 * r2 &
   + this%stoich_r3_no2 * r3  + this%stoich_r4_no2 * r4 &
   + this%stoich_r5_no2 * r5  + this%stoich_r6_no2 * r6 &
   + this%stoich_r7_no2 * r7  + this%stoich_r8_no2 * r8) &
   * C_c5h7o2n * L_water

 Rate_nh4 = (this%stoich_r1_nh4 * r1 + this%stoich_r2_nh4 * r2 &
   + this%stoich_r3_nh4 * r3 + this%stoich_r4_nh4 * r4 &
   + this%stoich_r5_nh4 * r5 + this%stoich_r6_nh4 * r6 &
   + this%stoich_r7_nh4 * r7 + this%stoich_r8_nh4 * r8) &
   * C_c5h7o2n * L_water

 Rate_n2o = (this%stoich_r1_n2o * r1 + this%stoich_r2_n2o * r2 &
   + this%stoich_r3_n2o * r3 + this%stoich_r4_n2o * r4 &
   + this%stoich_r5_n2o * r5 + this%stoich_r6_n2o * r6 &
   + this%stoich_r7_n2o * r7 + this%stoich_r8_n2o * r8) &
   * C_c5h7o2n * L_water

 Rate_n2 = (this%stoich_r1_n2 * r1 + this%stoich_r2_n2 * r2 &
   + this%stoich_r3_n2 * r3 + this%stoich_r4_n2 * r4 &
   + this%stoich_r5_n2 * r5 + this%stoich_r6_n2 * r6 &
   + this%stoich_r7_n2 * r7 + this%stoich_r8_n2 * r8) &
   * C_c5h7o2n * L_water

 Rate_co2 = (this%stoich_r1_co2 * r1 + this%stoich_r2_co2 * r2 &
  + this%stoich_r3_co2* r3 + this%stoich_r4_co2* r4 &
  + this%stoich_r5_co2* r5 + this%stoich_r6_co2* r6 &
  + this%stoich_r7_co2* r7 + this%stoich_r8_co2* r8) &
  * C_c5h7o2n * L_water

! Biomass Rate
! DOC (i.e., CH2O) includes -kdeg on biomass (decay) (based on expressions in White Paper)
 Rate_c5h7o2n=(this%stoich_r1_c5h7o2n * r1 + this%stoich_r2_c5h7o2n * r2 &
     + this%stoich_r3_c5h7o2n * r3 + this%stoich_r4_c5h7o2n * r4 &
     + this%stoich_r5_c5h7o2n * r5 + this%stoich_r6_c5h7o2n * r6 &
     + this%stoich_r7_c5h7o2n * r7 + this%stoich_r8_c5h7o2n * r8 &
     - this%k_deg) * C_c5h7o2n * L_water

! Residuals
 Residual(this%species_ch2o_id) = Residual(this%species_ch2o_id) - Rate_ch2o
 Residual(this%species_o2_id) = Residual(this%species_o2_id) - Rate_o2
 Residual(this%species_no3_id) = Residual(this%species_no3_id) - Rate_no3
 Residual(this%species_no2_id) = Residual(this%species_no2_id) - Rate_no2
 Residual(this%species_nh4_id) = Residual(this%species_nh4_id) - Rate_nh4
 Residual(this%species_n2o_id) = Residual(this%species_n2o_id) - Rate_n2o
 Residual(this%species_n2_id) = Residual(this%species_n2_id) - Rate_n2
 Residual(this%species_co2_id) = Residual(this%species_co2_id) - Rate_co2

 Residual(this%species_c5h7o2n_id + reaction%offset_immobile) = &
    Residual(this%species_c5h7o2n_id + reaction%offset_immobile) - Rate_c5h7o2n


  ! Description of subroutine arguments:

  ! Residual - 1D array storing residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entires in units kg water/sec
  !
  !  Jacobian [kg water/sec] * dc [mol/kg water] = -Res [mol/sec]
  !
  ! compute_derivative - Flag indicating whether analtical derivative should
  !   be calculated.  The user must provide either the analytical derivatives
  !   or a numerical approximation unless always running with
  !   NUMERICAL_JACOBIAN_RXN defined in input deck.  If the use of
  !   NUMERICAL_JACOBIAN_RXN is assumed, the user should provide an error
  !   message when compute_derivative is true.  E.g.
  !
  !   option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
  !                      'due to assumptions in Example'
  !   call PrintErrMsg(option)
  !
  ! rt_auxvar - Object holding chemistry information (e.g. concentrations,
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
  ! global_auxvar - Object holding information on flow (e.g. saturation,
  !   density, viscosity, temperature, etc)
  !
  !   Useful variables:
  !     global_auxvar%den(iphase) - fluid density [mol/m^3]
  !     global_auxvar%den_kg(iphase) - fluid density [kg/m^3]
  !     global_auxvar%sat(iphase) - saturation
  !     global_auxvar%temp - temperature [C]
  !
  ! porosity - effective porosity of grid cell [m^3 pore/m^3 bulk]
  ! volume - volume of grid cell [m^3]
  ! reaction - Provides access to variable describing chemistry.  E.g.
  !   reaction%ncomp - # chemical degrees of freedom (mobile and immobile)
  !   reaction%naqcomp - # chemical degrees of freedom on water
  !   reaction%primary_species_names(:) - names of primary species
  !
  ! option - Provides handle for controlling simulation, catching and
  !          reporting errors.

! 10. Add code for residual evaluation

  ! Unit of the residual must be in moles/second
  ! global_auxvar%sat(iphase) = saturation of cell
  ! 1.d3 converts m^3 water -> L water
  !L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
  !          material_auxvar%volume*1.d3
  ! always subtract contribution from residual


 !Residual(this%species_id) = Residual(this%species_id) - &
    !(-1.d0) * & ! negative stoichiometry
    !this%rate_constant * &  ! 1/sec
    !L_water * & ! L water
    !! rt_auxvar%total(this%species_id,iphase) = species total component
    !!   concentration
    !rt_auxvar%total(this%species_id,iphase) ! mol/L water



!  if (compute_derivative) then

!! 11. If using an analytical Jacobian, add code for Jacobian evaluation

    ! always add contribution to Jacobian
    ! units = (mol/sec)*(kg water/mol) = kg water/sec
    !Jacobian(this%species_id,this%species_id) = &
    !Jacobian(this%species_id,this%species_id) - &
    !  (-1.d0) * & ! negative stoichiometry
    !  this%rate_constant * & ! 1/sec
    !  L_water * & ! L water
                  ! kg water/L water
      ! rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) =
      !   derivative of total component concentration with respect to the
      !   free ion concentration of the same species.
    !  rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase)

  !endif

end subroutine ExampleEvaluate

! ************************************************************************** !

subroutine ExampleDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: John Doe
  ! Date: 00/00/00
  !

  implicit none

  class(reaction_sandbox_Ncycle_type) :: this

! 12. Add code to deallocate contents of the example object

end subroutine ExampleDestroy

end module Reaction_Sandbox_Ncycle_class
