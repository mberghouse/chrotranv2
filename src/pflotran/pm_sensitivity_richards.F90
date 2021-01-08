module PM_Sensitivity_Richards_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PM_Base_class
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module
  use Output_Sensitivity_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_sensitivity_richards_type
    class(realization_subsurface_type), pointer :: realization
    class(output_sensitivity_option_type), pointer :: output_sensitivity_option
  contains
    procedure, public :: SetRealization => PMSensitivityRichardsSetRealization
    procedure, public :: Setup => PMSensitivityRichardsSetup
    procedure, public :: ReadPMBlock => PMSensitivityRichardsReadPMBlock
    procedure, public :: InitializeRun => PMSensitivityRichardsInitializeRun
    procedure, public :: InitializeTimestep => PMSensitivityRichardsInitializeTimestep
    procedure, public :: FinalizeTimestep => PMSensitivityRichardsFinalizeTimestep
    procedure, public :: Solve => PMSensitivityRichardsSolve
    procedure, public :: Output => PMSensitivityRichardsOutput
    procedure, public :: InputRecord => PMSensitivityRichardsInputRecord
    procedure, public :: Destroy => PMSensitivityRichardsDestroy
  end type pm_sensitivity_richards_type

  public :: PMSensitivityRichardsCreate
  
contains

! *************************************************************************** !

function PMSensitivityRichardsCreate()
  !
  ! Creates and initializes the sensitivity Richards process model.
  !
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  implicit none
  
  class(pm_sensitivity_richards_type), pointer :: PMSensitivityRichardsCreate
  
  allocate(PMSensitivityRichardsCreate)
  call PMBaseInit(PMSensitivityRichardsCreate)
  nullify(PMSensitivityRichardsCreate%realization)
  PMSensitivityRichardsCreate%name = 'sensitivity Richards'
  PMSensitivityRichardsCreate%header = 'SENSITIVITY RICHARDS'
  allocate(PMSensitivityRichardsCreate%output_sensitivity_option)
  !PMSensitivityRichardsCreate%output_sensitivity_option => &
  call OutputSensitivityOptionInit( &
                         PMSensitivityRichardsCreate%output_sensitivity_option)
  
end function PMSensitivityRichardsCreate

! *************************************************************************** !

subroutine PMSensitivityRichardsSetRealization(this,realization)
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !

  use Realization_Subsurface_class

  implicit none
  
  class(pm_sensitivity_richards_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMSensitivityRichardsSetRealization

! *************************************************************************** !

subroutine PMSensitivityRichardsReadPMBlock(this,input)
  !
  ! Reads input file parameters for the sensitivity richards
  !
  ! Author: Moise Rousseau
  ! Date: 03/13/2017
  !
#if 0
SENSTIVITITY_RICHARDS
  FORMAT HDF5
  VARIABLES
    PRESSURE
    PERMEABILITY
    POROSITY
  /
  OUTPUT_TIME SYNC_WITH_SNAPSHOT_FILE
END
#endif
  
  use Input_Aux_module
  use Option_module
  use String_module
  use Output_Sensitivity_module
  
  implicit none
  
  class(pm_sensitivity_richards_type) :: this
  type(input_type), pointer :: input
  
  class(output_sensitivity_option_type), pointer :: output_sensitivity_option
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: ierr
  
  output_sensitivity_option => this%output_sensitivity_option
  option => this%option
  input%ierr = 0
  
  option%io_buffer = 'pflotran card:: SENSITIVITY RICHARDS'
  call PrintMsg(option)
  
  error_string = 'SENSITIVITY_RICHARDS'
  call InputPushBlock(input,option)
  
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    select case(trim(word))
    !-----------------------------------------
    !-----------------------------------------
      case('FORMAT') !HDF5 or petsc viewer (ascii, binary)
        call InputReadCard(input,option,word)
        call StringToUpper(word)
        call OutputSensitivityOptionSetOutputFormat(output_sensitivity_option,&
                                              option,word)
    !-----------------------------------------
      case('VARIABLES')
        call PMSensitivityRichardsReadOutputVariables(this, &
                                              output_sensitivity_option,input)
    !-----------------------------------------
      case('OUTPUT_TIME')
        call InputReadCard(input,option,word)
        call StringToUpper(word)
        call OutputSensitivityOptionSetOutputTimeOption( &
                                                  output_sensitivity_option,&
                                                  option,word)
    !-----------------------------------------
      case default
        call InputKeywordUnrecognized(input,word,'SENSITIVITY_RICHARDS', &
                                      option)
    !-----------------------------------------
    end select  
  enddo
  call InputPopBlock(input,option)
  
end subroutine PMSensitivityRichardsReadPMBlock

! *************************************************************************** !

subroutine PMSensitivityRichardsReadOutputVariables(this, &
                                            output_sensitivity_option,input)
  !
  ! Sets up the process model with external information.
  !
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !

  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  class(pm_sensitivity_richards_type) :: this
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  
  option => this%option
  
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'SENSITIVITY_RICHARDS, VARIABLES')
    call StringToUpper(word)
    select case(trim(word))
      case('PRESSURE')
      case('POROSITY')
      case('PERMEABILITY')
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'SENSITIVITY_RICHARDS, VARIABLES', &
                                      option)
    end select
    call OutputSensitivityAddVariable(output_sensitivity_option,word)
  end do
  
  call InputPopBlock(input,option)
  
end subroutine PMSensitivityRichardsReadOutputVariables

! *************************************************************************** !

subroutine PMSensitivityRichardsSetup(this)
  !
  ! Sets up the process model with external information.
  !
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !

  implicit none
  
  class(pm_sensitivity_richards_type) :: this
  
  ! TODO
  
end subroutine PMSensitivityRichardsSetup

! *************************************************************************** !

subroutine PMSensitivityRichardsInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation.
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  implicit none

  class(pm_sensitivity_richards_type) :: this
  
  PetscErrorCode :: ierr
  
  ierr = 0
  
  ! create the J matrix ?
  ! create the output and the ij structure ?
  
end subroutine PMSensitivityRichardsInitializeRun

! *************************************************************************** !

subroutine PMSensitivityRichardsInitializeTimestep(this)
  ! 
  ! Initializes the process model to take a time step in the simulation.
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  implicit none

  class(pm_sensitivity_richards_type) :: this
  
  !call PMBasePrintHeader(this)
  
  !if (this%option%time >= this%output_start_time) then
  !  call PMUFDBOutput(this)
  !endif

  
end subroutine PMSensitivityRichardsInitializeTimestep

! *************************************************************************** !

 subroutine PMSensitivityRichardsSolve(this,time,ierr)
  ! 
  ! Main driver for outputting the richards sensitivities
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  use Option_module
  use Sensitivity_Analysis_module
  use Output_Sensitivity_module
  
  implicit none
  
  class(pm_sensitivity_richards_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  type(option_type), pointer :: option
  PetscBool :: to_output
  
  output_sensitivity_option => this%output_sensitivity_option
  option => this%realization%option
  to_output = PETSC_FALSE
  
  
  option%io_buffer = "debug: solve"
  call PrintMsg(option)
  
  option%io_buffer = "print output time option"
  call PrintMsg(option)
  write(option%io_buffer,'(i8)') output_sensitivity_option%output_time_option
  call PrintMsg(option)
  
  ! check if it's the time to output
  call OutputSensitivityOptionIsTimeToOutput(output_sensitivity_option,&
                                             to_output)
  
  option%io_buffer = "solve"
  call PrintMsg(option)
  ! compute the sensitivities at the output time
  if (to_output) then
    option%io_buffer = "compute sensitivity"
    call PrintMsg(option)
  endif
  
  ! output it
  if (to_output) then
    option%io_buffer = "output sensitivity"
    call PrintMsg(option)
  endif

end subroutine PMSensitivityRichardsSolve

! ************************************************************************** !

subroutine PMSensitivityRichardsFinalizeTimestep(this)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017

  implicit none
  
  class(pm_sensitivity_richards_type) :: this
  
end subroutine PMSensitivityRichardsFinalizeTimestep

! *************************************************************************** !

subroutine PMSensitivityRichardsOutput(this)
  ! 
  ! Sets up output for the process model.
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  use Option_module
  use Output_Sensitivity_module
  
  implicit none

  class(pm_sensitivity_richards_type) :: this
  
  type(option_type), pointer :: option
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  
  ! TODO
  
end subroutine PMSensitivityRichardsOutput

! *************************************************************************** !

subroutine PMSensitivityRichardsInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  ! 
  
  implicit none
  
  class(pm_sensitivity_richards_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT
  
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  
end subroutine PMSensitivityRichardsInputRecord

! *************************************************************************** !

subroutine PMSensitivityRichardsStrip(this)
  ! 
  ! Strips the Sensitivity Richards process model.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/04/2021
  !
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  class(pm_sensitivity_richards_type) :: this
  
  call PMBaseDestroy(this)

  nullify(this%realization)
  
  ! TODO 
  ! nullify output_sensitivity_option

end subroutine PMSensitivityRichardsStrip

! ************************************************************************** !

subroutine PMSensitivityRichardsDestroy(this)
  ! 
  ! Strips and destroys the UFD Biosphere process model.
  ! 
  ! Author: Sensitivity Richards
  ! Date: 01/04/2021
  !

  implicit none
  
  class(pm_sensitivity_richards_type) :: this
  
  call PMSensitivityRichardsStrip(this)
  
end subroutine PMSensitivityRichardsDestroy

! ************************************************************************** !

end module PM_Sensitivity_Richards_class
