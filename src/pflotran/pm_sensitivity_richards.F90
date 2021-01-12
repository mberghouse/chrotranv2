module PM_Sensitivity_Richards_class

#include "petsc/finclude/petscsys.h"
  use petscsys
#include "petsc/finclude/petscsnes.h"
  use petscsnes
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
    !procedure, public :: Output => PMSensitivityRichardsOutput
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
  option%io_buffer = word
  call PrintMsg(option)

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
                                                  option,word,input)
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
  ! Just determine and update the output flag
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
  PetscBool :: timestep_flag, last_flag, sync_flag
  
  output_sensitivity_option => this%output_sensitivity_option
  option => this%realization%option
  timestep_flag = PETSC_FALSE
  last_flag = PETSC_FALSE
  sync_flag = PETSC_FALSE
  ierr = 0
  
  !update time for future calling of output
  output_sensitivity_option%time = time
  output_sensitivity_option%plot_number = &
                       output_sensitivity_option%plot_number + 1
  output_sensitivity_option%plot_flag = PETSC_FALSE
  
  !check for timestep
  if (floor(real(output_sensitivity_option%plot_number) / &
            output_sensitivity_option%output_every_timestep) /= 0) &
    timestep_flag = PETSC_TRUE
  !check for last
  ! TODO (moise)
  !check for sync
  ! TODO (moise)
  
  ! check if it's the time to output
  call OutputSensitivityOptionIsTimeToOutput(output_sensitivity_option,&
                                             timestep_flag,last_flag,sync_flag)
  
  call PMSensitivityRichardsOutput(this)

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
  use Sensitivity_Analysis_module
  use Discretization_module
  
  implicit none

  class(pm_sensitivity_richards_type) :: this
  
  type(option_type), pointer :: option
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  type(output_sensitivity_variable_type), pointer :: output_variable
  Mat :: J
  MatType :: J_mat_type
  PetscErrorCode :: ierr
  
  output_sensitivity_option => this%output_sensitivity_option
  option => this%realization%option
  
  if (output_sensitivity_option%plot_flag) then
    
    !prepare J matrix
    J_mat_type = MATBAIJ
    call DiscretizationCreateJacobian(this%realization%discretization, &
                                      option%nflowdof, J_mat_type, J, option)
    call MatSetOptionsPrefix(J,"Sensitivity_",ierr);CHKERRQ(ierr)
    
    output_variable => output_sensitivity_option%output_variables%first
    do 
      
      ! compute sensitivity using finite difference
      option%io_buffer = "SENSITIVITY ANALYSIS: Compute " // &
                         trim(output_variable%name) // " Sensitivity Matrix"
      call PrintMsg(option)
      
      call MatZeroEntries(J,ierr);CHKERRQ(ierr)
      call RichardsSensitivityInternalConn(J,this%realization,&
                                           output_variable%ivar,ierr)
      call RichardsSensitivityBoundaryConn(J,this%realization,&
                                           output_variable%ivar,ierr)
      !call RichardsSensitivitySourceSink(J,realization,ivar,ierr)
      !update here when porosity ok
      !call RichardsSensitivityAccumulation(J,realization,ivar,ierr)
      call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      
      !output it
      call OutputSensitivity(J,option,output_sensitivity_option,&
                             output_variable%name)
      
      if (.not.associated(output_variable%next)) exit
      output_variable => output_variable%next
    enddo
    
  endif 
  
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
