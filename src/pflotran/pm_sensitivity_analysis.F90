module PM_Sensitivity_Analysis_class

#include "petsc/finclude/petscsys.h"
  use petscsys
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module
  use Sensitivity_Aux_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_sensitivity_type
    class(realization_subsurface_type), pointer :: realization
    class(sensitivity_output_option_type), pointer :: sensitivity_output_option
    PetscBool :: sensitivity_flow !indicate to compute sensitivity for flow
    PetscBool :: sensitivity_transport
  contains
    procedure, public :: SetRealization => PMSensitivitySetRealization
    procedure, public :: Setup => PMSensitivitySetup
    procedure, public :: ReadPMBlock => PMSensitivityReadPMBlock
    procedure, public :: ReadSimulationOptionsBlock => &
                                         PMSensitivitysReadSimOptionsBlock
    procedure, public :: InitializeRun => PMSensitivityInitializeRun
    procedure, public :: InitializeTimestep => PMSensitivityInitializeTimestep
    procedure, public :: FinalizeTimestep => PMSensitivityFinalizeTimestep
    procedure, public :: Solve => PMSensitivitySolve
    !procedure, public :: Output => PMSensitivityOutput
    procedure, public :: InputRecord => PMSensitivityInputRecord
    procedure, public :: Destroy => PMSensitivityDestroy
  end type pm_sensitivity_type

  public :: PMSensitivityCreate
  
contains

! *************************************************************************** !

function PMSensitivityCreate()
  !
  ! Creates and initializes the sensitivity Richards process model.
  !
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  implicit none
  
  class(pm_sensitivity_type), pointer :: PMSensitivityCreate
  
  allocate(PMSensitivityCreate)
  call PMBaseInit(PMSensitivityCreate)
  nullify(PMSensitivityCreate%realization)
  PMSensitivityCreate%name = 'sensitivity'
  PMSensitivityCreate%header = 'SENSITIVITY ANALYSIS'
  PMSensitivityCreate%sensitivity_flow = PETSC_FALSE
  PMSensitivityCreate%sensitivity_transport = PETSC_FALSE
  allocate(PMSensitivityCreate%sensitivity_output_option)
  call SensitivityOutputOptionInit( &
                         PMSensitivityCreate%sensitivity_output_option)
  
end function PMSensitivityCreate

! *************************************************************************** !

subroutine PMSensitivitySetRealization(this,realization)
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !

  use Realization_Subsurface_class

  implicit none
  
  class(pm_sensitivity_type) :: this
  class(realization_subsurface_type), pointer :: realization
  
  this%realization => realization
  this%realization_base => realization

end subroutine PMSensitivitySetRealization

! *************************************************************************** !

subroutine PMSensitivitysReadSimOptionsBlock(this,input)
  ! 
  ! Reads input file parameters associated with the Richards process model
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/12/2021

  use Input_Aux_module
  use String_module
  use Option_module
  use Sensitivity_Aux_module
 
  implicit none
  
  class(pm_sensitivity_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found
  PetscReal :: tempreal

  option => this%option
  error_string = 'Sensitivity Analysis Options'
  input%ierr = 0
  word = ''
  
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('SENSITIVITY_FLOW')
        this%sensitivity_flow = PETSC_TRUE
      case('SENSITIVITY_TRANSPORT')
        this%sensitivity_transport = PETSC_TRUE
        option%io_buffer = 'SENSITIVITY_ANALYSIS for TRANSPORT process ' // &
                           'not yet implemented'
        call PrintErrMsg(option)
      case default
        option%io_buffer = 'Keyword ' // trim(word) // &
              ' not recognized for the ' // trim(error_string) // ' block.'
        call PrintErrMsg(option)
    end select
  enddo
  
  call InputPopBlock(input,option)
  
end subroutine PMSensitivitysReadSimOptionsBlock

! *************************************************************************** !

subroutine PMSensitivityReadPMBlock(this,input)
  !
  ! Reads input file parameters for the sensitivity flow
  !
  ! Author: Moise Rousseau
  ! Date: 03/13/2017
  !
#if 0
SENSTIVITITY_FLOW
  FORMAT HDF5
  VARIABLES
    PRESSURE
    PERMEABILITY
    POROSITY
  /
  OUTPUT SYNC_WITH_SNAPSHOT_FILE
END
#endif
  
  use Input_Aux_module
  use Option_module
  use String_module
  use Sensitivity_Aux_module
  
  implicit none
  
  class(pm_sensitivity_type) :: this
  type(input_type), pointer :: input
  
  class(sensitivity_output_option_type), pointer :: sensitivity_output_option
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: ierr
  
  sensitivity_output_option => this%sensitivity_output_option
  option => this%option
  input%ierr = 0
  
  option%io_buffer = 'pflotran card:: SENSITIVITY FLOW'
  call PrintMsg(option)
  
  error_string = 'sensitivity'
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
        select case(trim(word))
          case('HDF5')
            sensitivity_output_option%format = SENSITIVITY_OUTPUT_HDF5
          case('MATLAB')
            sensitivity_output_option%format = SENSITIVITY_OUTPUT_MATLAB
          case('BINARY')
            sensitivity_output_option%format = SENSITIVITY_OUTPUT_BINARY
          case('ASCII')
            sensitivity_output_option%format = SENSITIVITY_OUTPUT_ASCII
          case default
            call InputKeywordUnrecognized(input,word,'SENSITIVITY_FLOW,&
                                          &FORMAT',option)
        end select
    !-----------------------------------------
      case('VARIABLES')
        call PMSensitivityReadOutputVariables(this, &
                                              sensitivity_output_option,input)
    !-----------------------------------------
      case('OUTPUT')
        call InputReadCard(input,option,word)
        call StringToUpper(word)
        select case(trim(word))
          case('SYNC_WITH_SNAPSHOT_FILE')
            sensitivity_output_option%output_time_option = &
                                                SENSITIVITY_SYNC_OUTPUT
          case('PERIODIC_TIMESTEP')
            sensitivity_output_option%output_time_option = &
                                                SENSITIVITY_EVERY_X_TIMESTEP
            call InputReadInt(input,option,&
                              sensitivity_output_option%output_every_timestep)
          case('LAST_TIMESTEP')
            sensitivity_output_option%output_time_option = SENSITIVITY_LAST
          case default
            call InputKeywordUnrecognized(input,word,'SENSITIVITY_FLOW,&
                                          &OUTPUT_TIME',option)
        end select
    !-----------------------------------------
      case default
        call InputKeywordUnrecognized(input,word,'SENSITIVITY_FLOW',option)
    !-----------------------------------------
    end select  
  enddo
  call InputPopBlock(input,option)
  
end subroutine PMSensitivityReadPMBlock

! *************************************************************************** !

subroutine PMSensitivityReadOutputVariables(this, &
                                            sensitivity_output_option,input)
  !
  ! Sets up the process model with external information.
  !
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !

  use Input_Aux_module
  use Option_module
  use String_module
  use Sensitivity_Aux_module
  
  implicit none
  
  class(pm_sensitivity_type) :: this
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  type(input_type), pointer :: input
  
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  type(sensitivity_output_variable_type), pointer :: variable
  
  option => this%option
  
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'sensitivity, VARIABLES')
    call StringToUpper(word)
    
    select case(trim(word))
      case('PRESSURE')
        variable => SensitivityOutputVariableCreate()
        variable%name = 'Pressure'
        variable%units = '1.Pa-1' ! TODO (moise)
        variable%ivar = SENSITIVITY_PRESSURE
      case('PERMEABILITY')
        variable => SensitivityOutputVariableCreate()
        variable%name = 'Permeability'
        variable%units = '1.s-1' ! TODO (moise)
        variable%ivar = SENSITIVITY_PERMEABILITY
      case('POROSITY')
        variable => SensitivityOutputVariableCreate()
        variable%name = 'Porosity'
        variable%units = '' ! TODO (moise)
        variable%ivar = SENSITIVITY_POROSITY
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'SENSITIVITY_FLOW, VARIABLES', &
                                      option)
    end select
    call SensitivityAddOutputVariableToList(sensitivity_output_option, &
                                            variable)
  end do
  
  call InputPopBlock(input,option)
  
end subroutine PMSensitivityReadOutputVariables

! *************************************************************************** !

subroutine PMSensitivitySetup(this)
  !
  ! Sets up the process model with external information.
  ! i.e. check if the correct flow mode is selected
  !
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  use Option_module

  implicit none
  
  class(pm_sensitivity_type) :: this
  type(option_type), pointer :: option
  
  option => this%option
  
  select case(option%iflowmode)
    case(RICHARDS_MODE)
    case default
      option%io_buffer = "Flow mode not supported for carrying a Sensitivity &
                          &Analysis."
      call PrintErrMsg(option)
  end select
  
end subroutine PMSensitivitySetup

! *************************************************************************** !

subroutine PMSensitivityInitializeRun(this)
  ! 
  ! Initializes the process model for the simulation.
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  implicit none

  class(pm_sensitivity_type) :: this
  
  PetscErrorCode :: ierr
  
  ierr = 0
  
  ! create the J matrix ?
  ! create the output and the ij structure ?
  
end subroutine PMSensitivityInitializeRun

! *************************************************************************** !

subroutine PMSensitivityInitializeTimestep(this)
  ! 
  ! Initializes the process model to take a time step in the simulation.
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  implicit none

  class(pm_sensitivity_type) :: this
  
  !call PMBasePrintHeader(this)
  
  !if (this%option%time >= this%output_start_time) then
  !  call PMUFDBOutput(this)
  !endif

  
end subroutine PMSensitivityInitializeTimestep

! *************************************************************************** !

 subroutine PMSensitivitySolve(this,time,ierr)
  ! 
  ! Just determine and update the output flag
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  use Option_module
  use Sensitivity_Output_module
  use Sensitivity_Aux_module
  
  implicit none
  
  class(pm_sensitivity_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  type(option_type), pointer :: option
  PetscBool :: timestep_flag, last_flag, sync_flag
  
  sensitivity_output_option => this%sensitivity_output_option
  option => this%realization%option
  timestep_flag = PETSC_FALSE
  last_flag = PETSC_FALSE
  sync_flag = PETSC_FALSE
  ierr = 0
  
  !update time for future calling of output
  sensitivity_output_option%time = time
  sensitivity_output_option%plot_number = &
                       sensitivity_output_option%plot_number + 1
  sensitivity_output_option%plot_flag = PETSC_FALSE
  
  !check for timestep
  if (floor(real(sensitivity_output_option%plot_number) / &
            sensitivity_output_option%output_every_timestep) /= 0) &
    timestep_flag = PETSC_TRUE
  !check for last
  ! TODO (moise)
  !check for sync
  ! TODO (moise)
  
  ! check if it's the time to output
  call SensitivityOutputOptionIsTimeToOutput(sensitivity_output_option,&
                                             timestep_flag,last_flag,sync_flag)
  
  call PMSensitivityOutput(this)

end subroutine PMSensitivitySolve

! ************************************************************************** !

subroutine PMSensitivityFinalizeTimestep(this)
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017

  implicit none
  
  class(pm_sensitivity_type) :: this
  
end subroutine PMSensitivityFinalizeTimestep

! *************************************************************************** !

subroutine PMSensitivityOutput(this)
  ! 
  ! Sets up output for the process model.
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  !
  
  use Option_module
  use Output_Aux_module
  use Discretization_module
  use Sensitivity_Richards_module
  !use Sensitivity_TH_module
  !use Sensitivity_Transport_module
  use Sensitivity_Output_module
  
  implicit none

  class(pm_sensitivity_type) :: this
  
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  type(sensitivity_output_variable_type), pointer :: output_variable
  Mat :: J
  MatType :: J_mat_type
  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: word
  
  output_option => this%realization%output_option
  sensitivity_output_option => this%sensitivity_output_option
  option => this%realization%option
  
  if (sensitivity_output_option%plot_flag) then
    
    !prepare J matrix
    J_mat_type = MATBAIJ
    call DiscretizationCreateJacobian(this%realization%discretization, &
                                      option%nflowdof, J_mat_type, J, option)
    call MatSetOptionsPrefix(J,"Sensitivity_",ierr);CHKERRQ(ierr)
    
    output_variable => sensitivity_output_option%output_variables%first
    do 
      
      if (this%sensitivity_flow) then
        ! compute sensitivity using finite difference
        option%io_buffer = trim(this%header) //": Compute " // &
	                   trim(output_variable%name) // " Flow &
	                   &Sensitivity Matrix"
        call PrintMsg(option)
        
        call MatZeroEntries(J,ierr);CHKERRQ(ierr)
        select case(option%iflowmode)
	        case(RICHARDS_MODE)
	          call RichardsSensitivity(J,this%realization,output_variable%ivar)
	        !case(TH)
	        !  call THSensitivity()
	        case default
        end select
        call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
        call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
        
        !output it
        word = "flow"
        call OutputSensitivity(J,option,output_option, &
	                             sensitivity_output_option, &
	                             output_variable, word)
	    endif
	    
	    if (this%sensitivity_transport) then
	      ! TODO (moise)
	      option%io_buffer = trim(this%header) //": Compute " // &
	                   trim(output_variable%name) // " Transport &&
	                   Sensitivity Matrix"
        call PrintMsg(option)
	    endif
      
      output_variable => output_variable%next
      if (.not.associated(output_variable)) exit
    enddo
    
  endif 
  
  option%io_buffer = "END " // trim(this%header) // achar(10)
  call PrintMsg(option)
  
end subroutine PMSensitivityOutput

! *************************************************************************** !

subroutine PMSensitivityInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/13/2017
  ! 
  
  implicit none
  
  class(pm_sensitivity_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT
  
  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

  
end subroutine PMSensitivityInputRecord

! *************************************************************************** !

subroutine PMSensitivityStrip(this)
  ! 
  ! Strips the Sensitivity process model.
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/12/2021
  !
  
  implicit none
  
  class(pm_sensitivity_type) :: this
  
  call PMBaseDestroy(this)
  nullify(this%realization)
  call SensitivityOutputOptionDestroy(this%sensitivity_output_option)
  nullify(this%sensitivity_output_option)

end subroutine PMSensitivityStrip

! ************************************************************************** !

subroutine PMSensitivityDestroy(this)
  ! 
  ! Strips and destroys the UFD Biosphere process model.
  ! 
  ! Author: Sensitivity Richards
  ! Date: 01/04/2021
  !

  implicit none
  
  class(pm_sensitivity_type) :: this
  
  call PMSensitivityStrip(this)
  
end subroutine PMSensitivityDestroy

! ************************************************************************** !

end module PM_Sensitivity_Analysis_class
