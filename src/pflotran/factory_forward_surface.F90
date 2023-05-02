module Factory_Forward_Surface_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Simulation_Surface_class

  implicit none

  private

  public :: FactoryForwardSurfaceInitialize

contains

! ************************************************************************** !

subroutine FactoryForwardSurfaceInitialize(simulation,input_filename,option)
  !
  ! Sets up the forward surface simulation framework after PETSc initialization
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23

  use Driver_class
  use Option_module
  use Print_module
  use Logging_module
  use Input_Aux_module
  use String_module
  use Factory_Forward_module, only : FactoryForwardReadCommandLine

  implicit none

  class(simulation_surface_type), intent(in), pointer :: simulation
  character(len=MAXSTRINGLENGTH) :: input_filename
  type(option_type), pointer :: option

  class(driver_type), pointer :: driver
  character(len=MAXSTRINGLENGTH) :: filename
  PetscErrorCode :: ierr

  driver => option%driver

  option%input_filename = input_filename
  option%global_prefix = StringStripFilenameSuffix(input_filename)

  call FactoryForwardReadCommandLine(option)

  ! popped in SimulationBaseInitializeRun()
  call PetscLogStagePush(logging%stage(INIT_STAGE),ierr);CHKERRQ(ierr)
  call PetscLogEventBegin(logging%event_init,ierr);CHKERRQ(ierr)

  filename = trim(option%global_prefix) // trim(option%group_prefix) // '.out'
  if (OptionPrintToFile(option)) then
    if (option%fid_out <= 0) option%fid_out = FORWARD_OUT_UNIT
    open(option%fid_out, file=filename, action="write", status="unknown")
  endif

  call OptionPrintPFLOTRANHeader(option)
  call FactoryForwardSurfaceReadSimulationBlk(simulation,driver,option)
  !if (.not.associated(option%inversion)) then
  !  call FactoryForwardPrerequisite(simulation)
  !endif
  call InputCheckKeywordBlockCount(option)

  call LoggingSetupComplete()

end subroutine FactoryForwardSurfaceInitialize

! ************************************************************************** !

subroutine FactoryForwardSurfaceReadSimulationBlk(simulation,driver,option)
  !
  ! Sets up a forward surface flow simulation framework after PETSc initialization
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !
  use Driver_class
  use Option_module
  use Input_Aux_module
  use String_module

  use PM_Base_class
  use PMC_Base_class
  use Checkpoint_module
  use Output_Aux_module
  use Option_Checkpoint_module
  use Waypoint_module
  use Units_module
  use Factory_Surface_Read_module
  use Factory_Surface_module
  use PM_Surface_Flow_class

  implicit none

  class(simulation_surface_type), pointer :: simulation
  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: simulation_type

  character(len=MAXSTRINGLENGTH) :: prerequisite_filename ! sim_base%prereq...
  class(pm_base_type), pointer :: pm_master
  class(pm_base_type), pointer :: cur_pm
  type(checkpoint_option_type), pointer :: checkpoint_option
  type(waypoint_list_type), pointer :: checkpoint_waypoint_list

  class(pmc_base_type), pointer :: pmc_master

  PetscBool :: print_ekg

  nullify(pm_master)
  nullify(cur_pm)

  nullify(pmc_master)
  nullify(checkpoint_option)
  nullify(checkpoint_waypoint_list)
  print_ekg = PETSC_FALSE

  write(*,*)'In FactoryForwardSurfaceReadSimulationBlk'
  input => InputCreate(IN_UNIT,option%input_filename,option)

  simulation_type = ''
  prerequisite_filename = ''
  string = 'SIMULATION'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'PROCESS_MODEL','SIMULATION')

    call StringToUpper(word)

    select case(trim(word))
      case('SIMULATION_TYPE')
          call InputReadCard(input,option,simulation_type,PETSC_FALSE)
          call InputErrorMsg(input,option,'simulation_type', &
                             'SIMULATION')
      case('PROCESS_MODELS')
        call FactoryForwardSurfaceReadSimProcessModels(input,pm_master,option)
      case('CHECKPOINT')
        option%io_buffer = 'CHECKPOINT for surface flow has not been implemented.'
        call PrintErrMsg(option)
      case ('RESTART')
        option%io_buffer = 'RESTART for surface flow has not been implemented.'
        call PrintErrMsg(option)
      case('INPUT_RECORD_FILE')
        option%input_record = PETSC_TRUE
        call OpenAndWriteInputRecord(option)
      case default
        call InputKeywordUnrecognized(input,word,'SIMULATION',option)
    end select
  enddo
  call InputPopBlock(input,option)
  call InputDestroy(input)

  if (.not.associated(pm_master)) then
    option%io_buffer = 'No process models defined in SIMULATION block.'
    call PrintErrMsg(option)
  endif

  if (.not.associated(simulation)) then
    ! create the simulation objects
    select case(simulation_type)
      case('SURFACE')
        simulation => SimSurfaceCreate(driver,option)
      case default
        if (len_trim(simulation_type) == 0) then
          option%io_buffer = 'A SIMULATION_TYPE (e.g. "SIMULATION_TYPE &
            &SURFACE") must be specified within the SIMULATION block.'
          call PrintErrMsg(option)
        endif
        call InputKeywordUnrecognized(input,simulation_type, &
                       'SIMULATION,SIMULATION_TYPE',option)
    end select
  endif

  call WaypointListMerge(simulation%waypoint_list_outer, &
                         checkpoint_waypoint_list,option)

  simulation%process_model_list => pm_master

  select type(simulation)
    class is(simulation_surface_type)
      call FactorySurfaceInitialize(simulation)
  end select

end subroutine FactoryForwardSurfaceReadSimulationBlk

! ************************************************************************** !

subroutine FactoryForwardSurfaceReadSimProcessModels(input,pm_master,option)
  !
  ! Reads in the process models listed in simulation block
  !
  !
    use Option_module
    use Input_Aux_module
    use String_module
  
    use PM_Base_class
    use PM_Auxiliary_class
  
    use Factory_Surface_Read_module
  
    implicit none
  
    class(pm_base_type), pointer :: pm_master
    type(input_type), pointer :: input
    type(option_type), pointer :: option
  
    character(len=MAXWORDLENGTH) :: word
    character(len=MAXWORDLENGTH) :: pm_name
    class(pm_base_type), pointer :: cur_pm
    class(pm_base_type), pointer :: new_pm
  
    nullify(cur_pm)
    nullify(new_pm)
  
    write(*,*)'  FactoryForwardSurfaceReadSimProcessModels'
    call InputPushBlock(input,option)
    do
      call InputReadPflotranString(input,option)
      if (InputCheckExit(input,option)) exit
      call InputReadCard(input,option,word)
      call InputErrorMsg(input,option,'process_model', &
                         'SIMULATION,PROCESS_MODELS')
      call InputReadWord(input,option,pm_name,PETSC_TRUE)
      if (InputError(input)) then
        input%err_buf = 'Process Model Name'
        call InputDefaultMsg(input,option)
        pm_name = ''
      endif
      call StringToUpper(word)
      select case(trim(word))
        case('SURFACE_FLOW')
          call FactorySurfaceReadFlowPM(input,option,new_pm)
        case default
          call InputKeywordUnrecognized(input,word, &
                 'SIMULATION,PROCESS_MODELS',option)
      end select
      if (.not.associated(new_pm%option)) new_pm%option => option
      if (len_trim(pm_name) > 0) then
        new_pm%name = pm_name
      endif
      if (associated(cur_pm)) then
        cur_pm%next => new_pm
      else
        cur_pm => new_pm
      endif
      write(*,*)'associated(pm_master) ?',associated(pm_master)
      if (.not.associated(pm_master)) then
        pm_master => new_pm
      endif
      cur_pm => new_pm
      nullify(new_pm)
    enddo
    call InputPopBlock(input,option)
  
  end subroutine FactoryForwardSurfaceReadSimProcessModels
  
end module Factory_Forward_Surface_module