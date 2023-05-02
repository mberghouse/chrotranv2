module Init_Surface_Flow_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none


  public :: InitSurfaceFlowSetupRealization

contains

! ************************************************************************** !

subroutine InitSurfaceFlowSetupRealization(simulation)
  !
  ! Initializes material property data structres and assign them to the domain.
  !
  ! Author: Gautam Bisht
  ! Date: 04/03/23
  !
  use Condition_Control_module
  use Option_module
  use Patch_module
  use Realization_Surface_class
  use Simulation_Surface_class
  use PM_Base_class
  use PM_SWE_class
  use SWE_module

  implicit none

  class(simulation_surface_type) :: simulation

  class(pm_base_type), pointer :: pm_list
  class(realization_surface_type), pointer :: surface_realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch

  surface_realization => simulation%surface_realization
  option => surface_realization%option
  patch => surface_realization%patch

  pm_list => simulation%surface_flow_process_model_coupler%pm_list
  if (associated(pm_list%next)) then
    option%io_buffer = 'Select type block in InitSurfaceFlowSetupRealization &
      &Realization needs to be refactored since there is more than one &
      &process model in simulation%flow_process_model_coupler%pm_list.'
    call PrintErrMsg(option)
  endif

  select type(pm => pm_list)
    class is (pm_swe_type)
      call SWESetup(surface_realization)
    class default
      option%io_buffer = 'Unknown surface flow mode found during setup'
      call PrintErrMsg(option)
  end select

  call CondControlAssignFlowInitCondSurface(surface_realization)
  !call InitSubsurfFlowReadInitCond()     =  Not implemented

  select case(option%iflowmode)
  case (SWE_MODE)
    call SWEUpdateAuxVars(surface_realization)
  case default
    option%io_buffer = 'Unknown flowmode found during <Mode>UpdateAuxVars'
    call PrintErrMsg(option)
end select

end subroutine InitSurfaceFlowSetupRealization

end module Init_Surface_Flow_module
