module Factory_Surface_linage_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Simulation_Surface_class
  use PFLOTRAN_Constants_module
  use PM_Surface_Flow_class
  use PM_SWE_class

  implicit none

  private

  public :: &
            FactorySurfaceLinkExtractPMsFromPMList, &
            FactorySurfaceLinkSetupPMCLinages

contains

! ************************************************************************** !

subroutine FactorySurfaceLinkExtractPMsFromPMList(simulation,pm_surface_flow)
  !
  ! Extracts all possible PMs from the PM list
  !
  ! Author: Gautam Bisht
  ! Date: 06/05/18
  !

  use PM_Surface_Flow_class
  use PM_Base_class
  use Simulation_Surface_class
  use Option_module

  implicit none

  class(simulation_surface_type) :: simulation

  type(option_type), pointer :: option
  class(pm_surface_flow_type), pointer :: pm_surface_flow
  class(pm_base_type), pointer :: cur_pm, prev_pm

  option => simulation%option

  nullify(pm_surface_flow)

  cur_pm => simulation%process_model_list
  do
    if (.not.associated(cur_pm)) exit
    select type(cur_pm)
      class is(pm_surface_flow_type)
        pm_surface_flow => cur_pm
      class default
        option%io_buffer = &
        'PM Class unrecognized in FactorySurfaceLinkExtractPMsFromPMList.'
        call PrintErrMsg(option)
  end select

  prev_pm => cur_pm
  cur_pm => cur_pm%next

  ! we must destroy the linkage between pms so that they are in independent
  ! lists among pmcs
  nullify(prev_pm%next)

  enddo

end subroutine FactorySurfaceLinkExtractPMsFromPMList

! ************************************************************************** !

subroutine FactorySurfaceLinkSetupPMCLinages(simulation,pm_surface_flow)

  !
  ! Sets up all PMC linkages
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !

  use Input_Aux_module
  use PM_Base_class
  use PMC_Surface_class
  use Option_module
  use Realization_Surface_class
  use Timestepper_TS_class
  use Factory_Surface_Read_module

  class(simulation_surface_type) :: simulation
  class(pm_surface_flow_type), pointer :: pm_surface_flow

  class(pmc_surface_type), pointer :: pmc_surface
  type(option_type), pointer :: option
  class(realization_surface_type), pointer :: surface_realization
  type(input_type), pointer :: input

  pmc_surface => PMCSurfaceCreate()

  call FactorySurfaceLinkAddPMCSurfaceFlow(simulation,pm_surface_flow,'PMCSurfaceFlow')

  surface_realization => simulation%surface_realization
  option => surface_realization%option

  call pmc_surface%SetName('PMCSurfaceFlow')
  call pmc_surface%SetOption(option)
  call pmc_surface%SetWaypointList(simulation%waypoint_list_surface)

  pmc_surface%pm_list => pm_surface_flow
  pmc_surface%pm_ptr%pm => pm_surface_flow
  pmc_surface%surface_realization => surface_realization

  select type(pm_surface_flow)
    class is (pm_swe_type)
      pmc_surface%timestepper => TimestepperTSCreate()
    class default
      option%io_buffer = 'Unsupported PM in FactorySurfaceLinkSetupPMCLinages'
      call PrintErrMsg(option)
  end select
  pmc_surface%timestepper%name = 'SURFACE_FLOW'

  call pmc_surface%pm_list%InitializeSolver()

  simulation%surface_flow_process_model_coupler => pmc_surface
  simulation%surface_flow_process_model_coupler_list => &
    simulation%surface_flow_process_model_coupler

  input => InputCreate(IN_UNIT,option%input_filename,option)
  call FactorySurfaceReadRequiredCards(simulation,input)
  call FactorySurfaceReadInput(simulation,input)

  end subroutine FactorySurfaceLinkSetupPMCLinages

! ************************************************************************** !

subroutine FactorySurfaceLinkAddPMCSurfaceFlow(simulation,pm_surface_flow,pmc_name)

  !
  ! Sets up all PMC linkages
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !

  use PM_Base_class
  use PMC_Surface_class
  use Option_module
  use Realization_Surface_class
  use Timestepper_TS_class

  class(simulation_surface_type) :: simulation
  class(pm_surface_flow_type), pointer :: pm_surface_flow
  character(len=*) :: pmc_name

  class(pmc_surface_type), pointer :: pmc_surface
  type(option_type), pointer :: option
  class(realization_surface_type), pointer :: surface_realization

  surface_realization => simulation%surface_realization
  option => surface_realization%option

  pmc_surface => PMCSurfaceCreate()

  call pmc_surface%SetName(pmc_name)
  call pmc_surface%SetOption(option)
  call pmc_surface%SetWaypointList(simulation%waypoint_list_surface)

  pmc_surface%pm_list => pm_surface_flow
  pmc_surface%pm_ptr%pm => pm_surface_flow
  pmc_surface%surface_realization => surface_realization

  select type(pm_surface_flow)
    class is (pm_swe_type)
      pmc_surface%timestepper => TimestepperTSCreate()
    class default
      option%io_buffer = 'Unsupported PM in FactorySurfaceLinkSetupPMCLinages'
      call PrintErrMsg(option)
  end select
  pmc_surface%timestepper%name = 'SURFACE_FLOW'

  call pmc_surface%pm_list%InitializeSolver()
  
end subroutine FactorySurfaceLinkAddPMCSurfaceFlow

end module Factory_Surface_linage_module