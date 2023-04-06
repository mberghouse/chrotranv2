module Factory_Surface_module

#include "petsc/finclude/petscsys.h"

use petscsys
  use Simulation_Surface_class

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal

  implicit none

  private

  public :: FactorySurfaceInitialize

contains

! ************************************************************************** !

subroutine FactorySurfaceInitialize(simulation)
  !
  ! Sets up PFLOTRAN surface simulation
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !

  implicit none

  class(simulation_surface_type) :: simulation

  ! NOTE: PETSc must already have been initialized here!
  call FactorySurfaceInitPostPetsc(simulation)

end subroutine FactorySurfaceInitialize

! ************************************************************************** !

subroutine FactorySurfaceInitPostPetsc(simulation)
  !
  ! Sets up PFLOTRAN surface simulation
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !
  use Option_module
  use PM_Base_class
  use PM_Surface_Flow_class
  use Realization_Surface_class
  use Factory_Surface_linage_module

  implicit none

  class(simulation_surface_type) :: simulation

  type(option_type), pointer :: option
  class(pm_surface_flow_type), pointer :: pm_surface_flow
  class(realization_surface_type), pointer :: realization_surface

  option => simulation%option

  nullify(pm_surface_flow)

  write(*,*)'  call FactorySurfaceLinkExtractPMsFromPMList'
  call FactorySurfaceLinkExtractPMsFromPMList(simulation,pm_surface_flow)

  write(*,*)'  call FactorySurfaceSetFlowMode'
  call FactorySurfaceSetFlowMode(pm_surface_flow,option)

  realization_surface => RealizationSurfaceCreate(option)

  simulation%surface_realization => realization_surface

  call FactorySurfaceLinkSetupPMCLinages(simulation,pm_surface_flow)

  call FactorySurfaceInitSimulation(simulation)

  ! set first process model coupler as the master
  simulation%process_model_coupler_list%is_master = PETSC_TRUE

end subroutine FactorySurfaceInitPostPetsc

! ************************************************************************** !

subroutine FactorySurfaceSetFlowMode(pm_surface_flow,option)
  !
  ! Sets the surface flow mode
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !
  use Option_module
  use PM_Surface_Flow_class
  use PM_Base_class
  use PM_SWE_class
  use Option_module

  implicit none

  type(option_type), pointer :: option
  class(pm_surface_flow_type), pointer :: pm_surface_flow

  option%liquid_phase = 1

  if (.not.associated(pm_surface_flow)) then
    option%nphase = 1
    ! assume default isothermal when only transport
    option%use_isothermal = PETSC_TRUE
    return
  endif

  select type(pm_surface_flow)
    class is (pm_swe_type)
      write(*,*)'  In FactorySurfaceSetFlowMode: pm_swe_type'
      option%iflowmode = SWE_MODE
      option%nphase = 1
      option%nflowdof = 3
      option%nflowspec = 1
      option%use_isothermal = PETSC_TRUE

    class default
      option%io_buffer = 'Unsupported pm_surface_flow type in FactorySurfaceSetFlowMode'
      call PrintErrMsg(option)
  end select

  if (option%nflowdof == 0) then
    option%io_buffer = 'Number of flow degrees of freedom is zero.'
    call PrintErrMsg(option)
  endif
  if (option%nphase == 0) then
    option%io_buffer = 'Number of flow phases is zero.'
    call PrintErrMsg(option)
  endif
  if (option%nflowspec == 0) then
    option%io_buffer = 'Number of flow species is zero.'
    call PrintErrMsg(option)
  endif

end subroutine FactorySurfaceSetFlowMode

! ************************************************************************** !

subroutine FactorySurfaceInitSimulation(simulation)
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !

  use Init_Common_module
  use Init_Surface_Flow_module
  use Option_module
  use Realization_Surface_class

  implicit none

  class(simulation_surface_type) :: simulation

  class(realization_surface_type), pointer :: realization_surface
  type(option_type), pointer :: option

  write(*,*)'Add code in FactorySurfaceInitSimulation'
  call FactorySurfaceSetupRealization(simulation)

  realization_surface => simulation%surface_realization
  option => realization_surface%option
  call InitCommonAddOutputWaypoints(option,simulation%output_option, &
                                    simulation%waypoint_list_surface)

  call InitSurfaceFlowSetupRealization(simulation)
  call exit(0)

end subroutine FactorySurfaceInitSimulation

! ************************************************************************** !

subroutine FactorySurfaceSetupRealization(simulation)
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !
  use EOS_module
  use Option_module
  use Init_Common_module
  use Realization_Surface_class
  use Realization_Subsurface_class
  use Waypoint_module
  use Init_Surface_Flow_module
  use Init_Surface_module

  implicit none

  class(simulation_surface_type) :: simulation

  class(realization_surface_type), pointer :: realization_surface
  type(waypoint_list_type) :: waypoint_list
  
  type(option_type), pointer :: option

  realization_surface => simulation%surface_realization
  option => realization_surface%option

  ! set reference densities if not specified in input file.
  call EOSReferenceDensity(option)

  call RealizationSurfaceCreateDiscretization(realization_surface)
  realization_surface%discretization%grid%unstructured_grid%grid_type = TWO_DIM_GRID

  call InitCommonReadRegionFiles(realization_surface%patch,realization_surface%surf_region_list, &
                                 option)

  call RealizationLocalizeRegions(realization_surface%patch,realization_surface%surf_region_list, &
                                  option)

  call RealizationSurfacePassPtrsToPatches(realization_surface)
  call RealizationSurfaceProcessMatProp(realization_surface)
  call RealizationSurfaceProcessConditions(realization_surface)
  call RealizationSurfaceProcessCouplers(realization_surface)
  call SurfaceInitMaterialProperties(realization_surface)

  write(*,*)'Stopping in FactorySurfaceSetupRealization'
  call exit(0)

end subroutine FactorySurfaceSetupRealization

end module Factory_Surface_module