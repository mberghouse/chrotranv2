module Simulation_Surface_class
    
#include "petsc/finclude/petscvec.h"
use petscvec
use PFLOTRAN_Constants_module
use Simulation_Base_class
use Simulation_Aux_module
use Option_module
use Output_Aux_module
use PM_Base_class
use PMC_Base_class
use PMC_Surface_class
use Waypoint_module
use Regression_module
use Realization_Surface_class

  implicit none

  private

  type, public, extends(simulation_base_type) :: simulation_surface_type
    type(option_type), pointer :: option    
    PetscInt :: stop_flag
    type(output_option_type), pointer :: output_option
    class(pmc_base_type), pointer :: process_model_coupler_list
    class(pm_base_type), pointer :: process_model_list
    type(simulation_aux_type), pointer :: sim_aux
    class(pmc_base_type), pointer :: surface_flow_process_model_coupler_list
    class(pmc_surface_type), pointer :: surface_flow_process_model_coupler
    class(realization_surface_type), pointer :: surface_realization
    ! regression object
    type(regression_type), pointer :: regression
    type(waypoint_list_type), pointer :: waypoint_list_surface
    type(waypoint_list_type), pointer :: waypoint_list_outer ! outer sync loop
  contains
  end type simulation_surface_type

  public :: SimSurfaceCreate

contains

! ************************************************************************** !

function SimSurfaceCreate(driver,option)
  !
  ! Allocates and initializes a new simulation object
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !
  use Driver_class
  use Option_module

  implicit none

  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  class(simulation_surface_type), pointer :: SimSurfaceCreate

#ifdef DEBUG
  print *, 'SimSurfaceCreate'
#endif

  allocate(SimSurfaceCreate)
  call SimSurfaceInit(SimSurfaceCreate,driver,option)

end function SimSurfaceCreate

! ************************************************************************** !

subroutine SimSurfaceInit(this,driver,option)
  !
  ! Initializes simulation values
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !
  use Timestepper_Base_class, only : TS_CONTINUE
  use Waypoint_module
  use Driver_class
  use Option_module

  implicit none

  class(simulation_surface_type) :: this
  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

#ifdef DEBUG
  call PrintMsg(this%option,'SimSurfaceInit()')
#endif

  call SimulationBaseInit(this,driver)
  this%option => option
  this%output_option => OutputOptionCreate()
  nullify(this%process_model_coupler_list)
  nullify(this%process_model_list)
  this%sim_aux => SimAuxCreate()
  this%stop_flag = TS_CONTINUE
  nullify(this%surface_flow_process_model_coupler)
  nullify(this%surface_realization)
  nullify(this%regression)
  this%waypoint_list_surface => WaypointListCreate()
  this%waypoint_list_outer => WaypointListCreate()

end subroutine SimSurfaceInit

end module Simulation_Surface_class