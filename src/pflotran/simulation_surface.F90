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
    procedure, public :: InitializeRun => SimSurfaceInitializeRun
    procedure, public :: ExecuteRun => SimSurfaceExecuteRun
    procedure, public :: RunToTime => SimSurfaceRunToTime
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

subroutine SimSurfaceInitializeRun(this)
  !
  ! Initializes simulation
  !
  ! Author: Gautam Bisht
  ! Date: 05/09/23
  !
  use Timestepper_Base_class, only : TS_CONTINUE
  use Waypoint_module
  use Driver_class
  use Option_module

  implicit none

  class(simulation_surface_type) :: this

#ifdef DEBUG
  call PrintMsg(this%option,'SimSurfaceInitializeRun()')
#endif

  call SimulationBaseInitializeRun(this)

  call this%surface_flow_process_model_coupler_list%InitializeRun()

end subroutine SimSurfaceInitializeRun

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
  nullify(this%surface_flow_process_model_coupler_list)
  nullify(this%process_model_list)
  this%sim_aux => SimAuxCreate()
  this%stop_flag = TS_CONTINUE
  nullify(this%surface_flow_process_model_coupler)
  nullify(this%surface_realization)
  nullify(this%regression)
  this%waypoint_list_surface => WaypointListCreate()
  this%waypoint_list_outer => WaypointListCreate()

end subroutine SimSurfaceInit

! ************************************************************************** !

subroutine SimSurfaceExecuteRun(this)
  !
  ! Execute a simulation
  !
  ! Author: Gautam Bisht
  ! Date: 05/02/23
  !
  use Waypoint_module
  use Timestepper_Base_class, only : TS_CONTINUE
  use Checkpoint_module

  implicit none

  class(simulation_surface_type) :: this

  PetscReal :: final_time
  type(waypoint_type), pointer :: cur_waypoint
  character(len=MAXSTRINGLENGTH) :: append_name

#ifdef DEBUG
  call PrintMsg(this%option,'SimSurfaceExecuteRun()')
#endif

  if (.not.associated(this%surface_flow_process_model_coupler_list)) then
    return
  endif

  final_time = SimSurfaceGetFinalWaypointTime(this)
  cur_waypoint => this%waypoint_list_outer%first
  if (cur_waypoint%print_checkpoint) then
    append_name = &
         CheckpointAppendNameAtTime(this%surface_flow_process_model_coupler_list% &
                                        option%time, &
                                    this%surface_flow_process_model_coupler_list%option)
    call this%surface_flow_process_model_coupler_list%Checkpoint(append_name)
  endif
  call WaypointSkipToTime(cur_waypoint,this%option%time)
  do
    if (this%stop_flag /= TS_CONTINUE) exit ! end simulation
    if (.not.associated(cur_waypoint)) exit
    call this%RunToTime(min(final_time,cur_waypoint%time))
    cur_waypoint => cur_waypoint%next
  enddo
  append_name = '-restart'
  if (associated(this%option%checkpoint)) then
    call this%surface_flow_process_model_coupler_list%Checkpoint(append_name)
  endif

end subroutine SimSurfaceExecuteRun

! ************************************************************************** !

function SimSurfaceGetFinalWaypointTime(this)
  !
  ! Returns the earliest final waypoint time from the top layer of process
  ! model couplers.
  !
  ! Author: Gautam Bisht
  ! Date: 05/03/23
  !
  use Waypoint_module

  implicit none

  class(simulation_surface_type) :: this

  PetscReal :: SimSurfaceGetFinalWaypointTime

  class(pmc_base_type), pointer :: cur_process_model_coupler
  PetscReal :: final_time

  SimSurfaceGetFinalWaypointTime = 0.d0

  cur_process_model_coupler => this%surface_flow_process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    final_time = WaypointListGetFinalTime(cur_process_model_coupler% &
                                            waypoint_list)
    if (SimSurfaceGetFinalWaypointTime < 1.d-40 .or. &
        final_time < SimSurfaceGetFinalWaypointTime) then
          SimSurfaceGetFinalWaypointTime = final_time
    endif
    cur_process_model_coupler => cur_process_model_coupler%peer
  enddo

end function SimSurfaceGetFinalWaypointTime

! ************************************************************************** !

subroutine SimSurfaceRunToTime(this,target_time)
  !
  ! Executes simulation
  !
  ! Author: Gautam Bisht
  ! Date: 05/03/23
  !
  use Option_module
  use Simulation_Aux_module

  implicit none

  class(simulation_surface_type) :: this
  PetscReal :: target_time

#ifdef DEBUG
  call PrintMsg(this%option,'SimSurfaceRunToTime()')
#endif

  call this%surface_flow_process_model_coupler_list%RunToTime(target_time,this%stop_flag)

end subroutine SimSurfaceRunToTime

end module Simulation_Surface_class
