module PM_Surface_Flow_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PM_Base_class
  use Realization_Surface_class
  use Communicator_Base_class

  implicit none

  private

  type, public, extends (pm_base_type) :: pm_surface_flow_type
    class(realization_surface_type), pointer :: surface_realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Setup => PMSurfaceFlowStep
    procedure, public :: SetRealization => PMSurfaceSetRealization
    procedure, public :: UpdateSolution => PMSurfaceUpdateSolution
    procedure, public :: UpdateTimestep => PMSurfaceUpdateTimestep
    procedure, public :: InitializeRun => PMSurfaceInitializeRun
  end type

  public :: PMSurfaceFlowInit

contains

! ************************************************************************** !

subroutine PMSurfaceFlowInit(this)

  ! Intializes shared members of subsurface process models
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23

  implicit none

  class(pm_surface_flow_type) :: this

  call PMBaseInit(this)
  nullify(this%surface_realization)
  nullify(this%comm1)

end subroutine PMSurfaceFlowInit

! ************************************************************************** !

subroutine PMSurfaceFlowStep(this)

  ! Sets realization_surface_type 
  !
  ! Author: Gautam Bisht
  ! Date: 05/03/23

  use Realization_Subsurface_class

  implicit none

  class(pm_surface_flow_type) :: this

  this%comm1 => this%surface_realization%comm1

end subroutine PMSurfaceFlowStep


! ************************************************************************** !

subroutine PMSurfaceSetRealization(this,surface_realization)

  ! Sets realization_surface_type 
  !
  ! Author: Gautam Bisht
  ! Date: 05/02/23

  use Realization_Subsurface_class

  implicit none

  class(pm_surface_flow_type) :: this
  class(realization_surface_type), pointer :: surface_realization

  this%surface_realization => surface_realization
  this%realization_base => surface_realization

  this%solution_vec = surface_realization%field_surface%flow_xx
  this%residual_vec = surface_realization%field_surface%flow_r

end subroutine PMSurfaceSetRealization

! ************************************************************************** !

subroutine PMSurfaceUpdateSolution(this)
  !
  ! Updates the solution
  !
  ! Author: Gautam Bisht
  ! Date: 05/08/23
  !
  use petscvec

  implicit none

  class(pm_surface_flow_type) :: this

  PetscErrorCode :: ierr

  call VecCopy(this%surface_realization%field_surface%flow_xx, &
               this%surface_realization%field_surface%flow_yy,ierr);CHKERRQ(ierr)

end subroutine PMSurfaceUpdateSolution

! ************************************************************************** !
subroutine PMSurfaceUpdateTimestep(this,update_dt, &
                                dt,dt_min,dt_max,iacceleration, &
                                num_newton_iterations,tfac, &
                                time_step_max_growth_factor)
  !
  ! Updates the solution
  !
  ! Author: Gautam Bisht
  ! Date: 05/08/23
  !
  implicit none
  !
  class(pm_surface_flow_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

end subroutine PMSurfaceUpdateTimestep

! ************************************************************************** !

subroutine PMSurfaceInitializeRun(this)
  !
  ! Initializes the PM
  !
  ! Author: Gautam Bisht
  ! Date: 05/08/23
  !
  use petscvec

  implicit none

  class(pm_surface_flow_type) :: this

end subroutine PMSurfaceInitializeRun


end module PM_Surface_Flow_class