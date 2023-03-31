module PM_Surface_Flow_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PM_Base_class
  use Realization_Surface_class

  implicit none

  private

  type, public, extends (pm_base_type) :: pm_surface_flow_type
    class(realization_surface_type), pointer :: surface_realization
  contains
    !procedure, public :: Setup => PMSurfaceFlowStep
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

end subroutine PMSurfaceFlowInit

end module PM_Surface_Flow_class