module PM_SWE_class

#include "petsc/finclude/petscts.h"
  use petscts
  use PM_Base_class
  use PM_Surface_Flow_class

  implicit none

  private

  type, public, extends(pm_surface_flow_type) :: pm_swe_type
  contains
  end type pm_swe_type

  public :: PMSWECreate

contains

! ************************************************************************** !

function PMSWECreate()
  !
  ! Creates Shallow Water Equation process model shell
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !

  implicit none

  class(pm_swe_type), pointer :: PMSWECreate

  class(pm_swe_type), pointer :: this

  allocate(this)

  call PMSurfaceFlowInit(this)

  this%name = 'Shallow Water Equation'
  this%header = 'SWE'

  PMSWECreate => this

end function PMSWECreate


end module PM_SWE_class