module PMC_Surface_class

#include "petsc/finclude/petscts.h"
  use petscts
  use PMC_Base_class
  use Realization_Surface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends (pmc_base_type) :: pmc_surface_type
    class(realization_surface_type), pointer :: surface_realization
  contains
    procedure, public :: Init => PMCSurfaceInit
  end type pmc_surface_type

  public :: PMCSurfaceCreate

contains

! ************************************************************************** !

function PMCSurfaceCreate()
  ! 
  ! Allocates and initializes a new process model coupler object
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  ! 

  implicit none
  
  class(pmc_surface_type), pointer :: PMCSurfaceCreate
  
  class(pmc_surface_type), pointer :: pmc
  
  allocate(pmc)
  call pmc%Init()
  
  PMCSurfaceCreate => pmc  
  
end function PMCSurfaceCreate

! ************************************************************************** !

subroutine PMCSurfaceInit(this)
  ! 
  ! Initializes a new process model coupler object
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  ! 

  implicit none
  
  class(pmc_surface_type) :: this
  
  call PMCBaseInit(this)
  nullify(this%surface_realization)

end subroutine PMCSurfaceInit

end module PMC_Surface_class