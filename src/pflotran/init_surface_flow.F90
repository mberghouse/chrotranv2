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
  use Option_module
  use Patch_module
  use Realization_Surface_class
  use Simulation_Surface_class

  implicit none

  class(simulation_surface_type) :: simulation

  class(realization_surface_type), pointer :: surface_realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch

  surface_realization => simulation%surface_realization
  option => surface_realization%option
  patch => surface_realization%patch


  write(*,*)'stopping in InitSurfaceFlowSetupRealization'
  call exit(0)

end subroutine InitSurfaceFlowSetupRealization

end module Init_Surface_Flow_module
