module Surface_Auxiliary_module

  use Surface_Global_Aux_module
!  use Surface_Flow_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: surface_auxiliary_type
    type(surface_global_type), pointer :: SurfaceGlobal
  end type surface_auxiliary_type
  
  public :: SurfaceAuxInit, &
            SurfaceAuxDestroy

contains

! ************************************************************************** !

subroutine SurfaceAuxInit(surf_aux)
  ! 
  ! This routine initializes a surface-auxiliary object
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  ! 

  implicit none
  
  type(surface_auxiliary_type) :: surf_aux
  
  nullify(surf_aux%SurfaceGlobal)
  
end subroutine SurfaceAuxInit

! ************************************************************************** !

subroutine SurfaceAuxDestroy(surf_aux)
  ! 
  ! This routine deallocates pointers in a surface-auxiliary object
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  ! 

  implicit none
  
  type(surface_auxiliary_type) :: surf_aux
  
  call SurfaceGlobalAuxDestroy(surf_aux%SurfaceGlobal)

  nullify(surf_aux%SurfaceGlobal)

end subroutine SurfaceAuxDestroy

end module Surface_Auxiliary_module
