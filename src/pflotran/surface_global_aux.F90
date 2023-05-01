module Surface_Global_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: surface_global_auxvar_type
    PetscInt :: istate
    PetscReal :: h      ! [m]
    PetscReal :: temp   ! [C]
    PetscReal :: den_kg ! [kg/m^3]
    PetscBool :: is_dry
  end type surface_global_auxvar_type
  
  type, public :: surface_global_type
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(surface_global_auxvar_type), pointer :: auxvars(:)
    type(surface_global_auxvar_type), pointer :: auxvars_bc(:)
    type(surface_global_auxvar_type), pointer :: auxvars_ss(:)
  end type surface_global_type
  
  interface SurfaceGlobalAuxVarDestroy
    module procedure SurfaceGlobalAuxVarSingleDestroy
    module procedure SurfaceGlobalAuxVarArrayDestroy
  end interface SurfaceGlobalAuxVarDestroy

  public :: SurfaceGlobalAuxCreate, &
            SurfaceGlobalAuxDestroy, &
            SurfaceGlobalAuxVarInit, &
            SurfaceGlobalAuxVarCopy, &
            SurfaceGlobalAuxVarDestroy, &
            SurfaceGlobalAuxVarStrip

contains

! ************************************************************************** !

function SurfaceGlobalAuxCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  ! 

  use Option_module

  implicit none
  
  type(surface_global_type), pointer :: SurfaceGlobalAuxCreate
  
  type(surface_global_type), pointer :: aux

  allocate(aux) 
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)

  SurfaceGlobalAuxCreate => aux
  
end function SurfaceGlobalAuxCreate

! ************************************************************************** !

subroutine SurfaceGlobalAuxVarInit(auxvar,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  ! 
  use Option_module
  use PFLOTRAN_Constants_module, only : DUMMY_VALUE

  implicit none
  
  type(surface_global_auxvar_type) :: auxvar
  type(option_type) :: option
  
  auxvar%istate = 0
  auxvar%is_dry = PETSC_FALSE
  auxvar%h = 0.d0
  auxvar%temp = option%flow%reference_temperature
  auxvar%den_kg = 0.d0

end subroutine SurfaceGlobalAuxVarInit

! ************************************************************************** !

subroutine SurfaceGlobalAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  ! 

  use Option_module

  implicit none
  
  type(surface_global_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%istate = auxvar%istate
  auxvar2%is_dry = auxvar%is_dry
  auxvar2%h = auxvar%h
  auxvar2%temp = auxvar%temp
  auxvar2%den_kg = auxvar%den_kg

end subroutine SurfaceGlobalAuxVarCopy

! ************************************************************************** !

subroutine SurfaceGlobalAuxVarSingleDestroy(auxvar)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  ! 

  implicit none

  type(surface_global_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call SurfaceGlobalAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine SurfaceGlobalAuxVarSingleDestroy

! ************************************************************************** !

subroutine SurfaceGlobalAuxVarArrayDestroy(auxvars)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  ! 

  implicit none

  type(surface_global_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call SurfaceGlobalAuxVarStrip(auxvars(iaux))
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine SurfaceGlobalAuxVarArrayDestroy

! ************************************************************************** !

subroutine SurfaceGlobalAuxVarStrip(auxvar)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  ! 

  use Utility_module, only: DeallocateArray

  implicit none

  type(surface_global_auxvar_type) :: auxvar

end subroutine SurfaceGlobalAuxVarStrip

! ************************************************************************** !

subroutine SurfaceGlobalAuxDestroy(aux)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  ! 

  implicit none

  type(surface_global_type), pointer :: aux
  
  if (.not.associated(aux)) return
  
  call SurfaceGlobalAuxVarDestroy(aux%auxvars)
  call SurfaceGlobalAuxVarDestroy(aux%auxvars_bc)
  call SurfaceGlobalAuxVarDestroy(aux%auxvars_ss)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine SurfaceGlobalAuxDestroy

end module Surface_Global_Aux_module
