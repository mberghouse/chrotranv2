module SWE_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  type, public :: swe_auxvar_type
    PetscReal :: hu ! momentum in x-dir
    PetscReal :: hv ! momentum in y-dir
    PetscReal :: u  ! velocity in x-dir
    PetscReal :: v  ! velocity in y-dir
  end type swe_auxvar_type

  type, public :: swe_type
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(swe_auxvar_type), pointer :: auxvars(:)
    type(swe_auxvar_type), pointer :: auxvars_bc(:)
    type(swe_auxvar_type), pointer :: auxvars_ss(:)
  end type swe_type

  public :: SWEAuxCreate, &
            SWEAuxVarInit, &
            SWEAuxVarCompute

contains

! ************************************************************************** !
function SWEAuxCreate(option)
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(swe_type), pointer :: SWEAuxCreate

  type(swe_type), pointer :: aux

  allocate(aux)

  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0

  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)

  SWEAuxCreate => aux

end function SWEAuxCreate

! ************************************************************************** !
subroutine SWEAuxVarStrip(auxvar)
  !
  ! Delllocates an swe_auxvar_type object
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !

  implicit none

  type(swe_auxvar_type) :: auxvar

end subroutine SWEAuxVarStrip

! ************************************************************************** !
subroutine SWEAuxVarArrayDestroy(auxvars)
  !
  ! Delllocates an array of swe_auxvar_type objects
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !

  implicit none

  type(swe_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
  do iaux = 1, size(auxvars)
    call SWEAuxVarStrip(auxvars(iaux))
  enddo
  deallocate(auxvars)
  endif

end subroutine SWEAuxVarArrayDestroy

! ************************************************************************** !
subroutine SWEAuxDestroy(aux)
  !
  ! Delllocates an swe_type object
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !

  implicit none

  type(swe_type), pointer :: aux

  if (.not.associated(aux)) return

  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0

  call SWEAuxVarArrayDestroy(aux%auxvars)
  call SWEAuxVarArrayDestroy(aux%auxvars_bc)
  call SWEAuxVarArrayDestroy(aux%auxvars_ss)

  deallocate(aux)
  nullify(aux)

end subroutine SWEAuxDestroy

! ************************************************************************** !
subroutine SWEAuxVarCompute(x,auxvar,surf_global_auxvar,option)
  !
  ! Delllocates an swe_type object
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use Option_module
  use Surface_Global_Aux_module

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  type(swe_auxvar_type) :: auxvar
  type(surface_global_auxvar_type) :: surf_global_auxvar

end subroutine SWEAuxVarCompute

! ************************************************************************** !
subroutine SWEAuxVarInit(auxvar)
  !
  ! Initializes an swe_auxvar_type object
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !

  implicit none

  type(swe_auxvar_type) :: auxvar

  auxvar%hu = 0.d0
  auxvar%hv = 0.d0
  auxvar%u = 0.d0
  auxvar%v = 0.d0

end subroutine SWEAuxVarInit

end module SWE_Aux_module