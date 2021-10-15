module Inversion_Aux_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: inversion_auxvar_type
  contains
  end type inversion_auxvar_type

  type, public :: inversion_aux_type
    PetscInt, pointer :: imeasurement(:)
    PetscReal, pointer :: measurement(:)
    Mat :: Jsensitivity
    PetscInt :: num_aux
    class(inversion_auxvar_type), pointer :: auxvars(:)
  end type inversion_aux_type

  public :: InversionAuxCreate, &
            InversionAuxVarInit, &
            InversionAuxVarStrip, &
            InversionAuxDestroy

contains

! ************************************************************************** !

function InversionAuxCreate()
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !

  use Option_module

  implicit none

  type(inversion_aux_type), pointer :: InversionAuxCreate

  type(inversion_aux_type), pointer :: aux

  allocate(aux)
  nullify(aux%auxvars)

  aux%num_aux = 0
  nullify(aux%imeasurement)
  nullify(aux%measurement)

  aux%Jsensitivity = PETSC_NULL_MAT

  InversionAuxCreate => aux

end function InversionAuxCreate

! ************************************************************************** !

subroutine InversionAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  use Option_module

  implicit none

  class(inversion_auxvar_type) :: auxvar
  type(option_type) :: option

end subroutine InversionAuxVarInit

! ************************************************************************** !

subroutine InversionAuxVarStrip(auxvar)
  !
  ! Deallocates a inversion auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  use Utility_module, only : DeallocateArray

  implicit none

  class(inversion_auxvar_type) :: auxvar

end subroutine InversionAuxVarStrip

! ************************************************************************** !

subroutine InversionAuxDestroy(aux)
  !
  ! Deallocates a inversion auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 09/22/21
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(inversion_aux_type), pointer :: aux

  PetscInt :: iaux

  if (.not.associated(aux)) return

  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call InversionAuxVarStrip(aux%auxvars(iaux))
    enddo
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)

  ! these objects are destroyed elsewhere, do not destroy
  nullify(aux%imeasurement)
  nullify(aux%measurement)
  aux%Jsensitivity = PETSC_NULL_MAT

  deallocate(aux)
  nullify(aux)

end subroutine InversionAuxDestroy

end module Inversion_Aux_module
