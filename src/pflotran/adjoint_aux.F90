module Adjoint_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: adjoint_auxvar_type
    ! Simulated fields
    PetscReal :: pres ! liquid pressure
    PetscReal :: sat  ! liquid saturation
    PetscReal :: pc   ! capillary pressure
    PetscReal :: kr   ! relative permeability
    PetscReal :: effective_porosity
    PetscReal :: dpor_dp
    PetscReal :: dsat_dp ! derivative of saturation wrt pressure
    PetscReal :: dkr_dp  ! derivative of rel. perm. wrt pressure
    ! For adjoint fields
    PetscReal :: plambda           ! Adjoint pressure field
    PetscReal, pointer :: dAdk(:)  ! System matrix derivative dA/dKin
    PetscReal, pointer :: dBdk(:)  ! Accumulation matrix derivative dB/dKin
    PetscReal, pointer :: dcdk(:)  ! source vector derivative dc/dkin
  end type adjoint_auxvar_type

  type, public :: adjoint_type
    PetscInt :: num_aux, num_aux_bc
    type(adjoint_auxvar_type), pointer :: auxvars(:)
    type(adjoint_auxvar_type), pointer :: auxvars_bc(:)
  end type adjoint_type

contains

! ************************************************************************** !

function AdjointAuxCreate()
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 09/23/21
  !

  implicit none

  type(adjoint_type), pointer :: AdjointAuxCreate

  type(adjoint_type), pointer :: aux

  allocate(aux)
  aux%num_aux = 0
  aux%num_aux_bc = 0

  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)

  AdjointAuxCreate => aux

end function AdjointAuxCreate

! ************************************************************************** !

subroutine AdjointAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 09/23/21
  !

  use Option_module

  implicit none

  type(adjoint_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%pres = 0.d0
  auxvar%sat = 0.d0
  auxvar%pc = 0.d0
  auxvar%kr = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%dpor_dp = 0.d0
  auxvar%dsat_dp = 0.d0
  auxvar%dkr_dp = 0.d0

  auxvar%plambda = 0.d0

  nullify(auxvar%dAdk)
  nullify(auxvar%dBdk)
  nullify(auxvar%dcdk)

end subroutine AdjointAuxVarInit

! ************************************************************************** !

subroutine AdjointAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 09/23/21
  !

  use Option_module

  implicit none

  type(adjoint_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pres = auxvar%pres
  auxvar2%sat = auxvar%sat
  auxvar2%pc = auxvar%pc
  auxvar2%kr = auxvar%kr
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%dpor_dp = auxvar%dpor_dp
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dkr_dp = auxvar%dkr_dp

  auxvar2%plambda = auxvar%plambda

  ! pj - TODO: Check derivative matrices and vector if associated then copy

end subroutine AdjointAuxVarCopy

! ************************************************************************** !

subroutine AdjointAuxVarDestroy(auxvars)
  !
  ! Deallocates an auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 09/23/21
  !
  use Utility_module, only: DeallocateArray

  implicit none

  type(adjoint_auxvar_type), pointer :: auxvars(:)

  PetscInt :: iaux

  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call DeallocateArray(auxvars(iaux)%dAdk)
      call DeallocateArray(auxvars(iaux)%dBdk)
      call DeallocateArray(auxvars(iaux)%dcdk)
    enddo
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine AdjointAuxVarDestroy

! ************************************************************************** !

subroutine AdjointAuxDestroy(aux)
  !
  ! Deallocates an adjoint auxiliary object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 09/23/21
  !

  implicit none

  type(adjoint_type), pointer :: aux

  if (.not.associated(aux)) return

  call AdjointAuxVarDestroy(aux%auxvars)
  call AdjointAuxVarDestroy(aux%auxvars_bc)

  deallocate(aux)
  nullify(aux)

end subroutine AdjointAuxDestroy

end module Adjoint_Aux_module