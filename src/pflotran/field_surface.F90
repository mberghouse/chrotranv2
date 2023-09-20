module Field_Surface_module

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: field_surface_type

  Vec :: mannings0, mannings_loc

  Vec :: work, work_loc

  Vec :: area

  ! residual vectors
  Vec :: flow_r

  ! Solution vectors (yy = previous solution, xx = current iterate)
  Vec :: flow_xx, flow_xx_loc, flow_dxx, flow_yy, flow_accum

  end type field_surface_type

  public :: FieldSurfaceCreate, &
            FieldSurfaceDestroy

contains

! ************************************************************************** !

function FieldSurfaceCreate()
  !
  ! Allocates and initializes a new field_surface_type object
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !

  implicit none

  type(field_surface_type), pointer :: FieldSurfaceCreate

  type(field_surface_type), pointer :: field_surface

  allocate(field_surface)

  ! nullify PetscVecs
  field_surface%mannings0 = PETSC_NULL_VEC
  field_surface%mannings_loc = PETSC_NULL_VEC

  field_surface%work = PETSC_NULL_VEC
  field_surface%work_loc = PETSC_NULL_VEC

  field_surface%area = PETSC_NULL_VEC

  field_surface%flow_r = PETSC_NULL_VEC
  field_surface%flow_xx = PETSC_NULL_VEC
  field_surface%flow_xx_loc = PETSC_NULL_VEC
  field_surface%flow_dxx = PETSC_NULL_VEC
  field_surface%flow_yy = PETSC_NULL_VEC
  field_surface%flow_accum = PETSC_NULL_VEC


  FieldSurfaceCreate => field_surface

end function FieldSurfaceCreate

! ************************************************************************** !

subroutine FieldSurfaceDestroy(field_surface)
  !
  ! Deallocates a field object
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !

  implicit none

  type(field_surface_type), pointer :: field_surface

  PetscErrorCode :: ierr

  ! Destroy PetscVecs
  if (field_surface%mannings0 /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%mannings0,ierr);CHKERRQ(ierr)
  endif
  if (field_surface%mannings_loc /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%mannings_loc,ierr);CHKERRQ(ierr)
  endif

  if (field_surface%work /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%work,ierr);CHKERRQ(ierr)
  endif
  if (field_surface%work_loc  /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%work_loc,ierr);CHKERRQ(ierr)
  endif

  if (field_surface%area  /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%area,ierr);CHKERRQ(ierr)
  endif

  if (field_surface%flow_r /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%flow_r,ierr);CHKERRQ(ierr)
  endif
  if (field_surface%flow_xx /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%flow_xx,ierr);CHKERRQ(ierr)
  endif
  if (field_surface%flow_xx_loc /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%flow_xx_loc,ierr);CHKERRQ(ierr)
  endif
  if (field_surface%flow_dxx /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%flow_dxx,ierr);CHKERRQ(ierr)
  endif
  if (field_surface%flow_yy /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%flow_yy,ierr);CHKERRQ(ierr)
  endif
  if (field_surface%flow_accum /= PETSC_NULL_VEC) then
    call VecDestroy(field_surface%flow_accum,ierr);CHKERRQ(ierr)
  endif

  if (associated(field_surface)) deallocate(field_surface)
  nullify(field_surface)

end subroutine FieldSurfaceDestroy

end module Field_Surface_module