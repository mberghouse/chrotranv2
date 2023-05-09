module PM_SWE_class

#include "petsc/finclude/petscts.h"
  use petscts
  use PM_Base_class
  use PM_Surface_Flow_class

  implicit none

  private

  type, public, extends(pm_surface_flow_type) :: pm_swe_type
  contains
    procedure, public :: InitializeTimestep => PMSWEInitializeTimestep
    procedure, public :: PreSolve => PMSWEPreSolve
    procedure, public :: RHSFunction => PMSWERHSFunction
    procedure, public :: PostSolve => PMSWEPostSolve
    procedure, public :: FinalizeTimeStep => PMSWEFinalizeTimeSetup
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

! ************************************************************************** !

subroutine PMSWEInitializeTimestep(this)
  !
  !
  ! Author: Gautam Bisht
  ! Date: 04/05/23
  !

  implicit none

  class(pm_swe_type) :: this

end subroutine PMSWEInitializeTimestep

! ************************************************************************** !

subroutine PMSWEPreSolve(this)
  !
  !
  ! Author: Gautam Bisht
  ! Date: 04/05/23
  !

  implicit none

  class(pm_swe_type) :: this

end subroutine PMSWEPreSolve

! ************************************************************************** !

subroutine PMSWERHSFunction(this,ts,time,xx,ff,ierr)

  !
  !
  ! Author: Gautam Bisht
  ! Date: 04/05/23
  !
  use SWE_module

  implicit none

  class(pm_swe_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  PetscErrorCode :: ierr

  call SWERHSFunction(ts,time,xx,ff,this%surface_realization,ierr);CHKERRQ(ierr)

end subroutine PMSWERHSFunction

! ************************************************************************** !

subroutine PMSWEPostSolve(this)
  !
  !
  ! Author: Gautam Bisht
  ! Date: 05/08/23
  !

  implicit none

  class(pm_swe_type) :: this

end subroutine PMSWEPostSolve

! ************************************************************************** !

subroutine PMSWEFinalizeTimeSetup(this)
  !
  !
  ! Author: Gautam Bisht
  ! Date: 05/08/23
  !

  implicit none

  class(pm_swe_type) :: this

end subroutine PMSWEFinalizeTimeSetup


end module PM_SWE_class