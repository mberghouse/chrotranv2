module PM_DWave_class

#include "petsc/finclude/petscts.h"
  use petscts
  use PM_Base_class
  use PM_Surface_Flow_class

  implicit none

  private

  type, public, extends(pm_surface_flow_type) :: pm_dwave_type
  contains
    procedure, public :: InitializeTimestep => PMDWaveInitializeTimestep
    procedure, public :: PreSolve => PMDWavePreSolve
    procedure, public :: RHSFunction => PMDWaveRHSFunction
    procedure, public :: PostSolve => PMDWavePostSolve
    procedure, public :: FinalizeTimeStep => PMDWaveFinalizeTimeSetup
  end type pm_dwave_type

  public :: PMDWaveCreate

contains

! ************************************************************************** !

function PMDWaveCreate()
  !
  ! Creates diffusion wave process model shell
  !
  ! Author: Gautam Bisht
  ! Date: 05/19/23
  !

  implicit none

  class(pm_dwave_type), pointer :: PMDWaveCreate

  class(pm_dwave_type), pointer :: this

  allocate(this)

  call PMSurfaceFlowInit(this)

  this%name = 'Diffusion wave'
  this%header = 'DWave'

  PMDWaveCreate => this

end function PMDWaveCreate

! ************************************************************************** !

subroutine PMDWaveInitializeTimestep(this)
  !
  !
  ! Author: Gautam Bisht
  ! Date: 05/19/23
  !

  implicit none

  class(pm_dwave_type) :: this

  this%option%flow_dt = this%option%dt

end subroutine PMDWaveInitializeTimestep

! ************************************************************************** !

subroutine PMDWavePreSolve(this)
  !
  !
  ! Author: Gautam Bisht
  ! Date: 05/19/23
  !

  implicit none

  class(pm_dwave_type) :: this

end subroutine PMDWavePreSolve

! ************************************************************************** !

subroutine PMDWaveRHSFunction(this,ts,time,xx,ff,ierr)

  !
  !
  ! Author: Gautam Bisht
  ! Date: 05/19/23
  !
  use DWave_module

  implicit none

  class(pm_dwave_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  PetscErrorCode :: ierr

  call DWaveRHSFunction(ts,time,xx,ff,this%surface_realization,ierr);CHKERRQ(ierr)

end subroutine PMDWaveRHSFunction

! ************************************************************************** !

subroutine PMDWavePostSolve(this)
  !
  !
  ! Author: Gautam Bisht
  ! Date: 05/19/23
  !

  implicit none

  class(pm_dwave_type) :: this

end subroutine PMDWavePostSolve

! ************************************************************************** !

subroutine PMDWaveFinalizeTimeSetup(this)
  !
  !
  ! Author: Gautam Bisht
  ! Date: 05/19/23
  !

  implicit none

  class(pm_dwave_type) :: this

end subroutine PMDWaveFinalizeTimeSetup

end module PM_DWave_class