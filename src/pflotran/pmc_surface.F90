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
    procedure, public :: SetupSolvers => PMCSurfaceSetupSolvers
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

! ************************************************************************** !

subroutine PMCSurfaceSetupSolvers(this)
  !
  ! Sets up the PETSc solver
  !
  ! Author: Gautam Bisht
  ! Date: 04/05/23
  !
  use Option_module
  use Timestepper_TS_class

  implicit none

  class(pmc_surface_type) :: this

  type(option_type), pointer :: option

  option => this%option

  if (associated(this%timestepper)) then
    select type(ts => this%timestepper)
      class is(timestepper_TS_type)
        call PMCSurfaceSetupSolvers_TS(this)
      class default
        option%io_buffer = &
          'Unknown timestepper found in PMCSubsurfaceSetupSolvers '
        call PrintErrMsg(option)
    end select
  endif

end subroutine PMCSurfaceSetupSolvers

! ************************************************************************** !

subroutine PMCSurfaceSetupSolvers_TS(this)
  !
  ! Sets up the PETSc TS solver
  !
  ! Author: Gautam Bisht
  ! Date: 04/05/23
  !
  use Convergence_module
  use Discretization_module
  use Option_module
  use PMC_Base_class
  use PM_Base_Pointer_module
  use PM_Base_class
  use PM_Surface_Flow_class
  use PM_SWE_class
  use Solver_module
  use Timestepper_Base_class
  use Timestepper_TS_class

  implicit none

  class(pmc_surface_type) :: this

  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  PetscErrorCode :: ierr

  option => this%option

  select type(ts => this%timestepper)
    class is(timestepper_TS_type)
      solver => ts%solver
    class default
      call PrintErrMsg(option,"Attempting to set up PETSc TS when" // &
       " timestepper is not of TS type")
  end select

  call SolverCreateTS(solver,option%mycomm)

  call TSSetProblemType(solver%ts,TS_NONLINEAR,ierr);CHKERRQ(ierr)

  call TSSetType(solver%ts,TSEULER,ierr);CHKERRQ(ierr)

  select type(pm => this%pm_ptr%pm)

    class is (pm_swe_type)
      call PrintMsg(option,"  Beginning setup of SWE TS ")

      select case(option%iflowmode)
        case(SWE_MODE)
        case default
          option%io_buffer = 'Timestepper TS unsupported for mode: '// option%flowmode
          call PrintErrMsg(option)
      end select

      write(option%io_buffer,'(" number of dofs = ",i3,", number of &
                &phases = ",i3,i2)') option%nflowdof,option%nphase
      call PrintMsg(option)
      select case(option%iflowmode)
        case(SWE_MODE)
          option%io_buffer = " mode = SWE: h, hu, hv"
      end select
      call PrintMsg(option)

      call TSSetOptionsPrefix(solver%ts,"swe_",ierr);CHKERRQ(ierr)
      call TSSetFromOptions(solver%ts,ierr);CHKERRQ(ierr)

      call SolverCheckCommandLine(solver)

      call TSSetSolution(solver%ts,this%pm_ptr%pm%solution_vec,ierr);CHKERRQ(ierr)
      call TSSetRHSFunction(solver%ts,this%pm_ptr%pm%residual_vec, &
                            PMRHSFunctionPtr,this%pm_ptr,ierr);CHKERRQ(ierr)

      call PrintMsg(option,"  Finished setting up FLOW SNES ")

    class default

      option%io_buffer = 'Timestepper TS supported only for pm_swe_type '
      call PrintErrMsg(option)

  end select

end subroutine PMCSurfaceSetupSolvers_TS

end module PMC_Surface_class