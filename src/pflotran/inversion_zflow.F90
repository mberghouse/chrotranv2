module Inversion_ZFlow_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Inversion_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_subsurface_type) :: inversion_zflow_type


  contains
    procedure, public :: Init => InversionZFlowInit
    procedure, public :: ReadBlock => InversionZFlowReadBlock
    procedure, public :: Initialize => InversionZFlowInitialize
    procedure, public :: Step => InversionZFlowStep
    procedure, public :: Finalize => InversionZFlowFinalize
    procedure, public :: Strip => InversionZFlowStrip
  end type inversion_zflow_type

  public :: InversionZFlowCreate, &
            InversionZFlowFinalize, &
            InversionZFlowStrip, &
            InversionZFlowDestroy

contains

! ************************************************************************** !

function InversionZFlowCreate(driver)
  !
  ! Allocates and initializes a new zflow inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !

  use Driver_module

  class(driver_type), pointer :: driver

  class(inversion_zflow_type), pointer :: InversionZFlowCreate

  allocate(InversionZFlowCreate)
  call InversionZFlowCreate%Init(driver)

end function InversionZFlowCreate

! ************************************************************************** !

subroutine InversionZFlowInit(this,driver)
  !
  ! Initializes a new zflow inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !

  use Driver_module

  class(inversion_zflow_type) :: this
  class(driver_type), pointer :: driver

  call InversionSubsurfaceInit(this,driver)

end subroutine InversionZFlowInit

! ************************************************************************** !

subroutine InversionZFlowReadBlock(this,input,option)
  !
  ! Reads input file parameters associated an ZFlow inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !

  use Input_Aux_module
  use Option_module
  use String_module

  class(inversion_zflow_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'ZFlow Inversion'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_TRUE
    call InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine InversionZFlowReadBlock

! ************************************************************************** !

subroutine InversionZFlowInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !
  class(inversion_zflow_type) :: this

  PetscInt :: num_measurements

  call InversionBaseInitialize(this)

end subroutine InversionZFlowInitialize

! ************************************************************************** !

subroutine InversionZFlowStep(this)
  !
  ! Execute a simulation
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21

  use Option_module
  use Factory_Forward_module

  class(inversion_zflow_type) :: this

  type(option_type), pointer :: option

  option => OptionCreate()
  write(option%group_prefix,'(i6)') this%iteration
  option%group_prefix = 'Run' // trim(adjustl(option%group_prefix))
  call OptionSetDriver(option,this%driver)
  call FactoryForwardInitialize(this%forward_simulation, &
                                this%forward_simulation_filename,option)
  this%realization => this%forward_simulation%realization
  call this%Initialize()
  call this%forward_simulation%InitializeRun()
  call InversionZFlowSetup(this)
  if (option%status == PROCEED) then
    call this%forward_simulation%ExecuteRun()
  endif
  !call InversionZFlowInvert(this)
  call this%forward_simulation%FinalizeRun()
  call this%forward_simulation%Strip()
  deallocate(this%forward_simulation)
  nullify(this%forward_simulation)

  this%converg_flag = PETSC_FALSE
  if (this%iteration > this%maximum_iteration) this%converg_flag = PETSC_TRUE

end subroutine InversionZFlowStep

! ************************************************************************** !

subroutine InversionZFlowSetup(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hamond
  ! Date: 10/11/21
  !
  use Discretization_module
  use Material_module
  use PM_ZFlow_class
  use Inversion_Aux_module

  class(inversion_zflow_type) :: this

  PetscErrorCode :: ierr

  call InversionSubsurfaceSetup(this)

  call DiscretizationGlobalToLocal(this%realization%discretization, &
                                   this%quantity_of_interest, &
                                   this%realization%field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                               this%realization%field%work_loc, &
                               this%iqoi,ZERO_INTEGER)

  ! this code must follow InversionSubsurfaceSetup() as Jsensitivity is
  ! is created within it.
  select type(pm => &
                this%forward_simulation%flow_process_model_coupler%pm_list)
    class is(pm_zflow_type)
      pm%inversion_aux => InversionAuxCreate()
      pm%inversion_aux%Jsensitivity = this%Jsensitivity
      pm%inversion_aux%imeasurement => this%imeasurement
    class default
      call this%driver%PrintErrMsg('Unsupported process model in &
                                   &InversionZFlowSetup.')
  end select

end subroutine InversionZFlowSetup

! ************************************************************************** !

subroutine InversionZFlowFinalize(this)
  !
  ! Finalizes inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !
  class(inversion_zflow_type) :: this

  call InversionBaseStrip(this)

end subroutine InversionZFlowFinalize

! ************************************************************************** !

subroutine InversionZFlowStrip(this)
  !
  ! Deallocates members of inversion zflow
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !
  class(inversion_zflow_type) :: this

   call InversionSubsurfaceStrip(this)

end subroutine InversionZFlowStrip

! ************************************************************************** !

subroutine InversionZFlowDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !
  class(inversion_zflow_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionZFlowDestroy

end module Inversion_ZFlow_class