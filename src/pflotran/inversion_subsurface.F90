module Inversion_Subsurface_class

#include "petsc/finclude/petscmat.h"
  use petscmat

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Realization_Subsurface_class
  use Simulation_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_base_type) :: inversion_subsurface_type
    character(len=MAXSTRINGLENGTH) :: forward_simulation_filename
    class(simulation_subsurface_type), pointer :: forward_simulation
    class(realization_subsurface_type), pointer :: realization
    Mat :: Jsensitivity
    Mat :: JsensitivityT
    PetscReal, pointer :: measurement(:)
    PetscInt, pointer :: imeasurement(:)
    Vec :: quantity_of_interest
    PetscInt :: iqoi
    Vec :: ref_quantity_of_interest
    character(len=MAXWORDLENGTH) :: ref_qoi_dataset_name
  contains
    procedure, public :: Init => InversionSubsurfaceInit
  end type inversion_subsurface_type

  public :: InversionSubsurfaceInit, &
            InversionSubsurfReadSelectCase, &
            InversionSubsurfaceSetup, &
            InversionSubsurfaceStrip

contains

! ************************************************************************** !

subroutine InversionSubsurfaceInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY
  use Driver_module

  class(inversion_subsurface_type) :: this
  class(driver_type), pointer :: driver

  call InversionBaseInit(this,driver)

  this%quantity_of_interest = PETSC_NULL_VEC
  this%iqoi = UNINITIALIZED_INTEGER
  this%ref_quantity_of_interest = PETSC_NULL_VEC
  this%ref_qoi_dataset_name = ''
  this%forward_simulation_filename = ''

  this%Jsensitivity = PETSC_NULL_MAT
  this%JsensitivityT = PETSC_NULL_MAT
  nullify(this%measurement)
  nullify(this%imeasurement)

  nullify(this%forward_simulation)
  nullify(this%realization)

end subroutine InversionSubsurfaceInit

! ************************************************************************** !

subroutine InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                          error_string,option)

  use Input_Aux_module
  use Option_module
  use String_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, POROSITY
  use Utility_module

  class(inversion_subsurface_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, pointer :: tempint(:)
  PetscReal, pointer :: tempreal(:)

  found = PETSC_TRUE
  call InversionBaseReadSelectCase(this,input,keyword,found, &
                                   error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('FORWARD_SIMULATION_FILENAME')
      call InputReadFilename(input,option,this%forward_simulation_filename)
      call InputErrorMsg(input,option,keyword,error_string)
    case('QUANTITY_OF_INTEREST')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,keyword,error_string)
      call StringToUpper(word)
      select case(word)
        case('ELECTRICAL_CONDUCTIVITY')
          this%iqoi = ELECTRICAL_CONDUCTIVITY
        case('PERMEABILITY')
          this%iqoi = PERMEABILITY
        case('POROSITY')
          this%iqoi = POROSITY
        case default
          call InputKeywordUnrecognized(input,word,trim(error_string)// &
                                        & ',QUANTITY_OF_INTEREST',option)
      end select
    case('REFERENCE_QUANTITY_OF_INTEREST')
      call InputReadNChars(input,option,this%ref_qoi_dataset_name, &
                           MAXWORDLENGTH,PETSC_TRUE)
      call InputErrorMsg(input,option,'DATASET NAME', &
                         keyword)
    case('MEASUREMENTS')
      string = trim(error_string)//keyword
      i = 10
      allocate(tempint(i))
      tempint = UNINITIALIZED_INTEGER
      allocate(tempreal(i))
      tempreal = UNINITIALIZED_DOUBLE
      i = 0
      do
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,error_string)
        if (InputCheckExit(input,option)) exit
        i = i + 1
        if (i > size(tempint)) then
          call ReallocateArray(tempint)
          call ReallocateArray(tempreal)
        endif
        call InputReadInt(input,option,tempint(i))
        call InputErrorMsg(input,option,'cell id',string)
        call InputReadDouble(input,option,tempreal(i))
        call InputErrorMsg(input,option,'measurement',string)
      enddo
      allocate(this%imeasurement(i))
      this%imeasurement(:) = tempint(1:i)
      allocate(this%measurement(i))
      this%measurement(:) = tempreal(1:i)
      call DeallocateArray(tempint)
      call DeallocateArray(tempreal)
    case default
      found = PETSC_FALSE
  end select

end subroutine InversionSubsurfReadSelectCase

! ************************************************************************** !

subroutine InversionSubsurfaceSetup(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Discretization_module
  use Material_module

  class(inversion_subsurface_type) :: this

  PetscInt :: num_measurements
  PetscErrorCode :: ierr

    ! on first pass, store and set thereafter
  if (this%quantity_of_interest == PETSC_NULL_VEC) then
    num_measurements = size(this%imeasurement)
    call VecDuplicate(this%realization%field%work, &
                      this%quantity_of_interest,ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(this%realization%discretization, &
                                     this%realization%field%work_loc, &
                                     this%quantity_of_interest,ONEDOF)
    call MatCreateDense(this%driver%comm%mycomm, &
                        num_measurements, &
                        this%realization%patch%grid%nlmax, &
                        num_measurements, &
                        this%realization%patch%grid%nmax, &
                        PETSC_NULL_SCALAR, &
                        this%Jsensitivity,ierr);CHKERRQ(ierr)
    call MatCreateDense(this%driver%comm%mycomm, &
                        this%realization%patch%grid%nlmax, &
                        PETSC_DECIDE, &
                        this%realization%patch%grid%nmax, &
                        num_measurements, &
                        PETSC_NULL_SCALAR, &
                        this%JsensitivityT,ierr);CHKERRQ(ierr)
  endif

end subroutine InversionSubsurfaceSetup

! ************************************************************************** !

subroutine InversionSubsurfaceStrip(this)
  !
  ! Deallocates members of inversion Subsurface
  !
  ! Author: Glenn hammond
  ! Date: 09/20/21
  !
  use Utility_module

  class(inversion_subsurface_type) :: this

  PetscErrorCode :: ierr

  call InversionBaseStrip(this)

  nullify(this%realization)
  if (associated(this%forward_simulation)) then
    print *, 'Why is forward simulation still associated in &
             &InversionSubSurfStrip?'
    stop
  endif
  nullify(this%forward_simulation)
  if (this%quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%quantity_of_interest,ierr);CHKERRQ(ierr)
  endif
  if (this%ref_quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%ref_quantity_of_interest,ierr);CHKERRQ(ierr)
  endif

  call DeallocateArray(this%imeasurement)
  call DeallocateArray(this%measurement)

end subroutine InversionSubsurfaceStrip

end module Inversion_Subsurface_class
