module Factory_Surface_Read_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Simulation_Surface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: &
            FactorySurfaceReadFlowPM, &
            FactorySurfaceReadRequiredCards
contains

! ************************************************************************** !

subroutine FactorySurfaceReadFlowPM(input,option,pm)
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_SWE_class
  use Init_Common_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,SURFACE_FLOW'
  write(*,*)'++++ FactorySurfaceReadFlowPM'

  nullify(pm)
  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    write(*,*)'word: ',trim(word)
    select case(word)
      case('MODE')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'mode',error_string)
        call StringToUpper(word)

        select case(word)
          case('SWE')
            write(*,*)'word: ',trim(word)
            pm => PMSWECreate()
          case default
            error_string = trim(error_string) // ',MODE'
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
        pm%option => option
      case('OPTIONS')
        if (.not.associated(pm)) then
          option%io_buffer = 'MODE keyword must be read first under ' // &
                             trim(error_string)
          call PrintErrMsg(option)
        endif
        call pm%ReadSimulationOptionsBlock(input)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (.not.associated(pm)) then
    option%io_buffer = 'A flow MODE (card) must be included in the &
      &SURFACE_FLOW block in ' // trim(error_string) // '.'
    call PrintErrMsg(option)
  endif

end subroutine FactorySurfaceReadFlowPM

! ************************************************************************** !

subroutine FactorySurfaceReadRequiredCards(simulation,input)
  !
  ! Reads required cards from input file
  !
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !
  use Option_module
  use Discretization_module
  use Grid_module
  use Input_Aux_module
  use String_module
  use Patch_module
  use Realization_Surface_class
  use HDF5_Aux_module

  use Simulation_Surface_class
  use General_module
  use Reaction_module
  use Reaction_Aux_module
  use NW_Transport_Aux_module
  use Init_Common_module

  implicit none

  class(simulation_surface_type) :: simulation

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  class(realization_surface_type), pointer :: surface_realization
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  PetscBool :: found
  PetscBool :: qerr

  character(len = MAXSTRINGLENGTH) :: wname

  surface_realization => simulation%surface_realization
  patch => surface_realization%patch
  option => surface_realization%option
  discretization => surface_realization%discretization

  qerr  = PETSC_FALSE
  wname = '<missing>'
  found = PETSC_FALSE

  write(*,*)'>>>>>>>>>> In FactorySurfaceReadRequiredCards'
  call InputPushBlock(input,'SURFACE',option)

  ! GRID information - GRID is a required card for every simulation
  string = "GRID"
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call InputPushBlock(input,'GRID',option)
  call DiscretizationReadRequiredCards(discretization,input,option)

  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      patch => PatchCreate()
      patch%grid => discretization%grid
      if (.not.associated(surface_realization%patch_list)) then
        surface_realization%patch_list => PatchCreateList()
      endif
      call PatchAddToList(patch,surface_realization%patch_list)
      surface_realization%patch => patch
  end select
  call InputPopBlock(input,option)

  call InputPopBlock(input,option)

end subroutine FactorySurfaceReadRequiredCards

end module Factory_Surface_Read_module
