module SWE_module

#include "petsc/finclude/petscts.h"
  use petscts

  implicit none

  private

  public :: SWESetup

contains

! ************************************************************************** !

subroutine SWESetup(surface_realization)
  !
  ! Sets up the variable list for output and observation.
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use Realization_Surface_class
  use Output_Aux_module

  implicit none

  class (realization_surface_type) :: surface_realization
  type(output_variable_list_type), pointer :: list

  list => surface_realization%output_option%output_snap_variable_list
  call SWESetPlotVariables(list)

  list => surface_realization%output_option%output_obs_variable_list
  call SWESetPlotVariables(list)

end subroutine SWESetup

! ************************************************************************** !

subroutine SWESetPlotVariables(list)
  !
  ! Adds variables to be printed to list
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use Output_Aux_module
  use Variables_module
  use PFLOTRAN_Constants_module

  implicit none

  type(output_variable_list_type), pointer :: list

  character(len=MAXWORDLENGTH) :: name, units

  if (associated(list%first)) then
    return
  endif

  name = 'H'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_LIQUID_HEAD)

  name = 'hu'
  units = 'm^2/s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_FLOW_X_MOMENTUM)

  name = 'hv'
  units = 'm^2/s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_FLOW_Y_MOMENTUM)

end subroutine SWESetPlotVariables

end module SWE_module