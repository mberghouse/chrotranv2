module DWave_module

#include "petsc/finclude/petscts.h"
  use petscts

  implicit none

  private

  public :: DWaveSetup, &
            DWaveUpdateAuxVars, &
            DWaveRHSFunction

contains

! ************************************************************************** !

subroutine DWaveSetup(surface_realization)
  !
  ! Sets up the variable list for output and observation.
  !
  ! Author: Gautam Bisht
  ! Date: 05/19/23
  !
  use Realization_Surface_class
  use Output_Aux_module
  use Patch_module
  use Option_module
  use Grid_module
  use Coupler_module

  implicit none

  class (realization_surface_type) :: surface_realization
  type(output_variable_list_type), pointer :: list

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid

  list => surface_realization%output_option%output_snap_variable_list
  call DWaveSetPlotVariables(list)

  list => surface_realization%output_option%output_obs_variable_list
  call DWaveSetPlotVariables(list)

  option => surface_realization%option
  patch => surface_realization%patch
  grid => patch%grid

end subroutine DWaveSetup

! ************************************************************************** !

subroutine DWaveSetPlotVariables(list)
  !
  ! Adds variables to be printed to list
  !
  ! Author: Gautam Bisht
  ! Date: 05/19/23
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

end subroutine DWaveSetPlotVariables

! ************************************************************************** !

subroutine DWaveUpdateAuxVars(surface_realization)
  !
  ! Sets up the variable list for output and observation.
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use Realization_Surface_class
  use Patch_module
  use Option_module
  use Field_Surface_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Surface_Global_Aux_module

  implicit none

  class (realization_surface_type) :: surface_realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_surface_type), pointer :: field_surface
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_ss(:)

  PetscInt :: ghosted_id, istart, iend
  PetscReal, pointer :: xx_loc_p(:)
  PetscErrorCode :: ierr

  option => surface_realization%option
  patch => surface_realization%patch
  grid => patch%grid
  field_surface => surface_realization%field_surface

  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc
  surf_global_auxvars_ss => patch%surf_aux%SurfaceGlobal%auxvars_ss

  call VecGetArrayF90(field_surface%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif

    surf_global_auxvars(ghosted_id)%h = xx_loc_p(ghosted_id)

  enddo

  call VecRestoreArrayF90(field_surface%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

end subroutine DWaveUpdateAuxVars

! ************************************************************************** !

subroutine DWaveRHSFunction(ts,time,x,f,surface_realization,ierr)
  !
  ! Sets up the variable list for output and observation.
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use PFLOTRAN_Constants_module
  use Realization_Surface_class
  use Option_module
  use Field_Surface_module
  use Discretization_module

  implicit none

  TS :: ts
  PetscReal :: time
  Vec :: x
  Vec :: f
  class (realization_surface_type) :: surface_realization
  PetscReal :: max_courant_num
  PetscErrorCode :: ierr

  type(discretization_type), pointer :: discretization
  type(field_surface_type), pointer :: field_surface
  type(option_type), pointer :: option

  field_surface => surface_realization%field_surface
  discretization => surface_realization%discretization
  option => surface_realization%option

  call VecZeroEntries(f,ierr);CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(discretization,x,field_surface%flow_xx_loc,NFLOWDOF)

  call DWaveUpdateAuxVars(surface_realization)

  write(*,*) 'Stopping in DWaveRHSFunction. Add relevant physics.'
  call exit(0)

end subroutine DWaveRHSFunction

end module DWave_module