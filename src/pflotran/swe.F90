module SWE_module

#include "petsc/finclude/petscts.h"
  use petscts

  implicit none

  private

  public :: SWESetup, &
            SWEUpdateAuxVars, &
            SWERHSFunction

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
  use Patch_module
  use Option_module
  use Grid_module
  use SWE_Aux_module

  implicit none

  class (realization_surface_type) :: surface_realization
  type(output_variable_list_type), pointer :: list

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(swe_auxvar_type), pointer :: swe_auxvars(:)
  PetscInt :: ghosted_id

  list => surface_realization%output_option%output_snap_variable_list
  call SWESetPlotVariables(list)

  list => surface_realization%output_option%output_obs_variable_list
  call SWESetPlotVariables(list)

  option => surface_realization%option
  patch => surface_realization%patch
  patch%surf_aux%SWE => SWEAuxCreate(option)
  grid => patch%grid

  allocate(swe_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call SWEAuxVarInit(swe_auxvars(ghosted_id))
  enddo
  patch%surf_aux%SWE%auxvars => swe_auxvars
  patch%surf_aux%SWE%num_aux = grid%ngmax

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

! ************************************************************************** !

subroutine SWEUpdateAuxVars(surface_realization)
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
  use SWE_Aux_module

  implicit none

  class (realization_surface_type) :: surface_realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_surface_type), pointer :: field_surface
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(swe_auxvar_type), pointer :: swe_auxvars(:)
  type(swe_auxvar_type), pointer :: swe_auxvars_bc(:)
  type(swe_auxvar_type), pointer :: swe_auxvars_ss(:)
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

  swe_auxvars => patch%surf_aux%SWE%auxvars
  swe_auxvars_bc => patch%surf_aux%SWE%auxvars_bc
  swe_auxvars_ss => patch%surf_aux%SWE%auxvars_ss

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

    iend = ghosted_id*option%nflowdof
    istart = iend - option%nflowdof + 1

    call SWEAuxVarCompute(xx_loc_p(istart:iend), swe_auxvars(ghosted_id), &
                          surf_global_auxvars(ghosted_id), option)

  enddo


  call VecRestoreArrayF90(field_surface%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

end subroutine SWEUpdateAuxVars

! ************************************************************************** !

subroutine SWERHSFunctionInternalConn(f,surface_realization,ierr)
  !
  ! Sets up the variable list for output and observation.
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use Realization_Surface_class
  use Option_module
  use Field_Surface_module
  use Field_Surface_module
  use Discretization_module
  use Patch_module
  use Grid_module
  use Connection_module

  implicit none

  Vec :: f
  class (realization_surface_type) :: surface_realization
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(discretization_type), pointer :: discretization
  type(field_surface_type), pointer :: field_surface
  type(option_type), pointer :: option
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn, sum_connection
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscReal, pointer :: f_p(:)

  field_surface => surface_realization%field_surface
  discretization => surface_realization%discretization
  option => surface_realization%option
  patch => surface_realization%patch
  grid => patch%grid

  write(*,*)'In SWERHSFunctionInternalConn'

  call VecGetArrayF90(f,f_p,ierr);CHKERRQ(ierr)

  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

    enddo
  enddo

  call VecRestoreArrayF90(f,f_p,ierr);CHKERRQ(ierr)

  write(*,*)'stopping in SWERHSFunctionInternalConn'
  call exit(0)

end subroutine SWERHSFunctionInternalConn

! ************************************************************************** !

subroutine SWERHSFunction(ts,time,x,f,surface_realization,ierr)
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
  use Field_Surface_module
  use Discretization_module

  implicit none

  TS :: ts
  PetscReal :: time
  Vec :: x
  Vec :: f
  class (realization_surface_type) :: surface_realization
  PetscErrorCode :: ierr

  type(discretization_type), pointer :: discretization
  type(field_surface_type), pointer :: field_surface
  type(option_type), pointer :: option

  field_surface => surface_realization%field_surface
  discretization => surface_realization%discretization
  option => surface_realization%option

  write(*,*)'In SWERHSFunction'
  call VecZeroEntries(f,ierr);CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(discretization,x,field_surface%flow_xx_loc,NFLOWDOF)

  call SWERHSFunctionInternalConn(f,surface_realization,ierr);CHKERRQ(ierr)

  write(*,*)'stopping in SWERHSFunction'
  call exit(0)

end subroutine SWERHSFunction

end module SWE_module
