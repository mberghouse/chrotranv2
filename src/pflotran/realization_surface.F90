module Realization_Surface_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Condition_module
  use Dataset_Base_class
  use Field_Surface_module
  use Option_module
  use PFLOTRAN_Constants_module
  use Realization_Base_class
  use Region_module
  use Material_Surface_module

  implicit none

  private

  type, public, extends(realization_base_type) :: realization_surface_type

    type(region_list_type), pointer :: surf_region_list
    type(condition_list_type),pointer :: surf_flow_conditions

    type(field_surface_type), pointer :: field_surface
    type(material_surface_property_type), pointer :: surf_material_properties

    character(len=MAXSTRINGLENGTH) :: surf_filename
    character(len=MAXSTRINGLENGTH) :: subsurf_filename

    class(dataset_base_type), pointer :: datasets
  
  end type realization_surface_type

  public :: RealizationSurfaceCreate, &
            RealizationSurfaceCreateDiscretization, &
            RealizationSurfacePassPtrsToPatches, &
            RealizationSurfaceProcessMatProp, &
            RealizationSurfaceProcessConditions, &
            RealizationSurfaceProcessCouplers

contains

! ************************************************************************** !

function RealizationSurfaceCreate(option)
  ! 
  ! This routine allocates and initializes a new SurfaceRealization object
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  ! 
  implicit none

  type(option_type), pointer :: option
  class(realization_surface_type),pointer :: RealizationSurfaceCreate
  class(realization_surface_type),pointer :: surf_realization
  
  allocate(surf_realization)
  call RealizationBaseInit(surf_realization,option)
  surf_realization%option => option

  allocate(surf_realization%surf_region_list)
  call RegionInitList(surf_realization%surf_region_list)
  
  allocate(surf_realization%surf_flow_conditions)
  call FlowConditionInitList(surf_realization%surf_flow_conditions)

  surf_realization%field_surface => FieldSurfaceCreate()

  nullify(surf_realization%surf_material_properties)
  nullify(surf_realization%datasets)

  RealizationSurfaceCreate => surf_realization

end function RealizationSurfaceCreate

! ************************************************************************** !

subroutine RealizationSurfaceCreateDiscretization(surf_realization)
  ! 
  ! This routine creates grid
  ! 
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  ! 
  use Grid_module
  use Grid_Unstructured_Aux_module, only : UGridMapIndices
  use Grid_Unstructured_module, only     : UGridEnsureRightHandRule
  use Coupler_module
  use Discretization_module
  use Grid_Unstructured_Cell_module
  use DM_Kludge_module
  
  implicit none

  class(realization_surface_type) :: surf_realization
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_surface_type), pointer :: field_surface
  type(option_type), pointer :: option
  PetscErrorCode :: ierr

  option => surf_realization%option
  field_surface => surf_realization%field_surface
  discretization => surf_realization%discretization

  write(*,*)'calling DiscretizationCreateDMs()'
  call DiscretizationCreateDMs(discretization, option%nflowdof, &
                               ZERO_INTEGER, ZERO_INTEGER, &
                               ZERO_INTEGER, ZERO_INTEGER, &
                               option)
  write(*,*)'calling DiscretizationCreateDMs() done'

  ! n degree of freedom, global
  call DiscretizationCreateVector(discretization,NFLOWDOF,field_surface%flow_xx,GLOBAL,option)
  call VecSet(field_surface%flow_xx,0.d0,ierr);CHKERRQ(ierr)
  write(*,*)'field_surface%flow_xx created'

  call DiscretizationDuplicateVector(discretization,field_surface%flow_xx,field_surface%flow_yy)
  call DiscretizationDuplicateVector(discretization,field_surface%flow_xx,field_surface%flow_dxx)
  call DiscretizationDuplicateVector(discretization,field_surface%flow_xx,field_surface%flow_r)
  call DiscretizationDuplicateVector(discretization,field_surface%flow_xx,field_surface%flow_accum)
  call DiscretizationDuplicateVector(discretization,field_surface%flow_xx,field_surface%work)

  ! 1 degree of freedom, global
  call DiscretizationCreateVector(discretization,ONEDOF,field_surface%mannings0,GLOBAL,option)
  call VecSet(field_surface%mannings0,0.d0,ierr);CHKERRQ(ierr)
  call DiscretizationDuplicateVector(discretization,field_surface%mannings0,field_surface%area)

  ! n degrees of freedom, local
  call DiscretizationCreateVector(discretization,NFLOWDOF,field_surface%flow_xx_loc,LOCAL,option)
  call VecSet(field_surface%flow_xx_loc,0.d0,ierr);CHKERRQ(ierr)
  call DiscretizationDuplicateVector(discretization,field_surface%flow_xx_loc,field_surface%work_loc)

  ! 1-dof degrees of freedom, local
  call DiscretizationCreateVector(discretization,ONEDOF,field_surface%mannings_loc,LOCAL,option)
  call VecSet(field_surface%mannings_loc,0.d0,ierr);CHKERRQ(ierr)

  grid => discretization%grid

  select case(discretization%itype)
    case(UNSTRUCTURED_GRID)
      ! set up nG2L, NL2G, etc.
      call GridMapIndices(grid, &
                          discretization%dm_1dof, &
                          discretization%stencil_type,&
                          option)
      call GridComputeCoordinates(grid,discretization%origin_global,option, &
                                  discretization%dm_1dof%ugdm) 
      call UGridEnsureRightHandRule(grid%unstructured_grid,grid%x, &
                                    grid%y,grid%z,grid%nG2A,grid%nL2G,option)
    case default
      option%io_buffer = 'Surface flow is only supported for an unstructured grid.'
      call PrintErrMsg(option)
  end select

end subroutine RealizationSurfaceCreateDiscretization

! ************************************************************************** !

subroutine RealizationSurfacePassPtrsToPatches(surf_realization)
  !
  ! Set patch%field => realization%field
  !
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  !

  implicit none

  class(realization_surface_type) :: surf_realization

  surf_realization%patch%field_surface => surf_realization%field_surface

end subroutine RealizationSurfacePassPtrsToPatches

! ************************************************************************** !

subroutine RealizationSurfaceProcessMatProp(surf_realization)
  !
  ! This routine sets up linkeage between surface material properties
  !
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  !

  use Option_module
  use Patch_module

  implicit none

  class(realization_surface_type) :: surf_realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch

  option => surf_realization%option
  patch => surf_realization%patch

  patch%surface_material_properties => surf_realization%surf_material_properties
  call MaterialSurfacePropConvertListToArray( &
                                patch%surface_material_properties, &
                                patch%surface_material_property_array, &
                                option)

  call MaterialSurfaceCreateIntToExtMapping(patch%surface_material_property_array, &
                                            patch%imat_internal_to_external)

end subroutine RealizationSurfaceProcessMatProp

! ************************************************************************** !

subroutine RealizationSurfaceProcessConditions(surf_realization)
  !
  ! Sets linkages of conditions. Presently, only flow condition is supported.
  !
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  !

  use Option_module
  use Patch_module

  implicit none

  class(realization_surface_type) :: surf_realization

  if (surf_realization%option%nflowdof > 0) then
    call RealizationSurfaceProcessFlowConditions(surf_realization)
  endif

end subroutine RealizationSurfaceProcessConditions

! ************************************************************************** !

subroutine RealizationSurfaceProcessFlowConditions(surf_realization)
  !
  ! Sets linkage of flow conditions to dataset
  !
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  !

  use Option_module
  use Patch_module
  use Dataset_Base_class
  use Dataset_module
  use Condition_module

  implicit none

  class(realization_surface_type) :: surf_realization

  type(flow_condition_type), pointer :: cur_surf_flow_condition
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  
  option => surf_realization%option
  
  ! loop over flow conditions looking for linkage to datasets
  cur_surf_flow_condition => surf_realization%surf_flow_conditions%first
  do
    if (.not.associated(cur_surf_flow_condition)) exit
    string = 'flow_condition ' // trim(cur_surf_flow_condition%name)
    ! find datum dataset
    call DatasetFindInList(surf_realization%datasets, &
                           cur_surf_flow_condition%datum, &
                           cur_surf_flow_condition%default_time_storage, &
                           string,option)
    select case(option%iflowmode)
      case(SWE_MODE)
        do i = 1, size(cur_surf_flow_condition%sub_condition_ptr)
           ! find dataset
          call DatasetFindInList(surf_realization%datasets, &
                 cur_surf_flow_condition%sub_condition_ptr(i)%ptr%dataset, &
                 cur_surf_flow_condition%default_time_storage, &
                 string,option)
          ! find gradient dataset
          call DatasetFindInList(surf_realization%datasets, &
                 cur_surf_flow_condition%sub_condition_ptr(i)%ptr%gradient, &
                 cur_surf_flow_condition%default_time_storage, &
                 string,option)
        enddo
      case default
        option%io_buffer='RealizSurfProcessFlowConditions not implemented in this mode'
        call PrintErrMsg(option)
    end select
    cur_surf_flow_condition => cur_surf_flow_condition%next
  enddo

end subroutine RealizationSurfaceProcessFlowConditions

! ************************************************************************** !

subroutine RealizationSurfaceProcessCouplers(surf_realization)
  !
  ! Sets connectivity and pointers for couplers
  !
  ! Author: Gautam Bisht
  ! Date: 04/04/23
  !

  use Patch_module
  use Coupler_module
  use Strata_module
  use Option_module
  use Material_Surface_module

  implicit none

  class(realization_surface_type) :: surf_realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(condition_list_type), pointer :: flow_conditions
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  PetscInt :: temp_int

  patch => surf_realization%patch
  option => surf_realization%option
  flow_conditions => surf_realization%surf_flow_conditions

  ! boundary conditions
  coupler => patch%boundary_condition_list%first
  do
    if (.not.associated(coupler)) exit
    write(*,*)'BC coupler ',trim(coupler%name)
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%region_list)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in boundary condition "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call PrintErrMsg(option)
    endif
    if (associated(patch%grid%structured_grid)) then
      if (coupler%region%num_cells > 0 .and. &
          (coupler%region%iface == 0 .and. &
           .not.associated(coupler%region%faces))) then
        option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '", which is tied to a boundary condition, has not &
                 &been assigned a face in the structured grid. '
        call PrintErrMsg(option)
      endif
    endif
    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name, &
                                      flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in boundary condition "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call PrintErrMsg(option)
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in &
                           &BOUNDARY_CONDITION: ' // trim(coupler%name) // '.'
        call PrintErrMsg(option)
      endif
    endif

    coupler => coupler%next
  enddo


  ! initial conditions
  coupler => patch%initial_condition_list%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%region_list)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in initial condition "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call PrintErrMsg(option)
    endif
    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name, &
                                      flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in initial condition "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call PrintErrMsg(option)
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in ' // &
                           'INITIAL_CONDITION: ' // trim(coupler%name) // '.'
        call PrintErrMsg(option)
      endif
    endif

    coupler => coupler%next
  enddo

  ! source/sinks
  coupler => patch%source_sink_list%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%region_list)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in source/sink "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call PrintErrMsg(option)
    endif

    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name, &
                                      flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in source/sink "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call PrintErrMsg(option)
        endif
        ! check to ensure that a rate subcondition exists
        if (.not.associated(coupler%flow_condition%rate) .and. &
              .not.associated(coupler%flow_condition%well)) then
          temp_int = 0
          if (associated(coupler%flow_condition%general)) then
            if (associated(coupler%flow_condition%general%rate)) then
              temp_int = 1
            endif
          endif
          if (associated(coupler%flow_condition%hydrate)) then
            if (associated(coupler%flow_condition%hydrate%rate)) then
              temp_int = 1
            endif
          endif
          if (temp_int == 0) then
            option%io_buffer = 'FLOW_CONDITIONs associated with &
              &SOURCE_SINKs must have a RATE or WELL expression within them.'
            call PrintErrMsg(option)
          endif
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in &
                           &SOURCE_SINK: ' // trim(coupler%name) // '.'
        call PrintErrMsg(option)
      endif
    endif

    coupler => coupler%next
  enddo

  ! strata
  ! connect pointers from strata to regions
  strata => patch%strata_list%first
  do
    if (.not.associated(strata)) exit
    ! pointer to region
    if (len_trim(strata%region_name) > 0) then
      strata%region => RegionGetPtrFromList(strata%region_name, &
                                                  patch%region_list)
      if (.not.associated(strata%region)) then
        option%io_buffer = 'Region "' // trim(strata%region_name) // &
                 '" in strata not found in region list'
        call PrintErrMsg(option)
      endif
      if (strata%active) then
        ! pointer to material
        strata%material_surface_property => &
          MaterialSurfacePropGetPtrFromArray(strata%material_property_name, &
          patch%surface_material_property_array)
        if (.not.associated(strata%material_surface_property)) then
          option%io_buffer = 'Surface Material "' // &
                            trim(strata%material_property_name) // &
                            '" not found in material list'
          call PrintErrMsg(option)
        endif
      endif
    else
      nullify(strata%region)
      nullify(strata%material_surface_property)
    endif
    strata => strata%next
  enddo

end subroutine RealizationSurfaceProcessCouplers

end module Realization_Surface_class