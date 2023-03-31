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

  implicit none

  private

  type, public, extends(realization_base_type) :: realization_surface_type

    type(region_list_type), pointer :: surf_regions
    type(condition_list_type),pointer :: surf_flow_conditions

    type(field_surface_type), pointer :: field_surface

    character(len=MAXSTRINGLENGTH) :: surf_filename
    character(len=MAXSTRINGLENGTH) :: subsurf_filename

    class(dataset_base_type), pointer :: datasets
  
  end type realization_surface_type

  public :: RealizationSurfaceCreate, &
            RealizationSurfaceCreateDiscretization

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

  allocate(surf_realization%surf_regions)
  call RegionInitList(surf_realization%surf_regions)
  
  allocate(surf_realization%surf_flow_conditions)
  call FlowConditionInitList(surf_realization%surf_flow_conditions)

  surf_realization%field_surface => FieldSurfaceCreate()

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
  type(dm_ptr_type), pointer :: dm_ptr
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

end module Realization_Surface_class