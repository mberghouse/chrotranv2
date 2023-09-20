module Init_Surface_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  public :: SurfaceInitMaterialProperties

contains

! ************************************************************************** !

subroutine SurfaceInitMaterialProperties(surface_realization)
  !
  ! Sets connectivity and pointers for couplers
  !
  ! Author: Gautam Bisht
  ! Date: 04/06/23
  !

  use petscvec
  use Realization_Surface_class
  use Discretization_module
  use Strata_module
  use Region_module
  use Material_module
  use Option_module
  use Grid_module
  use Field_module
  use Patch_module
  use Field_Surface_module
  use Material_Surface_module
  use HDF5_module

  implicit none

  class(realization_surface_type) :: surface_realization

  PetscReal, pointer :: man0_p(:)

  PetscInt :: icell, local_id, ghosted_id, surf_material_id
  PetscInt :: istart, iend
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(field_surface_type), pointer :: field_surface
  type(strata_type), pointer :: strata
  type(patch_type), pointer :: cur_patch

  type(material_surface_property_type), pointer :: material_surface_property
  type(material_surface_property_type), pointer :: null_material_surface_property
  type(region_type), pointer :: region
  PetscBool :: update_ghosted_material_ids

  option => surface_realization%option
  discretization => surface_realization%discretization
  field_surface => surface_realization%field_surface

  ! loop over all patches and allocation material id arrays
  cur_patch => surface_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    if (.not.associated(cur_patch%imat)) then
      allocate(cur_patch%imat(cur_patch%grid%ngmax))
      ! initialize to "unset"
      cur_patch%imat = UNINITIALIZED_INTEGER
      ! also allocate saturation function id
      allocate(cur_patch%cc_id(cur_patch%grid%ngmax))
      cur_patch%cc_id = UNINITIALIZED_INTEGER
      allocate(cur_patch%cct_id(cur_patch%grid%ngmax))
      cur_patch%cct_id = UNINITIALIZED_INTEGER
    endif
    cur_patch => cur_patch%next
  enddo

  ! if material ids are set based on region, as opposed to being read in
  ! we must communicate the ghosted ids.  This flag toggles this operation.
  update_ghosted_material_ids = PETSC_FALSE
  cur_patch => surface_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    grid => cur_patch%grid
    strata => cur_patch%strata_list%first
    do
      if (.not.associated(strata)) exit
      ! Read in cell by cell material ids if they exist
      if (.not.associated(strata%region) .and. strata%active) then
        option%io_buffer = 'Reading of material prop from file for' // &
          ' surface flow is not implemented.'
        call PrintErrMsgByRank(option)
        !call readMaterialsFromFile(realization,strata%realization_dependent, &
        !                           strata%material_property_filename)
      ! Otherwise, set based on region
      else if (strata%active) then
        update_ghosted_material_ids = PETSC_TRUE
        region => strata%region
        material_surface_property => strata%material_surface_property
        if (associated(region)) then
          istart = 1
          iend = region%num_cells
        else
          istart = 1
          iend = grid%nlmax
        endif
        do icell=istart, iend
          if (associated(region)) then
            local_id = region%cell_ids(icell)
          else
            local_id = icell
          endif
          ghosted_id = grid%nL2G(local_id)
          cur_patch%imat(ghosted_id) = material_surface_property%internal_id
        enddo
      endif
      strata => strata%next
    enddo
    cur_patch => cur_patch%next
  enddo

  if (update_ghosted_material_ids) then
    ! update ghosted material ids
    call RealizationSurfaceLocalToLocalWithArray(surface_realization,MATERIAL_ID_ARRAY)
  endif

  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_material_surface_property => MaterialSurfacePropertyCreate()
  cur_patch => surface_realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    call VecGetArrayF90(field_surface%mannings0,man0_p,ierr);CHKERRQ(ierr)

    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      surf_material_id = cur_patch%imat(ghosted_id)
      if (surf_material_id == 0) then ! accomodate inactive cells
        material_surface_property = null_material_surface_property
      else if ( surf_material_id > 0 .and. &
                surf_material_id <= &
                size(cur_patch%surface_material_property_array)) then
        material_surface_property => &
          cur_patch%surface_material_property_array(surf_material_id)%ptr
        if (.not.associated(material_surface_property)) then
          write(dataset_name,*) surf_material_id
          option%io_buffer = 'No material property for surface material id ' // &
                              trim(adjustl(dataset_name)) &
                              //  ' defined in input file.'
          call PrintErrMsgByRank(option)
        endif
      else if (Uninitialized(surf_material_id)) then
        write(dataset_name,*) grid%nG2A(ghosted_id)
        option%io_buffer = 'Uninitialized surface material id in patch at cell ' // &
                            trim(adjustl(dataset_name))
        call PrintErrMsgByRank(option)
      else if (surf_material_id > size(cur_patch%surface_material_property_array)) then
        write(option%io_buffer,*) surf_material_id
        option%io_buffer = 'Unmatched surface material id in patch:' // &
          adjustl(trim(option%io_buffer))
        call PrintErrMsgByRank(option)
      else
        option%io_buffer = 'Something messed up with surface material ids. ' // &
          ' Possibly material ids not assigned to all grid cells. ' // &
          ' Contact Glenn!'
        call PrintErrMsgByRank(option)
      endif
      man0_p(local_id) = material_surface_property%mannings
    enddo ! local_id - loop

    call VecRestoreArrayF90(field_surface%mannings0,man0_p,ierr);CHKERRQ(ierr)

    cur_patch => cur_patch%next
  enddo ! looping over patches

  call MaterialSurfacePropertyDestroy(null_material_surface_property)
  nullify(null_material_surface_property)

  call DiscretizationGlobalToLocal(discretization,field_surface%mannings0, &
                                   field_surface%mannings_loc,ONEDOF)


end subroutine SurfaceInitMaterialProperties

end module Init_Surface_module