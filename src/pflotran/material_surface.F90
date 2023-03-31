module Material_Surface_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
  type, public :: material_surface_property_type
    
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: external_id
    PetscInt :: internal_id
    PetscReal :: mannings
    
    type(material_surface_property_type), pointer :: next
  end type material_surface_property_type
  
  type, public :: surface_material_property_ptr_type
    type(material_surface_property_type), pointer :: ptr
  end type surface_material_property_ptr_type

  public :: MaterialSurfacePropertyCreate, &
            MaterialSurfacePropertyDestroy, &
            MaterialSurfacePropertyAddToList, &
            MaterialSurfacePropertyRead, &
            MaterialSurfacePropConvertListToArray, &
            MaterialSurfacePropGetPtrFromArray, &
            MaterialSurfaceGetMaxExternalID, &
            MaterialSurfaceCreateIntToExtMapping, &
            MaterialSurfaceCreateExtToIntMapping, &
            MaterialSurfaceApplyMapping

  contains

! ************************************************************************** !

function MaterialSurfacePropertyCreate()
  ! 
  ! This routine creates a surface material property
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/09/12
  ! 

  implicit none
  
  type(material_surface_property_type), pointer :: MaterialSurfacePropertyCreate
  type(material_surface_property_type), pointer :: material_surface_property
  
  allocate(material_surface_property)

  material_surface_property%name        = ''
  material_surface_property%internal_id = 0
  material_surface_property%external_id = 0
  material_surface_property%mannings    = 0.d0
  
  nullify(material_surface_property%next)
  
  MaterialSurfacePropertyCreate => material_surface_property

end function MaterialSurfacePropertyCreate

! ************************************************************************** !

subroutine MaterialSurfacePropertyRead(material_surface_property,input,option)
  ! 
  ! This routine reads in contents of a surface material property
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/09/12
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  
  implicit none
  
  type(material_surface_property_type) :: material_surface_property
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: string

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    
    if (InputCheckExit(input,option)) exit
  
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword','SURFACE_MATERIAL_PROPERTY')
    call StringToUpper(keyword)
    
    select case(trim(keyword))
      case('ID')
        call InputReadInt(input,option,material_surface_property%external_id)
        call InputErrorMsg(input,option,'id','SURFACE_MATERIAL_PROPERTY')
      case('MANNINGS_N')
        call InputReadDouble(input,option,material_surface_property%mannings)
        call InputErrorMsg(input,option,'MANNINGS','SURFACE_MATERIAL_PROPERTY')
      case default
        call InputKeywordUnrecognized(input,keyword, &
                                      'SURFACE_MATERIAL_PROPERTY',option)
      end select
  enddo
  call InputPopBlock(input,option)
  
end subroutine MaterialSurfacePropertyRead

! ************************************************************************** !

subroutine MaterialSurfacePropertyAddToList(material_surface_property,list)
  ! 
  ! This routine adds a surface material property to a linked list
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/09/12
  ! 

  implicit none
  
  type(material_surface_property_type), pointer :: material_surface_property
  type(material_surface_property_type), pointer :: list
  type(material_surface_property_type), pointer :: cur_material_surface_property
  
  if (associated(list)) then
    cur_material_surface_property => list
    ! loop to end of list
    do
      if (.not.associated(cur_material_surface_property%next)) exit
      cur_material_surface_property => cur_material_surface_property%next
    enddo
    cur_material_surface_property%next => material_surface_property
    material_surface_property%internal_id = cur_material_surface_property%internal_id + 1
  else
    list => material_surface_property
    material_surface_property%internal_id = 1
  endif
  
end subroutine MaterialSurfacePropertyAddToList

! ************************************************************************** !

recursive subroutine MaterialSurfacePropertyDestroy(material_surface_property)
  ! 
  ! This routine destroys a surface material property
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/09/12
  ! 

  implicit none
  
  type(material_surface_property_type), pointer :: material_surface_property
  
  if (.not.associated(material_surface_property)) return
  
  call MaterialSurfacePropertyDestroy(material_surface_property%next)
  
  deallocate(material_surface_property)
  nullify(material_surface_property)
  
end subroutine MaterialSurfacePropertyDestroy

! ************************************************************************** !

subroutine MaterialSurfacePropConvertListToArray(list,array,option)
  ! 
  ! This routine creates an array of pointers to the surface_material_properties
  ! in the list (similar to subroutine MaterialPropConvertListToArray)
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/11/12
  ! 
  use Option_module
  use String_module

  implicit none

  type(material_surface_property_type), pointer :: list
  type(surface_material_property_ptr_type), pointer :: array(:)
  type(option_type) :: option

  type(material_surface_property_type), pointer :: cur_material_property
  type(material_surface_property_type), pointer :: prev_material_property
  type(material_surface_property_type), pointer :: next_material_property
  PetscInt :: i, j, length1,length2, max_internal_id, max_external_id
  PetscInt, allocatable :: id_count(:)
  PetscBool :: error_flag
  character(len=MAXSTRINGLENGTH) :: string

  ! check to ensure that max internal id is equal to the number of
  ! material properties and that internal ids are contiguous
  max_internal_id = 0
  max_external_id = 0
  cur_material_property => list
  do
    if (.not.associated(cur_material_property)) exit
    max_internal_id = max_internal_id + 1
    max_external_id = max(max_external_id,cur_material_property%external_id)
    if (max_internal_id /= cur_material_property%internal_id) then
      write(string,*) cur_material_property%external_id
      option%io_buffer = 'Non-contiguous internal material id for ' // &
        'material named "' // trim(cur_material_property%name) // &
        '" with external id "' // trim(adjustl(string)) // '" '
      write(string,*) cur_material_property%internal_id
      option%io_buffer = trim(option%io_buffer) // &
        'and internal id "' // trim(adjustl(string)) // '".'
      call PrintErrMsg(option)
    endif
    cur_material_property => cur_material_property%next
  enddo

  allocate(array(max_internal_id))
  do i = 1, max_internal_id
    nullify(array(i)%ptr)
  enddo

  ! use id_count to ensure that an id is not duplicated
  allocate(id_count(max_external_id))
  id_count = 0

  cur_material_property => list
  do
    if (.not.associated(cur_material_property)) exit
    id_count(cur_material_property%external_id) = &
      id_count(cur_material_property%external_id) + 1
    array(cur_material_property%internal_id)%ptr => cur_material_property
    cur_material_property => cur_material_property%next
  enddo

  ! check to ensure that an id is not duplicated
  error_flag = PETSC_FALSE
  do i = 1, max_external_id
    if (id_count(i) > 1) then
      write(string,*) i
      option%io_buffer = 'Material ID ' // trim(adjustl(string)) // &
        ' is duplicated in input file.'
      call PrintMsg(option)
      error_flag = PETSC_TRUE
    endif
  enddo

  deallocate(id_count)

  if (error_flag) then
    option%io_buffer = 'Duplicate Material IDs.'
    call PrintErrMsg(option)
  endif

  ! ensure unique material names
  error_flag = PETSC_FALSE
  do i = 1, size(array)
    if (associated(array(i)%ptr)) then
      length1 = len_trim(array(i)%ptr%name)
      do j = 1, i-1
        if (associated(array(j)%ptr)) then
          length2 = len_trim(array(j)%ptr%name)
          if (length1 /= length2) cycle
          if (StringCompare(array(i)%ptr%name,array(j)%ptr%name,length1)) then
            option%io_buffer = 'Material name "' // &
              trim(adjustl(array(i)%ptr%name)) // &
              '" is duplicated in input file.'
            call PrintMsg(option)
            error_flag = PETSC_TRUE
          endif
        endif
      enddo
    endif
  enddo

  if (error_flag) then
    option%io_buffer = 'Duplicate Material names.'
    call PrintErrMsg(option)
  endif

end subroutine MaterialSurfacePropConvertListToArray

! ************************************************************************** !

function MaterialSurfacePropGetPtrFromArray(material_surface_property_name, &
                                            material_surface_property_array)
  ! 
  ! This routine returns a pointer to the surface material property matching
  ! surface_material_propertry_name (similar to subroutine
  ! MaterialPropGetPtrFromArray)
  ! 
  ! Author: Gautam Bisht, ORNL
  ! Date: 02/11/12
  ! 

  use String_module

  implicit none

  type(material_surface_property_type), pointer :: MaterialSurfacePropGetPtrFromArray
  type(surface_material_property_ptr_type), pointer :: material_surface_property_array(:)
  character(len=MAXWORDLENGTH) :: material_surface_property_name
  PetscInt :: length
  PetscInt :: imaterial_surface_property

  nullify(MaterialSurfacePropGetPtrFromArray)

  do imaterial_surface_property = 1, size(material_surface_property_array)
    length = len_trim(material_surface_property_name)
    if (.not.associated(material_surface_property_array(imaterial_surface_property)%ptr)) cycle
    if (length == &
        len_trim(material_surface_property_array(imaterial_surface_property)%ptr%name) .and. &
        StringCompare(material_surface_property_array(imaterial_surface_property)%ptr%name, &
                        material_surface_property_name,length)) then
      MaterialSurfacePropGetPtrFromArray => &
        material_surface_property_array(imaterial_surface_property)%ptr
      return
    endif
  enddo

end function MaterialSurfacePropGetPtrFromArray

! ************************************************************************** !

function MaterialSurfaceGetMaxExternalID(material_surface_property_array)
  !
  ! Maps internal material ids to external for I/O, etc. [copy of
  ! MaterialGetMaxExternalID()]
  !
  ! Author: Gautam Bisht
  ! Date: 08/05/14
  !
  implicit none

  type(surface_material_property_ptr_type) :: material_surface_property_array(:)

  PetscInt :: MaterialSurfaceGetMaxExternalID

  PetscInt :: i

  MaterialSurfaceGetMaxExternalID = UNINITIALIZED_INTEGER
  do i = 1, size(material_surface_property_array)
    MaterialSurfaceGetMaxExternalID = max(MaterialSurfaceGetMaxExternalID, &
                                         (material_surface_property_array(i)%ptr%external_id))
  enddo

end function MaterialSurfaceGetMaxExternalID

! ************************************************************************** !

subroutine MaterialSurfaceCreateIntToExtMapping(material_surface_property_array,mapping)
  !
  ! Maps internal material ids to external for I/O, etc. [copy of
  ! MaterialCreateIntToExtMapping()]
  !
  ! Author: Gautam Bisht.
  ! Date: 08/08/14
  !
  implicit none

  type(surface_material_property_ptr_type) :: material_surface_property_array(:)
  PetscInt, pointer :: mapping(:)

  PetscInt :: i

  allocate(mapping(size(material_surface_property_array)))
  mapping = UNINITIALIZED_INTEGER

  do i = 1, size(material_surface_property_array)
    mapping(material_surface_property_array(i)%ptr%internal_id) = &
      material_surface_property_array(i)%ptr%external_id
  enddo

end subroutine MaterialSurfaceCreateIntToExtMapping

! ************************************************************************** !

subroutine MaterialSurfaceCreateExtToIntMapping(material_surface_property_array,mapping)
  !
  ! Maps external material ids to internal for setup. This array should be
  ! temporary and never stored for the duration of the simulation.
  ! [copy of MaterialCreateExtToIntMapping()]
  !
  ! Author: Gautam Bisht
  ! Date: 08/08/14
  !
  implicit none

  type(surface_material_property_ptr_type) :: material_surface_property_array(:)
  PetscInt, pointer :: mapping(:)

  PetscInt :: i

  allocate(mapping(MaterialSurfaceGetMaxExternalID(material_surface_property_array)))
  mapping = -888

  do i = 1, size(material_surface_property_array)
    mapping(material_surface_property_array(i)%ptr%external_id) = &
      material_surface_property_array(i)%ptr%internal_id
  enddo

end subroutine MaterialSurfaceCreateExtToIntMapping

! ************************************************************************** !

subroutine MaterialSurfaceApplyMapping(mapping,array)
  !
  ! Maps internal material ids to external for I/O, etc. [copy of
  ! MaterialApplyMapping()]
  !
  ! Author: Gautam Bisht
  ! Date: 08/08/14
  !
  implicit none

  PetscInt :: mapping(:)
  PetscInt :: array(:)

  PetscInt :: i
  PetscInt :: mapping_size
  PetscInt :: mapped_id

  mapping_size = size(mapping)
  do i = 1, size(array)
    if (array(i) <= mapping_size) then
      mapped_id = mapping(array(i))
    else
      mapped_id = -888 ! indicates corresponding mapped value does not exist.
    endif
    if (mapped_id == -888) then ! negate material id to indicate not found
      mapped_id = -1*array(i)
    endif
    array(i) = mapped_id
  enddo

end subroutine MaterialSurfaceApplyMapping

end module Material_Surface_module
