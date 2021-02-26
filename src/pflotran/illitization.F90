module Illitization_module
  !
  ! Illitization models
  !
  ! Author: Alex Salazar III
  ! Date: 02/25/21
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  !---------------------------------------------------------------------------
  type, public :: illitization_base_type
    PetscReal :: ilt_threshold ! temperature threshold to begin illitization
    PetscReal :: ilt_fs  ! fraction of smectite in material
    PetscReal :: ilt_fi  ! fraction of illite in material
    PetscReal :: ilt_fs0 ! initial fraction of smectite in material
    PetscReal :: ilt_fi0 ! initial fraction of illite in material
    PetscReal :: ilt_ds  ! track accumulated changes in smectite
  contains
    procedure, public :: Verify => ILTBaseVerify
    procedure, public :: Test => ILTBaseTest
    procedure, public :: CalculateILT => ILTBaseIllitization
  end type illitization_base_type
  !---------------------------------------------------------------------------
  type, public, extends(illitization_base_type) :: ILT_default_type
    ! Model by Huang et al., 1993
    PetscReal :: ilt_rate ! temperature-dependent illitization rate in sec^-1
    PetscReal :: ilt_ea   ! activation energy in J/mol
    PetscReal :: ilt_freq ! frequency term in L/mol-sec
    PetscReal :: ilt_K_conc ! molar concentration of potassium
    PetscReal :: ilt_shift_perm ! permeability shift factor for illite fraction
  contains
    procedure, public :: Verify => ILTDefaultVerify
    procedure, public :: CalculateILT => ILTDefaultIllitization
  end type ILT_default_type
  !---------------------------------------------------------------------------
  type, public :: illitization_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
    PetscBool :: test
    class(illitization_base_type), pointer :: illitization_function
    class(illitization_type), pointer :: next
  end type illitization_type
  !---------------------------------------------------------------------------
  type, public :: illitization_ptr_type
    class(illitization_type), pointer :: ptr
  end type illitization_ptr_type
  !---------------------------------------------------------------------------

  public :: IllitizationCreate, &
            IllitizationGetID, &
            IllitizationRead, &
            IllitizationAddToList, &
            IllitizationConvertListToArray, &
            IllitizationInputRecord, &
            IllitizationDestroy, &
            ILTDefaultCreate, &
            ILTAssignDefault

contains

! ************************************************************************** !

subroutine ILTBaseVerify(this,name,option)

  use Option_module

  implicit none

  class(illitization_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  if (Uninitialized(this%ilt_threshold)) then
    option%io_buffer = 'Illitization temperature threshold must be specified ' &
                     //'for function "'//trim(name)//'".'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%ilt_fs0)) then
    option%io_buffer = 'Initial smectite fraction must be specified ' &
                     //'for function "'//trim(name)//'".'
    call PrintErrMsg(option)
  endif

end subroutine ILTBaseVerify

! ************************************************************************** !

subroutine ILTBaseIllitization(this,temperature,illitization,option)

  use Option_module

  implicit none

  class(illitization_base_type) :: this
  PetscReal, intent(in) :: temperature
  PetscReal, intent(out) :: illitization
  type(option_type), intent(inout) :: option

  illitization = 0.d0

  option%io_buffer = 'ILTBaseIllitization must be extended.'
  call PrintErrMsg(option)

end subroutine ILTBaseIllitization

! ************************************************************************** !

subroutine ILTBaseTest(this,ilt_name,option)

  use Option_module

  implicit none

  class(illitization_base_type) :: this
  character(len=MAXWORDLENGTH) :: ilt_name
  type(option_type), intent(inout) :: option

end subroutine ILTBaseTest

! ************************************************************************** !

subroutine ILTDestroy(ilf)

  implicit none

  class(illitization_base_type), pointer :: ilf

  if (.not.associated(ilf)) return
  deallocate(ilf)
  nullify(ilf)

end subroutine ILTDestroy

! ************************************************************************** !

function ILTDefaultCreate()

  implicit none

  class(ILT_default_type), pointer :: ILTDefaultCreate

  allocate(ILTDefaultCreate)

  ILTDefaultCreate%ilt_fs         = 1.0d0
  ILTDefaultCreate%ilt_fi         = 0.0d0
  ILTDefaultCreate%ilt_fs0        = 1.0d0
  ILTDefaultCreate%ilt_fi0        = 0.0d0
  ILTDefaultCreate%ilt_shift_perm = 1.0d0
  ILTDefaultCreate%ilt_threshold  = 0.0d0
  ILTDefaultCreate%ilt_ds         = 0.0d0
  ILTDefaultCreate%ilt_rate   = UNINITIALIZED_DOUBLE
  ILTDefaultCreate%ilt_ea     = UNINITIALIZED_DOUBLE
  ILTDefaultCreate%ilt_freq   = UNINITIALIZED_DOUBLE
  ILTDefaultCreate%ilt_K_conc = UNINITIALIZED_DOUBLE

end function ILTDefaultCreate

! ************************************************************************** !

subroutine ILTDefaultVerify(this,name,option)

  use Option_module

  implicit none

  class(ILT_default_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'ILLITIZATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'ILLITIZATION_FUNCTION, DEFAULT'
  endif
  call ILTBaseVerify(this,string,option)
  if (Uninitialized(this%ilt_ea)) then
    option%io_buffer = 'Illitization activation energy must be specified in' &
                     //' function "'//trim(string)//'".'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%ilt_freq)) then
    option%io_buffer = 'Illitization frequency term must be specified in' &
                     //' function "'//trim(string)//'".'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%ilt_K_conc)) then
    option%io_buffer = 'Illitization potassium concentration must be ' &
                     //'specified in function "'//trim(string)//'".'
    call PrintErrMsg(option)
  endif

end subroutine ILTDefaultVerify

! ************************************************************************** !

subroutine ILTDefaultIllitization(this,temperature,illitization,option)

  use Option_module

  implicit none

  class(ILT_default_type) :: this
  PetscReal, intent(in) :: temperature
  PetscReal, intent(out) :: illitization
  type(option_type), intent(inout) :: option

  PetscReal :: tempreal

  ! Model based on Huang et al., 1993
 
  tempreal = temperature
  ! Check if we are above the temperature threshold for acknowledging effect
  if(tempreal >= this%ilt_threshold) then

    ! Use Kelvin
    tempreal = tempreal + 273.15d0

    ! Illitization rate - Arrhenius-type model from Huang et al., 1993 [L/mol-s]
    this%ilt_rate = this%ilt_freq * &
      exp(-1.0d0 * this%ilt_ea / (IDEAL_GAS_CONSTANT * tempreal))

    ! Modify rate with potassium concentration and initial fraction [1/s]
    this%ilt_rate = this%ilt_rate * (this%ilt_fs0**2) * this%ilt_K_conc

    ! Log accumulated changes in smectite
    this%ilt_ds = this%ilt_ds + this%ilt_rate * option%dt

    ! Change in smectite
    this%ilt_fs = this%ilt_fs / (1.0d0 + this%ilt_ds)
                            ! (1.0d0 + this%ilt_rate * option%dt)

    if (this%ilt_fs > 1.0d0) then
      this%ilt_fs = 1.0d0
    elseif (this%ilt_fs < 0.0d0) then
      this%ilt_fs = 0.0d0
    endif

    ! Fraction illitized
    this%ilt_fi = 1 - this%ilt_fs
    
    illitization = this%ilt_fi

  else
    illitization = this%ilt_fi
  endif

end subroutine ILTDefaultIllitization

! ************************************************************************** !

function IllitizationCreate()

  implicit none

  class(illitization_type), pointer :: IllitizationCreate
  class(illitization_type), pointer :: illitization_function

  allocate(illitization_function)
  illitization_function%name = ''
  illitization_function%print_me = PETSC_FALSE
  illitization_function%test = PETSC_FALSE
  nullify(illitization_function%illitization_function)
  nullify(illitization_function%next)

  IllitizationCreate => illitization_function

end function IllitizationCreate

! ************************************************************************** !

subroutine IllitizationRead(this,input,option)

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(illitization_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string, verify_string
  class(illitization_base_type), pointer :: illitization_function_ptr

  nullify(illitization_function_ptr)

  input%ierr = 0
  error_string = 'ILLITIZATION'
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
    case('ILLITIZATION_FUNCTION')
      call InputReadCard(input,option,word)
      call InputErrorMsg(input,option, &
           'ILLITIZATION_FUNCTION',error_string)
      call StringToUpper(word)
      select case(word)
        case('DEFAULT')
          this%illitization_function => ILTDefaultCreate()
        case default
          call InputKeywordUnrecognized(input,word, &
               'ILLITIZATION_FUNCTION',option)
      end select
      call ILTRead(this%illitization_function,input,option)
    case('TEST')
      this%test = PETSC_TRUE
    case default
      call InputKeywordUnrecognized(input,keyword,'ILLITIZATION',option)
    end select
  enddo
  call InputPopBlock(input,option)

  verify_string = 'ILLITIZATION(' // trim(this%name) // '),'

  if (associated(this%illitization_function)) then
    call this%illitization_function%Verify(verify_string,option)
  else
    option%io_buffer = 'A illitization function has &
         &not been set under ILLITIZATION "' // &
         trim(this%name) // '". An ILLITIZATION_FUNCTION &
         &block must be specified.'
  endif

end subroutine IllitizationRead

! ************************************************************************** !

subroutine ILTAssignDefault(illitization_function,&
                            fs,fi,fs0,fi0,ea,freq,Kc,shift,thresh,option)

  use Option_module

  implicit none

  class(illitization_base_type) :: illitization_function
  PetscReal :: fs,fi,fs0,fi0,ea,freq,Kc,shift,thresh
  type(option_type) :: option

  select type(ilf => illitization_function)
      !------------------------------------------
    class is(ILT_default_type)
      ilf%ilt_fs = fs
      ilf%ilt_fi = fi
      ilf%ilt_fs0 = fs0
      ilf%ilt_fi0 = fi0
      ilf%ilt_ea = ea
      ilf%ilt_ds = 0.0
      ilf%ilt_freq = freq
      ilf%ilt_K_conc = Kc
      ilf%ilt_shift_perm = shift
      ilf%ilt_threshold = thresh
  end select

end subroutine ILTAssignDefault

! ************************************************************************** !

subroutine ILTRead(illitization_function,input,option)
  !
  ! Reads in contents of a ILLITIZATION_FUNCTION block
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(illitization_base_type) :: illitization_function
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string

  input%ierr = 0
  error_string = 'ILLITIZATION_FUNCTION,'
  select type(ilf => illitization_function)
  class is(ILT_default_type)
    error_string = trim(error_string) // 'DEFAULT'
  end select

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select type(ilf => illitization_function)
    !------------------------------------------
    class is(ILT_default_type)
      select case(trim(keyword))
        case('EA')
          ! Activation energy in Arrhenius term
          call InputReadDouble(input,option,ilf%ilt_ea)
          call InputErrorMsg(input,option,'activation energy', &
                             'ILLITIZATION, DEFAULT')
          call InputReadAndConvertUnits(input,ilf%ilt_ea, &
                      'J/mol','ILLITIZATION, DEFAULT, activation energy',option)
        case('FREQ')
          ! Frequency factor (scaling constant of Arrhenius term)
          call InputReadDouble(input,option,ilf%ilt_freq)
          call InputErrorMsg(input,option,'frequency term', &
                             'ILLITIZATION, DEFAULT')
          call InputReadAndConvertUnits(input,ilf%ilt_freq, &
                       'L/s-mol','ILLITIZATION, DEFAULT, frequency term',option)
        case('K_CONC')
          ! Concentration of Potassium cation
          call InputReadDouble(input,option,ilf%ilt_K_conc)
          call InputErrorMsg(input,option,'potassium concentration', &
                             'ILLITIZATION, DEFAULT')
          call InputReadAndConvertUnits(input,ilf%ilt_K_conc,&
                    'M','ILLITIZATION, DEFAULT, potassium concentration',option)
        case('SMECTITE_INITIAL')
          ! Initial fraction of smectite in the smectite/illite mixture
          call InputReadDouble(input,option,ilf%ilt_fs0)
          call InputErrorMsg(input,option,'initial smectite fraction', &
                             'ILLITIZATION, DEFAULT')
        case('SHIFT_PERM')
          ! Factor modifying permeability per fraction illitized
          call InputReadDouble(input,option,ilf%ilt_shift_perm)
          call InputErrorMsg(input,option,'permeability shift factor', &
                             'ILLITIZATION, DEFAULT')
        case default
          call ILTBaseRead(ilf,input,keyword,error_string,'default',option)
      end select
    !------------------------------------------
    class default
      option%io_buffer = 'Read routine not implemented for ' &
           // trim(error_string) // '.'
      call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine ILTRead

! ************************************************************************** !

subroutine ILTBaseRead(ilf,input,keyword,error_string,kind,option)
  !
  ! Reads in contents of ILLITIZATION_FUNCTION block for dervived
  ! types of illitization base class
  !
  use Option_module
  use Input_Aux_module
  use String_module

  class(ILT_default_type) :: ilf
  type(input_type) :: input
  character(len=MAXWORDLENGTH)   :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=*)  :: kind
  type(option_type) :: option

  select case(keyword)
    case('THRESHOLD')
      ! Specifies the temperature threshold for activating illitization
      call InputReadDouble(input,option, &
                           ilf%ilt_threshold)
      call InputErrorMsg(input,option,'temperature threshold', &
                         'ILLITIZATION, DEFAULT')
      call InputReadAndConvertUnits(input,ilf%ilt_threshold, &
                  'C','ILLITIZATION, DEFAULT, temperature threshold',option)
    case default
       call InputKeywordUnrecognized(input,keyword, &
            'illitization function ('//trim(kind)//')',option)
  end select

end subroutine

! ************************************************************************** !

subroutine IllitizationAddToList(new_ilf,list)

  implicit none

  class(illitization_type), pointer :: new_ilf
  class(illitization_type), pointer :: list

  class(illitization_type), pointer :: cur_ilf

  if (associated(list)) then
    cur_ilf => list
    ! loop to end of list
    do
      if (.not.associated(cur_ilf%next)) exit
      cur_ilf => cur_ilf%next
    enddo
    cur_ilf%next => new_ilf
  else
    list => new_ilf
  endif

end subroutine IllitizationAddToList

! ************************************************************************** !

subroutine IllitizationConvertListToArray(list,array,option)

  use String_module
  use Option_module

  implicit none

  class(illitization_type), pointer :: list
  type(illitization_ptr_type), pointer :: array(:)
  type(option_type) :: option

  class(illitization_type), pointer :: cur_ilf
  PetscInt :: count

  count = 0
  cur_ilf => list
  do
    if (.not.associated(cur_ilf)) exit
    count = count + 1
    cur_ilf => cur_ilf%next
  enddo

  if (associated(array)) deallocate(array)
  allocate(array(count))

  count = 0
  cur_ilf => list
  do
    if (.not.associated(cur_ilf)) exit
    count = count + 1
    array(count)%ptr => cur_ilf
    if (cur_ilf%test) then
      call OptionSetBlocking(option,PETSC_FALSE)
      if (option%myrank == option%io_rank) then
        if (associated(cur_ilf%illitization_function)) then
          ! call cur_ilf%illitization_function%Test( &
          !      cur_ilf%name,option)
        endif
      endif
      call OptionSetBlocking(option,PETSC_TRUE)
      call OptionCheckNonBlockingError(option)
    endif
    cur_ilf => cur_ilf%next
  enddo

end subroutine IllitizationConvertListToArray

! ************************************************************************** !

function IllitizationGetID(illitization_array, &
     illitization_name, material_property_name, option)

  use Option_module
  use String_module

  type(illitization_ptr_type), pointer :: illitization_array(:)
  character(len=MAXWORDLENGTH) :: illitization_name
  character(len=MAXWORDLENGTH) :: test1, test2
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: iid, IllitizationGetID
  PetscInt :: i, j

  do i = 1, size(illitization_array)
      test1 = illitization_array(i)%ptr%name
      do j = 1, size(illitization_array)
        if (i == j) cycle
        test2 = illitization_array(j)%ptr%name
        if (test1 == test2) then
          option%io_buffer = 'Duplicate illitization function '//&
                             trim(test2)//&
                             ' has been detected.'
          call PrintErrMsg(option)
        endif
      enddo
  enddo

  IllitizationGetID = 0
  do iid = 1, size(illitization_array)
    if (StringCompare(illitization_name,illitization_array(iid)%ptr%name)) then
      IllitizationGetID = iid
      return
    endif
  enddo
  
  ! IllitizationGetID = UNINITIALIZED_INTEGER
  option%io_buffer = 'Illitization function "' // &
       trim(illitization_name) // &
       '" in material property "' // &
       trim(material_property_name) // &
       '" not found among available illitization functions.'
  call PrintErrMsg(option)

end function IllitizationGetID

! ************************************************************************** !

subroutine IllitizationInputRecord(illitization_list)

  implicit none

  class(illitization_type), pointer :: illitization_list

  class(illitization_type), pointer :: cur_ilf
  character(len=MAXWORDLENGTH) :: word1
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
       &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'ILLITIZATION FUNCTIONS'

  cur_ilf => illitization_list
  do
    if (.not.associated(cur_ilf)) exit

    write(id,'(a29)',advance='no') 'illitization function name: '
    write(id,'(a)') adjustl(trim(cur_ilf%name))

    if (associated(cur_ilf%illitization_function)) then
      write(id,'(a29)',advance='no') 'model: '
      select type (ilf => cur_ilf%illitization_function)
        !---------------------------------
      class is (ILT_default_type)
        write(id,'(a)') 'Huang et al., 1993'
        write(id,'(a29)',advance='no') 'frequency: '
        write(word1,'(es12.5)') ilf%ilt_freq
        write(id,'(a)') adjustl(trim(word1))//' L/mol-s'
        write(id,'(a29)',advance='no') 'activation energy: '
        write(word1,'(es12.5)') ilf%ilt_ea
        write(id,'(a)') adjustl(trim(word1))//' J/mol'
        write(id,'(a29)',advance='no') 'K+ concentration: '
        write(word1,'(es12.5)') ilf%ilt_K_conc
        write(id,'(a)') adjustl(trim(word1))//' M'
        write(id,'(a29)',advance='no') 'initial smectite: '
        write(word1,'(es12.5)') ilf%ilt_fs0
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'temperature threshold: '
        write(word1,'(es12.5)') ilf%ilt_threshold
        write(id,'(a)') adjustl(trim(word1))//' C'
        write(id,'(a29)',advance='no') 'shift (permeability): '
        write(word1,'(es12.5)') ilf%ilt_shift_perm
        write(id,'(a)') adjustl(trim(word1))
      end select
    endif

    write(id,'(a29)') '---------------------------: '
    cur_ilf => cur_ilf%next
  enddo

end subroutine IllitizationInputRecord

! ************************************************************************** !

recursive subroutine IllitizationDestroy(ilt)

  implicit none

  class(illitization_type), pointer :: ilt

  if (.not.associated(ilt)) return

  call IllitizationDestroy(ilt%next)

  call ILTDestroy(ilt%illitization_function)

  deallocate(ilt)
  nullify(ilt)

end subroutine IllitizationDestroy

end module Illitization_module
