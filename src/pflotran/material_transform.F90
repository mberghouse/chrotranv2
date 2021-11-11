module Material_Transform_module
  !
  ! Models to transform material properties
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
    PetscReal :: ilt_fs0 ! initial fraction of smectite in material
    PetscReal :: ilt_shift_perm ! permeability shift factor for illite fraction
    class(ilt_kd_effects_type), pointer :: ilt_shift_kd_list
  contains
    procedure, public :: Verify => ILTBaseVerify
    procedure, public :: Test => ILTBaseTest
    procedure, public :: CalculateILT => ILTBaseIllitization
    procedure, public :: ShiftKd => ILTBaseShiftSorption
    procedure, public :: CheckElements => ILTBaseCheckElements
  end type illitization_base_type
  !---------------------------------------------------------------------------
  type, public, extends(illitization_base_type) :: ILT_default_type
    ! Model by W-L Huang, J.M. Longo, & D.R. Pevear, 1993
    PetscReal :: ilt_ea   ! activation energy in J/mol
    PetscReal :: ilt_freq ! frequency term (in L/mol-sec for default)
    PetscReal :: ilt_K_conc ! molar concentration of potassium
  contains
    procedure, public :: Verify => ILTDefaultVerify
    procedure, public :: CalculateILT => ILTDefaultIllitization
    procedure, public :: ShiftKd => ILTShiftSorption
    procedure, public :: CheckElements => ILTCheckElements
  end type ILT_default_type
  !---------------------------------------------------------------------------
  type, public, extends(ILT_default_type) :: ILT_general_type
    ! Generalized model by J. Cuadros and J. Linares, 1996
    PetscReal :: ilt_K_exp ! exponent of potassium concentration
    PetscReal :: ilt_exp   ! exponent of smectite fraction
  contains
    procedure, public :: Verify => ILTGeneralVerify
    procedure, public :: CalculateILT => ILTGeneralIllitization
  end type ILT_general_type
  !---------------------------------------------------------------------------
  type, public :: material_transform_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
    PetscBool :: test
    class(illitization_base_type), pointer :: illitization_function
    class(material_transform_type), pointer :: next
  end type material_transform_type
  !---------------------------------------------------------------------------
  type, public :: material_transform_ptr_type
    class(material_transform_type), pointer :: ptr
  end type material_transform_ptr_type
  !---------------------------------------------------------------------------
  type :: ilt_kd_effects_type
    PetscInt :: num_elements
    PetscReal, pointer :: f_kd(:,:) ! factors for modifying the kd value
    character(len=MAXWORDLENGTH), pointer :: f_kd_mode(:) ! function type
    character(len=MAXWORDLENGTH), pointer :: f_kd_element(:) ! element affected
  end type
  !---------------------------------------------------------------------------

  public :: MaterialTransformCreate, &
            MaterialTransformGetID, &
            MaterialTransformCheckILT, &
            MaterialTransformAddToList, &
            MaterialTransformConvertListToArray, &
            MaterialTransformDestroy, &
            MaterialTransformInputRecord, &
            MaterialTransformRead, &
            ILTBaseCreate, &
            ILTDefaultCreate

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
  else
    if (this%ilt_fs0 <= 0.0d0 .or. this%ilt_fs0 > 1.0d0) then
      option%io_buffer = 'Initial smectite fraction for function "' &
        //trim(name)//'" must be nonzero positive number up to 1.'
      call PrintErrMsg(option)
    endif
  endif
  if (Initialized(this%ilt_shift_perm)) then
    if (this%ilt_shift_perm <= -1.0d0) then
      option%io_buffer = 'If negative, the illitization permeability shift ' &
                       //'factor must be greater than -1 ' &
                       //'in function "'//trim(name)//'".'
      call PrintErrMsg(option)
    endif
  endif

end subroutine ILTBaseVerify

! ************************************************************************** !

subroutine ILTBaseIllitization(this,fs,temperature,dt,fi,scale,option)

  use Option_module

  implicit none

  class(illitization_base_type) :: this
  PetscReal, intent(inout) :: fs
  PetscReal, intent(in) :: temperature
  PetscReal, intent(in) :: dt
  PetscReal, intent(out) :: fi
  PetscReal, intent(out) :: scale
  type(option_type), intent(inout) :: option

  fi = 0.0d+0
  scale = 0.0d+0

end subroutine ILTBaseIllitization

! ************************************************************************** !

subroutine ILTBaseTest(this,name,option)

  use Option_module

  implicit none

  class(illitization_base_type) :: this
  character(len=MAXWORDLENGTH) :: name
  type(option_type), intent(inout) :: option

  ! Test with pertubrations to initial smectite and temperature over time
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: ns = 10
  PetscInt, parameter :: nt = 61
  PetscInt, parameter :: np = 27
  PetscReal, parameter :: perturbation = 1.0d-6
  PetscReal :: deltaSmec
  PetscReal :: deltaTemp
  PetscReal :: smec_vec(ns)
  PetscReal :: temp_vec(nt)
  PetscReal :: time_vec(np)
  PetscReal :: fi(ns,nt,np)
  PetscReal :: dfi_dtemp_numerical(ns,nt,np)
  PetscReal :: perturbed_temp
  PetscReal :: fi_temp_pert
  PetscReal :: smec_min, smec_max
  PetscReal :: temp_min, temp_max
  PetscReal :: dt,scale,fs0_original
  PetscReal :: fs, fsp
  PetscInt :: i,j,k

  ! thermal conductivity as a function of temp. and liq. sat.
  smec_min = 1.0d-1 ! Minimum fraction smectite
  smec_max = 1.0d+0 ! Maximum fraction smectite
  temp_min = 2.0d+1 ! Minimum temperature in Celcius
  temp_max = 2.6d+2 ! Maximum temperature in Celcius

  deltaSmec = (smec_max - smec_min)/(ns - 1)
  deltaTemp = (temp_max - temp_min)/(nt - 1)

  smec_vec = [(smec_min + i*deltaSmec, i=0,ns-1)]
  temp_vec = [(temp_min + i*deltaTemp, i=0,nt-1)]
  time_vec = (/0.,1.,2.5,5.,7.5,10.,25.,50.,75.,100.,250.,500.,750.,1000., &
               2500.,5000.,7500.,10000.,20000.,30000.,40000.,50000.,60000.,&
               70000.,80000.,90000.,100000./)

  fs0_original = this%ilt_fs0

  do i = 1,ns
    do j = 1,nt
      ! reset base variables to initial
      this%ilt_fs0 = smec_vec(i)
      fs  = smec_vec(i)
      do k = 2,np

        ! get change in time
        dt = time_vec(k) - time_vec(k-1) ! years
        dt = dt*(365.25*24*60*60)        ! convert to seconds

        ! base case with analytical derivatives
        fsp = fs
        call this%CalculateILT(fs,temp_vec(j),dt,fi(i,j,k),scale,option)

        ! calculate numerical derivatives via finite differences
        perturbed_temp = temp_vec(j) * (1.d0 + perturbation)
        call this%CalculateILT(fsp,perturbed_temp,dt,fi_temp_pert,scale,option)

        dfi_dtemp_numerical(i,j,k) = (fi_temp_pert - fi(i,j,k))/ &
                                      (temp_vec(j)*perturbation)
      enddo
    enddo
  enddo

  write(string,*) name
  string = trim(name) // '_ilt_vs_time_and_temp.dat'
  open(unit=86,file=string)
  write(86,*) '"initial smectite [-]", "temperature [C]", &
               "time [yr]", "illite [-]", "dillite/dT [1/yr]"'
  do i = 1,ns
    do j = 1,nt
      do k = 2,np
        write(86,'(5(ES14.6))') smec_vec(i), temp_vec(j), time_vec(k), &
             fi(i,j,k),dfi_dtemp_numerical(i,j,k)
      enddo
    enddo
  enddo
  close(86)

  ! reset to original values
  this%ilt_fs0 = fs0_original

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

function ILTBaseCreate()

  implicit none

  class(illitization_base_type), pointer :: ILTBaseCreate

  allocate(ILTBaseCreate)

  ILTBaseCreate%ilt_threshold  = 0.0d0
  ILTBaseCreate%ilt_fs0        = 0.0d0
  ILTBaseCreate%ilt_shift_perm = 0.0d0
  nullify(ILTBaseCreate%ilt_shift_kd_list)

end function ILTBaseCreate

! ************************************************************************** !

function ILTDefaultCreate()

  implicit none

  class(ILT_default_type), pointer :: ILTDefaultCreate

  allocate(ILTDefaultCreate)

  ILTDefaultCreate%ilt_threshold  = 0.0d0
  ILTDefaultCreate%ilt_fs0        = 1.0d0
  ILTDefaultCreate%ilt_shift_perm = 0.0d0
  ILTDefaultCreate%ilt_ea     = UNINITIALIZED_DOUBLE
  ILTDefaultCreate%ilt_freq   = UNINITIALIZED_DOUBLE
  ILTDefaultCreate%ilt_K_conc = UNINITIALIZED_DOUBLE
  nullify(ILTDefaultCreate%ilt_shift_kd_list)

end function ILTDefaultCreate

! ************************************************************************** !

function ILTGeneralCreate()

  implicit none

  class(ILT_general_type), pointer :: ILTGeneralCreate

  allocate(ILTGeneralCreate)

  ILTGeneralCreate%ilt_threshold  = 0.0d0
  ILTGeneralCreate%ilt_fs0        = 1.0d0
  ILTGeneralCreate%ilt_shift_perm = 0.0d0
  ILTGeneralCreate%ilt_freq       = 1.0d0 ! Default of 1.0 in general model
  ILTGeneralCreate%ilt_ea     = UNINITIALIZED_DOUBLE
  ILTGeneralCreate%ilt_K_conc = UNINITIALIZED_DOUBLE
  ILTGeneralCreate%ilt_K_exp  = UNINITIALIZED_DOUBLE
  ILTGeneralCreate%ilt_exp    = UNINITIALIZED_DOUBLE
  nullify(ILTGeneralCreate%ilt_shift_kd_list)

end function ILTGeneralCreate

! ************************************************************************** !

function ILTKdEffectsCreate()
  ! 
  ! Creates object for modifying Kd values from smectite/illite transition
  ! 
  ! Author: Alex Salazar III
  ! Date: 10/06/2021

  implicit none

  class(ilt_kd_effects_type), pointer :: ILTKdEffectsCreate

  allocate(ILTKdEffectsCreate)
  
  ILTKdEffectsCreate%num_elements = UNINITIALIZED_INTEGER
  nullify(ILTKdEffectsCreate%f_kd)
  nullify(ILTKdEffectsCreate%f_kd_element)

end function ILTKdEffectsCreate

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

subroutine ILTGeneralVerify(this,name,option)

  use Option_module

  implicit none

  class(ILT_general_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'ILLITIZATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'ILLITIZATION_FUNCTION, GENERAL'
  endif
  call ILTBaseVerify(this,string,option)
  call ILTDefaultVerify(this,string,option)
  if (Uninitialized(this%ilt_exp)) then
    option%io_buffer = 'Illitization smectite exponent must be specified in' &
                     //' function "'//trim(string)//'".'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%ilt_K_exp)) then
    option%io_buffer = 'Illitization postassium exponent must be specified in' &
                     //' function "'//trim(string)//'".'
    call PrintErrMsg(option)
  endif

end subroutine ILTGeneralVerify

! ************************************************************************** !

subroutine ILTDefaultIllitization(this,fs,temperature,dt,fi,scale,option)

  use Option_module

  implicit none

  class(ILT_default_type) :: this
  PetscReal, intent(inout) :: fs
  PetscReal, intent(in) :: temperature
  PetscReal, intent(in) :: dt
  PetscReal, intent(out) :: fi
  PetscReal, intent(out) :: scale
  type(option_type), intent(inout) :: option

  PetscReal :: ds   ! change in smectite
  PetscReal :: T    ! temperature in Kelvin
  PetscReal :: rate ! temperature-dependent illitization rate in sec^-1

  ! Model based on W-L Huang, J.M. Longo, & D.R. Pevear,
  !  "An Experimentally Derived Kinetic Model for Smectite-to-Illite Conversion
  !  and its Use as a Geothermometer," Clay and Clay Minerals, vol. 41, no 2., 
  !  pp. 162-177, 1993

  ! Use Kelvin to calculate rate
  T = temperature + 273.15d0

  ! Check if temperature is above threshold for illitization
  if (temperature >= this%ilt_threshold) then
    ! Negative of illitization rate [L/mol-s]
    rate = this%ilt_K_conc * this%ilt_freq * &
      exp(-1.0d0 * this%ilt_ea / (IDEAL_GAS_CONSTANT * T))
  else
    rate = 0.0d0
  endif

  ! Log change in smectite as time proceeds
  ds = rate * dt

  ! Fraction smectite
  fs = fs / (1.0d0 + (fs * ds))

  if (fs > 1.0d0) then
    fs = 1.0d0
  elseif (fs < 0.0d0) then
    fs = 0.0d0
  endif

  ! Fraction illite
  fi = 1.0d0 - fs
  
  ! Calculate scale factor
  scale = ((fi - (1.0d+0 - this%ilt_fs0)) / this%ilt_fs0)

end subroutine ILTDefaultIllitization

! ************************************************************************** !

subroutine ILTGeneralIllitization(this,fs,temperature,dt,fi,scale,option)

  use Option_module

  implicit none

  class(ILT_general_type) :: this
  PetscReal, intent(inout) :: fs
  PetscReal, intent(in) :: temperature
  PetscReal, intent(in) :: dt
  PetscReal, intent(out) :: fi
  PetscReal, intent(out) :: scale
  type(option_type), intent(inout) :: option

  PetscReal :: T    ! temperature in Kelvin
  PetscReal :: rate ! temperature-dependent illitization rate in sec^-1

  ! Model based on J. Cuadros & J. Linares, "Experimental Kinetic Study of the
  !   Smectite-to-Illite Transformation," Geochimica et Cosmochimica Acta, 
  !   vol. 60, no. 3, pp. 439-453, 1996

  ! Use Kelvin to calculate rate
  T = temperature + 273.15d0

  ! Check if temperature is above threshold for illitization
  if (temperature >= this%ilt_threshold) then
    ! Negative of illitization rate [L/mol-s]
    rate = this%ilt_freq * &
      exp(-1.0d0 * this%ilt_ea / (IDEAL_GAS_CONSTANT * T))
  else
    rate = 0.0d0
  endif

  ! Fraction smectite - pivot solution based on choice of exponent
  if (this%ilt_exp == 1.0d0) then
    ! n = 1
    fs = fs * exp(-1.0d0 * rate * (this%ilt_K_conc**this%ilt_K_exp) * dt)
  else
    ! n != 1
    fs = (rate * (this%ilt_K_conc**this%ilt_K_exp) * &
         (this%ilt_exp - 1.0d0) * dt + &
         fs**(1.0d0 - this%ilt_exp))**(1.0d0/(1.0d0 - this%ilt_exp))
  endif

  if (fs > 1.0d0) then
    fs = 1.0d0
  elseif (fs < 0.0d0) then
    fs = 0.0d0
  endif

  ! Fraction illite
  fi = 1.0d0 - fs
  
  ! Calculate scale factor
  scale = ((fi - (1.0d+0 - this%ilt_fs0)) / this%ilt_fs0)

end subroutine ILTGeneralIllitization

! ************************************************************************** !

subroutine ILTBaseShiftSorption(this,kd0,ele,material_auxvar,option)

  use Option_module
  use Material_Aux_class

  implicit none

  class(illitization_base_type) :: this
  PetscReal, intent(inout) :: kd0
  character(len=MAXWORDLENGTH), intent(in) :: ele
  class(material_auxvar_type), intent(in) :: material_auxvar
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'Illitization function must be extended to modify ' &
                   //'the kd values of elements in UFD Decay.'
  call PrintErrMsgByRank(option)
  
end subroutine ILTBaseShiftSorption

! ************************************************************************** !

subroutine ILTShiftSorption(this,kd0,ele,material_auxvar,option)

  use Option_module
  use Material_Aux_class

  implicit none

  class(ILT_default_type) :: this
  PetscReal, intent(inout) :: kd0
  character(len=MAXWORDLENGTH), intent(in) :: ele
  class(material_auxvar_type), intent(in) :: material_auxvar
  type(option_type), intent(inout) :: option
  
  class(ilt_kd_effects_type), pointer :: kdl
  character(len=MAXWORDLENGTH) :: fkdele
  character(len=MAXWORDLENGTH) :: fkdmode
  PetscReal, allocatable :: fkd(:)
  PetscInt :: i, j, k
  PetscReal :: scale
  
  if (.not. associated(this%ilt_shift_kd_list)) return
  if (.not. associated(material_auxvar%iltf)) return

  ! Find element and functional properties
  kdl => this%ilt_shift_kd_list
  do i = 1, kdl%num_elements
    fkdele = kdl%f_kd_element(i)
    ! If elements match, proceed
    if (trim(fkdele) == trim(ele)) then
      ! Identify function
      fkdmode = kdl%f_kd_mode(i)
      ! Allocate vector of function values
      select case(fkdmode)
        case ('DEFAULT')
          j = 1
        case default
          option%io_buffer = 'Sorption modification function "' &
                           // trim(fkdmode) &
                           //'" was not found among the available options.'
        call PrintErrMsgByRank(option)
      end select
      allocate(fkd(j))
      ! Populate local vector of function values
      do k = 1, j
        fkd(k) = kdl%f_kd(i,k)
      enddo
      ! Done
      exit
    endif
  enddo
  
  ! Apply function to modify kd
  select case(fkdmode)
    case ('DEFAULT')
      scale = material_auxvar%iltf%ilt_scale
      
      kd0 = kd0 * (1.0d0 + scale*fkd(1))
    case default
      kd0 = kd0
  end select
  
  if (allocated(fkd)) deallocate(fkd)

end subroutine ILTShiftSorption

! ************************************************************************** !

subroutine ILTBaseCheckElements(this,pm_ufd_elements,num,option)

  use Option_module

  implicit none

  class(illitization_base_type) :: this
  PetscInt, intent(in) :: num
  character(len=MAXWORDLENGTH), intent(in) :: pm_ufd_elements(num)
  type(option_type), intent(inout) :: option
  
  return
  
end subroutine ILTBaseCheckElements

! ************************************************************************** !

subroutine ILTCheckElements(this,pm_ufd_elements,num,option)

  use Option_module

  implicit none

  class(ILT_default_type) :: this
  PetscInt, intent(in) :: num
  character(len=MAXWORDLENGTH), intent(in) :: pm_ufd_elements(num)
  type(option_type), intent(inout) :: option
  
  class(ilt_kd_effects_type), pointer :: kdl
  character(len=MAXWORDLENGTH) :: fkdele1, fkdele2
  PetscInt :: i, j
  PetscBool :: found
  
  if (.not. associated(this%ilt_shift_kd_list)) return
  
  kdl => this%ilt_shift_kd_list
  do i = 1, kdl%num_elements
    ! Element specified in illitization function
    fkdele1 = kdl%f_kd_element(i)
    found = PETSC_FALSE
    do j = 1, num
      ! Element specified in UFD Decay
      fkdele2 = pm_ufd_elements(j)
      ! If elements match, proceed to next in list
      if (trim(fkdele1) == trim(fkdele2)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not. found) then
      option%io_buffer = 'Element "'// trim(fkdele1) &
                       //'" listed for kd modification was not found among ' &
                       //'the elements in UFD Decay.'
      call PrintErrMsgByRank(option)
    endif
  enddo
  
end subroutine ILTCheckElements

! ************************************************************************** !

function MaterialTransformCreate()

  implicit none

  class(material_transform_type), pointer :: MaterialTransformCreate
  class(material_transform_type), pointer :: material_transform

  allocate(material_transform)
  material_transform%name = ''
  material_transform%print_me = PETSC_FALSE
  material_transform%test = PETSC_FALSE
  nullify(material_transform%illitization_function)
  nullify(material_transform%next)

  MaterialTransformCreate => material_transform

end function MaterialTransformCreate

! ************************************************************************** !

subroutine MaterialTransformRead(this,input,option)

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(material_transform_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string

  input%ierr = 0
  error_string = 'MATERIAL_TRANSFORM "'//trim(this%name)//'"'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      !------------------------------------------
      case('ILLITIZATION')
        call IllitizationRead(this,input,option)
      !------------------------------------------
      case default
        call InputKeywordUnrecognized(input,keyword, &
               'MATERIAL_TRANSFORM "'//trim(this%name)//'"',option)
    end select
    
  enddo
  
  call InputPopBlock(input,option)

end subroutine MaterialTransformRead

! ************************************************************************** !

subroutine IllitizationRead(this,input,option)

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(material_transform_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string, verify_string
  class(illitization_base_type), pointer :: illitization_function_ptr

  nullify(illitization_function_ptr)

  input%ierr = 0
  error_string = 'ILLITIZATION'
  
  if (associated(this%illitization_function)) then
    option%io_buffer = 'There may only be one instance of '// &
                       trim(error_string) // &
                       ' in MATERIAL_TRANSFORM "'//trim(this%name)//'".'
    call PrintErrMsg(option)
  endif
  
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
      !------------------------------------------
      case('ILLITIZATION_FUNCTION')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option, &
             'ILLITIZATION_FUNCTION',error_string)
        call StringToUpper(word)
        select case(word)
          !-------------------------------------
          case('DEFAULT','HUANG')
            this%illitization_function => ILTDefaultCreate()
          !-------------------------------------
          case('GENERAL','CUADROS_AND_LINARES')
            this%illitization_function => ILTGeneralCreate()
          !-------------------------------------
          case default
            call InputKeywordUnrecognized(input,word, &
                 'ILLITIZATION_FUNCTION',option)
        end select
        call ILTRead(this%illitization_function,input,option)
      !------------------------------------------
      case('TEST')
        this%test = PETSC_TRUE
      !------------------------------------------
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
    class is(ILT_general_type)
      error_string = trim(error_string) // 'GENERAL'
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
          case default
            call ILTDefaultRead(ilf,input,keyword,error_string,'DEFAULT',option)
        end select
      !------------------------------------------
      class is(ILT_general_type)
        select case(trim(keyword))
          case('K_EXP')
            ! Exponent of potassium cation concentration
            call InputReadDouble(input,option,ilf%ilt_K_exp)
            call InputErrorMsg(input,option,'potassium concentration exponent',&
                               'ILLITIZATION, GENERAL')
          case('SMECTITE_EXP')
            ! Exponent of smectite fraction
            call InputReadDouble(input,option,ilf%ilt_exp)
            call InputErrorMsg(input,option,'smectite exponent', &
                               'ILLITIZATION, GENERAL')
          case default
            call ILTDefaultRead(ilf,input,keyword,error_string,'GENERAL',option)
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
  ! Reads in contents of ILLITIZATION_FUNCTION block for the illitization 
  !   base class
  !
  use Option_module
  use Input_Aux_module
  use String_module

  class(illitization_base_type) :: ilf
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH)   :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=*)  :: kind
  type(option_type) :: option
  
  PetscInt :: i,j,jmax
  PetscInt, parameter :: MAX_KD_SIZE = 100
  character(len=MAXWORDLENGTH) :: word
  class(ilt_kd_effects_type), pointer :: shift_kd_list
  PetscReal :: f_kd(MAX_KD_SIZE,10)
  PetscInt :: f_kd_mode_size(MAX_KD_SIZE)
  character(len=MAXWORDLENGTH) :: f_kd_element(MAX_KD_SIZE)
  character(len=MAXWORDLENGTH) :: f_kd_mode(MAX_KD_SIZE)

  select case(keyword)
    case('SMECTITE_INITIAL')
      ! Initial fraction of smectite in the smectite/illite mixture
      call InputReadDouble(input,option,ilf%ilt_fs0)
      call InputErrorMsg(input,option,'initial smectite fraction', &
                         'ILLITIZATION, '//trim(kind)//'')
    case('THRESHOLD_TEMPERATURE')
      ! Specifies the temperature threshold for activating illitization
      call InputReadDouble(input,option, &
                           ilf%ilt_threshold)
      call InputErrorMsg(input,option,'temperature threshold', &
                         'ILLITIZATION, '//trim(kind)//'')
      call InputReadAndConvertUnits(input,ilf%ilt_threshold,'C', &
                                    'ILLITIZATION, '//trim(kind)// &
                                    ', temperature threshold',option)
    case('SHIFT_PERM')
      ! Factor modifying permeability as a function of the illite fraction
      call InputReadDouble(input,option,ilf%ilt_shift_perm)
      call InputErrorMsg(input,option,'permeability shift factor', &
                         'ILLITIZATION, '//trim(kind)//'')
    case('SHIFT_KD')
      ! Factors modifying selected kd values as a function of the
      !   illite fraction
      shift_kd_list => ILTKdEffectsCreate()
      i = 0
      j = 0
      f_kd_mode_size(:) = 0
      f_kd(:,:) = UNINITIALIZED_DOUBLE
      f_kd_mode(:) = ''
      f_kd_element(:) = ''
      
      call InputPushBlock(input,option)
      
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        i = i + 1
        if (i > MAX_KD_SIZE) then
          write(word,*) i-1
          option%io_buffer = 'f_kd array in ILLITIZATION must be' &
            //'allocated larger than ' // trim(adjustl(word)) &
            //' under SHIFT_KD in' // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        
        ! Element
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'f_kd element symbol', &
                           error_string)
        f_kd_element(i) = word
        
        ! Function type
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'f_kd function type', &
                           error_string)
        f_kd_mode(i) = word
        
        ! Function parameters
        select case(f_kd_mode(i))
          case('DEFAULT')
            j = 1
            f_kd_mode_size(i) = j
            call InputReadDouble(input,option,f_kd(i,1))
            call InputErrorMsg(input,option,'f_kd, LINEAR',error_string)
            
            ! Check user values
            if (f_kd(i,1) < -1.0d+0) then
              option%io_buffer = 'Functional parameter #1 in "' &
                               // trim(f_kd_mode(i)) // '" for element "'&
                               // trim(f_kd_element(i)) &
                               //'" must not be less than -1 in ' &
                               //'in ILLITIZATION, '//trim(kind)//'.'
              call PrintErrMsg(option)
              
            endif
          case default
            option%io_buffer = 'Sorption modification function "' &
                             // trim(f_kd_mode(i)) // '" for element "'&
                             // trim(f_kd_element(i)) &
                             //'" was not found among the available options ' &
                             //'in ILLITIZATION, '//trim(kind)//'.'
            call PrintErrMsg(option)
        end select
      enddo
      
      call InputPopBlock(input,option)
      
      if (i == 0) then
        option%io_buffer = 'No f_kd/element combinations specified &
          &under SHIFT_KD in ' // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif
      
      jmax = maxval(f_kd_mode_size)
      allocate(shift_kd_list%f_kd(i,jmax))
      shift_kd_list%f_kd = f_kd(1:i,1:jmax)
      allocate(shift_kd_list%f_kd_element(i))
      shift_kd_list%f_kd_element = f_kd_element(1:i)
      allocate(shift_kd_list%f_kd_mode(i))
      shift_kd_list%f_kd_mode = f_kd_mode(1:i)
      shift_kd_list%num_elements = i
      
      ilf%ilt_shift_kd_list => shift_kd_list
      
      nullify(shift_kd_list)
    case default
      call InputKeywordUnrecognized(input,keyword, &
           'illitization function ('//trim(kind)//')',option)
  end select

end subroutine ILTBaseRead

! ************************************************************************** !

subroutine ILTDefaultRead(ilf,input,keyword,error_string,kind,option)
  !
  ! Reads in contents of ILLITIZATION_FUNCTION block for illitization
  !   default class
  !
  use Option_module
  use Input_Aux_module
  use String_module

  class(ILT_default_type) :: ilf
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH)   :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=*)  :: kind
  type(option_type) :: option
  
  select case(keyword)
    case('EA')
      ! Activation energy in Arrhenius term
      call InputReadDouble(input,option,ilf%ilt_ea)
      call InputErrorMsg(input,option,'activation energy', &
                         'ILLITIZATION, '//trim(kind)//'')
      call InputReadAndConvertUnits(input,ilf%ilt_ea, &
                                    'J/mol','ILLITIZATION, '//trim(kind)// &
                                    ', activation energy',option)
    case('FREQ')
      ! Frequency factor (scaling constant of Arrhenius term)
      call InputReadDouble(input,option,ilf%ilt_freq)
      call InputErrorMsg(input,option,'frequency term', &
                         'ILLITIZATION, '//trim(kind)//'')
      call InputReadAndConvertUnits(input,ilf%ilt_freq, &
                                    'L/s-mol','ILLITIZATION, '//trim(kind)// &
                                    ', frequency term',option)
    case('K_CONC')
      ! Concentration of potassium cation
      call InputReadDouble(input,option,ilf%ilt_K_conc)
      call InputErrorMsg(input,option,'potassium concentration', &
                         'ILLITIZATION, '//trim(kind)//'')
      call InputReadAndConvertUnits(input,ilf%ilt_K_conc,'M',&
                                    'ILLITIZATION, ' //trim(kind)// &
                                    ', potassium concentration',option)
    case default
      call ILTBaseRead(ilf,input,keyword,error_string,kind,option)
  end select

end subroutine ILTDefaultRead

! ************************************************************************** !

subroutine MaterialTransformAddToList(new_mtf,list)

  implicit none

  class(material_transform_type), pointer :: new_mtf
  class(material_transform_type), pointer :: list

  class(material_transform_type), pointer :: cur_mtf

  if (associated(list)) then
    cur_mtf => list
    ! loop to end of list
    do
      if (.not.associated(cur_mtf%next)) exit
      cur_mtf => cur_mtf%next
    enddo
    cur_mtf%next => new_mtf
  else
    list => new_mtf
  endif

end subroutine MaterialTransformAddToList

! ************************************************************************** !

subroutine MaterialTransformConvertListToArray(list,array,option)

  use String_module
  use Option_module

  implicit none

  class(material_transform_type), pointer :: list
  type(material_transform_ptr_type), pointer :: array(:)
  type(option_type) :: option

  class(material_transform_type), pointer :: cur_mtf
  PetscInt :: count

  count = 0
  cur_mtf => list
  do
    if (.not.associated(cur_mtf)) exit
    count = count + 1
    cur_mtf => cur_mtf%next
  enddo

  if (associated(array)) deallocate(array)
  allocate(array(count))

  count = 0
  cur_mtf => list
  do
    if (.not.associated(cur_mtf)) exit
    count = count + 1
    array(count)%ptr => cur_mtf
    if (cur_mtf%test) then
      call OptionSetBlocking(option,PETSC_FALSE)
      if (OptionIsIORank(option)) then
        if (associated(cur_mtf%illitization_function)) then
          call cur_mtf%illitization_function%Test( &
               cur_mtf%name,option)
        endif
      endif
      call OptionSetBlocking(option,PETSC_TRUE)
      call OptionCheckNonBlockingError(option)
    endif
    cur_mtf => cur_mtf%next
  enddo

end subroutine MaterialTransformConvertListToArray

! ************************************************************************** !

function MaterialTransformGetID(material_transform_array, &
           material_transform_name, material_property_name, option)

  use Option_module
  use String_module

  type(material_transform_ptr_type), pointer :: material_transform_array(:)
  character(len=MAXWORDLENGTH) :: material_transform_name
  character(len=MAXWORDLENGTH) :: test1, test2
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: iid, MaterialTransformGetID
  PetscInt :: i, j

  do i = 1, size(material_transform_array)
      test1 = material_transform_array(i)%ptr%name
      do j = 1, size(material_transform_array)
        if (i == j) cycle
        test2 = material_transform_array(j)%ptr%name
        if (test1 == test2) then
          option%io_buffer = 'Duplicate material transform function '//&
                             trim(test2)//&
                             ' has been detected.'
          call PrintErrMsg(option)
        endif
      enddo
  enddo

  MaterialTransformGetID = 0
  do iid = 1, size(material_transform_array)
    if (StringCompare(material_transform_name, &
                      material_transform_array(iid)%ptr%name)) then
      MaterialTransformGetID = iid
      return
    endif
  enddo

  ! MaterialTransformGetID = UNINITIALIZED_INTEGER
  option%io_buffer = 'Material transform function "' // &
                     trim(material_transform_name) // &
                     '" specified in material property "' // &
                     trim(material_property_name) // &
                     '" not found among available functions.'
  call PrintErrMsg(option)

end function MaterialTransformGetID

! ************************************************************************** !

function MaterialTransformCheckILT(material_transform_array,id)
  
  type(material_transform_ptr_type), pointer :: material_transform_array(:)
  PetscInt, intent(in) :: id

  PetscBool :: MaterialTransformCheckILT
  type(illitization_base_type), pointer :: ilt
  
  MaterialTransformCheckILT = PETSC_FALSE
  
  if (associated(material_transform_array(id)%ptr%illitization_function)) then
    select type(ilt => material_transform_array(id)%ptr%illitization_function)
      ! Type must be extended
      class is(ILT_default_type)
        MaterialTransformCheckILT = PETSC_TRUE
    end select
  endif

end function MaterialTransformCheckILT

! ************************************************************************** !

subroutine MaterialTransformInputRecord(material_transform_list)

  implicit none

  class(material_transform_type), pointer :: material_transform_list

  class(material_transform_type), pointer :: cur_mtf
  character(len=MAXWORDLENGTH) :: word1
  PetscInt :: id = INPUT_RECORD_UNIT
  
  class(ilt_kd_effects_type), pointer :: kdl
  PetscInt :: i, j, k
  PetscBool :: inactive
  
  inactive = PETSC_TRUE
  
  cur_mtf => material_transform_list
  do
    if (.not.associated(cur_mtf)) exit
    if (associated(cur_mtf%illitization_function)) then
      select type (ilf => cur_mtf%illitization_function)
        class is (ILT_default_type)
          inactive = PETSC_FALSE
          exit
        end select
    endif
    cur_mtf => cur_mtf%next
  enddo

  if (inactive) return

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
       &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'MATERIAL TRANSFORM FUNCTIONS'

  cur_mtf => material_transform_list
  do
    if (.not.associated(cur_mtf)) exit
    
    if (associated(cur_mtf%illitization_function)) then
      select type (ilf => cur_mtf%illitization_function)
        type is (illitization_base_type)
          exit
        end select
    endif

    write(id,'(a29)',advance='no') 'material transform name: '
    write(id,'(a)') adjustl(trim(cur_mtf%name))
    
    if (associated(cur_mtf%illitization_function)) then
      write(id,'(a29)') '--------------: '
      write(id,'(a29)',advance='no') 'illitization model: '
      select type (ilf => cur_mtf%illitization_function)
        !---------------------------------
        class is (ILT_default_type)
          write(id,'(a)') 'Huang et al., 1993'
          write(id,'(a29)',advance='no') 'initial smectite: '
          write(word1,'(es12.5)') ilf%ilt_fs0
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'frequency: '
          write(word1,'(es12.5)') ilf%ilt_freq
          write(id,'(a)') adjustl(trim(word1))//' L/mol-s'
          write(id,'(a29)',advance='no') 'activation energy: '
          write(word1,'(es12.5)') ilf%ilt_ea
          write(id,'(a)') adjustl(trim(word1))//' J/mol'
          write(id,'(a29)',advance='no') 'K+ concentration: '
          write(word1,'(es12.5)') ilf%ilt_K_conc
          write(id,'(a)') adjustl(trim(word1))//' M'
          write(id,'(a29)',advance='no') 'temperature threshold: '
          write(word1,'(es12.5)') ilf%ilt_threshold
          write(id,'(a)') adjustl(trim(word1))//' C'
          write(id,'(a29)',advance='no') 'shift (permeability): '
          write(word1,'(es12.5)') ilf%ilt_shift_perm
          write(id,'(a)') adjustl(trim(word1))
          if (associated(ilf%ilt_shift_kd_list)) then
            write(id,'(a29)') 'shift (kd): '
            kdl => ilf%ilt_shift_kd_list
            do i = 1, kdl%num_elements
              write(id,'(a29)',advance='no') " "
              write(word1,'(a)') kdl%f_kd_element(i)
              write(id,'(a)',advance='no') adjustl(trim(word1))//" "
              write(word1,'(a)') kdl%f_kd_mode(i)
              write(id,'(a)',advance='no') adjustl(trim(word1))//" "
              select case(kdl%f_kd_mode(i))
                case ('DEFAULT')
                  j = 1
              end select
              do k = 1, j
                write(word1,'(es12.5)') kdl%f_kd(i,k)
                write(id,'(a)',advance='no') adjustl(trim(word1))//" "
              enddo
              write(id,'(a)')
            enddo
          endif
        !---------------------------------
        class is (ILT_general_type)
          write(id,'(a)') 'Cuadros and Linares, 1996'
          write(id,'(a29)',advance='no') 'initial smectite: '
          write(word1,'(es12.5)') ilf%ilt_fs0
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'smectite exponent: '
          write(word1,'(es12.5)') ilf%ilt_exp
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'frequency: '
          write(word1,'(es12.5)') ilf%ilt_freq
          write(id,'(a)') adjustl(trim(word1))//' '
          write(id,'(a29)',advance='no') 'activation energy: '
          write(word1,'(es12.5)') ilf%ilt_ea
          write(id,'(a)') adjustl(trim(word1))//' J/mol'
          write(id,'(a29)',advance='no') 'K+ concentration: '
          write(word1,'(es12.5)') ilf%ilt_K_conc
          write(id,'(a)') adjustl(trim(word1))//' M'
          write(id,'(a29)',advance='no') 'K+ conc. exponent: '
          write(word1,'(es12.5)') ilf%ilt_K_exp
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'temperature threshold: '
          write(word1,'(es12.5)') ilf%ilt_threshold
          write(id,'(a)') adjustl(trim(word1))//' C'
          write(id,'(a29)',advance='no') 'shift (permeability): '
          write(word1,'(es12.5)') ilf%ilt_shift_perm
          write(id,'(a)') adjustl(trim(word1))
          if (associated(ilf%ilt_shift_kd_list)) then
            write(id,'(a29)') 'shift (kd): '
            kdl => ilf%ilt_shift_kd_list
            do i = 1, kdl%num_elements
              write(id,'(a29)',advance='no') " "
              write(word1,'(a)') kdl%f_kd_element(i)
              write(id,'(a)',advance='no') adjustl(trim(word1))//" "
              write(word1,'(a)') kdl%f_kd_mode(i)
              write(id,'(a)',advance='no') adjustl(trim(word1))//" "
              select case(kdl%f_kd_mode(i))
                case ('DEFAULT')
                  j = 1
              end select
              do k = 1, j
                write(word1,'(es12.5)') kdl%f_kd(i,k)
                write(id,'(a)',advance='no') adjustl(trim(word1))//" "
              enddo
              write(id,'(a)')
            enddo
          endif
      end select
    endif

    write(id,'(a29)') '---------------------------: '
    cur_mtf => cur_mtf%next
  enddo

end subroutine MaterialTransformInputRecord

! ************************************************************************** !

recursive subroutine MaterialTransformDestroy(mtf)

  implicit none

  class(material_transform_type), pointer :: mtf

  if (.not. associated(mtf)) return

  call MaterialTransformDestroy(mtf%next)

  if (associated(mtf%illitization_function)) then
    call ILTDestroy(mtf%illitization_function)
  endif

  deallocate(mtf)
  nullify(mtf)

end subroutine MaterialTransformDestroy

end module Material_Transform_module
