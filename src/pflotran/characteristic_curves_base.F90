module Characteristic_Curves_Base_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private
  
  PetscReal, parameter, public :: DEFAULT_PCMAX = 1.d9
  
  type, public :: polynomial_type
    PetscReal :: low
    PetscReal :: high
    PetscReal :: coefficients(4)
  end type polynomial_type  
  
!-----------------------------------------------------------------------------
!-- Saturation Functions -----------------------------------------------------
!-----------------------------------------------------------------------------
  type, public :: sat_func_base_type
    type(polynomial_type), pointer :: sat_poly
    type(polynomial_type), pointer :: pres_poly
    PetscReal :: Sr
    PetscReal :: pcmax
    PetscBool :: analytical_derivative_available
    PetscBool :: calc_int_tension
    PetscBool :: unit_test
    character(len=MAXWORDLENGTH) :: input_filename
  contains
    procedure, public :: Init => SFBaseInit
    procedure, public :: Verify => SFBaseVerify
    procedure, public :: Test => SFBaseTest
    procedure, public :: UnitTest => SFBaseUnitTest
    procedure, public :: SetupPolynomials => SFBaseSetupPolynomials
    procedure, public :: CapillaryPressure => SFBaseCapillaryPressure
    procedure, public :: Saturation => SFBaseSaturation
    procedure, public :: D2SatDP2 => SFBaseD2SatDP2
    procedure, public :: CalcInterfacialTension => SFBaseSurfaceTension
  end type sat_func_base_type

!-----------------------------------------------------------------------------
!-- Relative Permeability Functions ------------------------------------------
!-----------------------------------------------------------------------------  
  type, public :: rel_perm_func_base_type
    type(polynomial_type), pointer :: poly
    PetscReal :: Sr
    PetscReal :: Srg
    PetscBool :: analytical_derivative_available
    PetscBool :: unit_test
    character(len=MAXWORDLENGTH) :: input_filename
  contains
    procedure, public :: Init => RPFBaseInit
    procedure, public :: Verify => RPFBaseVerify
    procedure, public :: Test => RPF_Base_Test
    procedure, public :: UnitTest => RPF_Base_UnitTest
    procedure, public :: SetupPolynomials => RPFBaseSetupPolynomials
    procedure, public :: RelativePermeability => RPF_Base_RelPerm
  end type rel_perm_func_base_type
  
  public :: PolynomialCreate, &
            PolynomialDestroy, &
            SFBaseInit, &
            SFBaseVerify, &
            SFBaseTest, &
            SFBaseCapillaryPressure, &
            SFBaseSaturation, &
            RPFBaseInit, &
            RPFBaseVerify, &
            RPF_Base_Test, &
            RPF_Base_UnitTest, &
            RPF_Base_RelPerm, &
            SaturationFunctionDestroy, &
            PermeabilityFunctionDestroy

contains

! ************************************************************************** !

function PolynomialCreate()

  implicit none
  
  type(polynomial_type), pointer :: PolynomialCreate  

  allocate(PolynomialCreate)
  PolynomialCreate%low = 0.d0
  PolynomialCreate%high = 0.d0
  PolynomialCreate%coefficients(:) = 0.d0
  
end function PolynomialCreate

! ************************************************************************** !
! ************************************************************************** !

subroutine SFBaseInit(this)

  implicit none
  
  class(sat_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%sat_poly)
  nullify(this%pres_poly)
  this%Sr = UNINITIALIZED_DOUBLE
  this%pcmax = DEFAULT_PCMAX
  this%analytical_derivative_available = PETSC_FALSE
  this%calc_int_tension = PETSC_FALSE
  this%unit_test = PETSC_FALSE
  this%input_filename = ''
  
end subroutine SFBaseInit

! ************************************************************************** !

subroutine SFBaseVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  
  
  if (Uninitialized(this%Sr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                            name)
    call PrintErrMsg(option)
  endif
  
  if ((.not.this%analytical_derivative_available) .and. &
      (.not.option%flow%numerical_derivatives)) then
    option%io_buffer = 'Analytical derivatives are not available for the &
      &capillary pressure - saturation function chosen: ' // &
      trim(name)
    call PrintErrMsg(option)
  endif
  
end subroutine SFBaseVerify

! ************************************************************************** !

subroutine RPFBaseInit(this)

  implicit none
  
  class(rel_perm_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%poly)
  this%Sr = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  this%analytical_derivative_available = PETSC_FALSE
  this%unit_test = PETSC_FALSE
  this%input_filename = ''
  
end subroutine RPFBaseInit

! ************************************************************************** !

subroutine RPFBaseVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rel_perm_func_base_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  if (Uninitialized(this%Sr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                            name)
    call PrintErrMsg(option)
  endif
  
  if ((.not.this%analytical_derivative_available) .and. &
      (.not.option%flow%numerical_derivatives)) then
    option%io_buffer = 'Analytical derivatives are not available for the &
      &relative permeability function chosen: ' // trim(name)
    call PrintErrMsg(option)
  endif
  
end subroutine RPFBaseVerify

! ************************************************************************** !

subroutine SFBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing saturation functions

  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  option%io_buffer = 'SF Smoothing not supported for ' // trim(error_string)
  call PrintErrMsg(option)
  
end subroutine SFBaseSetupPolynomials

! ************************************************************************** !

subroutine RPFBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing relative permeability functions

  use Option_module
  
  implicit none
  
  class(rel_perm_func_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  option%io_buffer = 'RPF Smoothing not supported for ' // trim(error_string)
  call PrintErrMsg(option)
  
end subroutine RPFBaseSetupPolynomials

! ************************************************************************** !

subroutine SFBaseCapillaryPressure(this,liquid_saturation, & 
                                   capillary_pressure,dpc_dsatl,option)
  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseCapillaryPressure must be extended.'
  call PrintErrMsg(option)
  
end subroutine SFBaseCapillaryPressure

! ************************************************************************** !

subroutine SFBaseSaturation(this,capillary_pressure, &
                            liquid_saturation,dsat_dpres,option)
  use Option_module

  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseSaturation must be extended.'
  call PrintErrMsg(option)
  
end subroutine SFBaseSaturation

! ************************************************************************** !

subroutine SFBaseD2SatDP2(this,pc,d2s_dp2,option)

  use Option_module

  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: pc
  PetscReal, intent(out) :: d2s_dp2
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseD2SatDP2 must be extended.'
  call PrintErrMsg(option)
  
end subroutine SFBaseD2SatDP2

! ************************************************************************** !

subroutine SFBaseTest(this,cc_name,option)

  use Option_module
  use Material_Aux_class

  implicit none
  
  class(sat_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: num_values = 101
  PetscReal :: pc, pc_increment
  PetscReal :: capillary_pressure(num_values)
  PetscReal :: liquid_saturation(num_values)
  PetscReal :: dpc_dsatl(num_values)
  PetscReal :: dpc_dsatl_numerical(num_values)
  PetscReal :: dsat_dpres(num_values)
  PetscReal :: dsat_dpres_numerical(num_values)
  PetscReal :: capillary_pressure_pert
  PetscReal :: liquid_saturation_pert
  PetscReal :: perturbation
  PetscReal :: pert
  PetscReal :: dummy_real
  PetscInt :: count, i
 
  ! calculate saturation as a function of capillary pressure
  ! start at 1 Pa up to maximum capillary pressure
  pc = 1.d0
  pc_increment = 1.d0
  perturbation = 1.d-6
  count = 0
  do
    if (pc > this%pcmax) exit
    count = count + 1
    call this%Saturation(pc,liquid_saturation(count),dsat_dpres(count),option)
    capillary_pressure(count) = pc
    ! calculate numerical derivative dsat_dpres_numerical
    capillary_pressure_pert = pc + pc*perturbation
    call this%Saturation(capillary_pressure_pert,liquid_saturation_pert, &
                         dummy_real,option)
    dsat_dpres_numerical(count) = (liquid_saturation_pert - &
         & liquid_saturation(count))/(pc*perturbation)*(-1.d0) ! dPc/dPres
    ! get next value for pc
    if (pc > 0.99d0*pc_increment*10.d0) pc_increment = pc_increment*10.d0
    pc = pc + pc_increment
  enddo

  write(string,*) cc_name
  string = trim(cc_name) // '_sat_from_pc.dat'
  open(unit=86,file=string)
  write(86,*) '"capillary pressure", "saturation", "dsat/dpres", &
              &"dsat/dpres_numerical"'
  do i = 1, count
    write(86,'(4es14.6)') capillary_pressure(i), liquid_saturation(i), &
                          dsat_dpres(i), dsat_dpres_numerical(i)
  enddo
  close(86)

 ! calculate capillary pressure as a function of saturation
  do i = 1, num_values
    liquid_saturation(i) = dble(i-1)*0.01d0
    if (liquid_saturation(i) < 1.d-7) then
      liquid_saturation(i) = 1.d-7
    else if (liquid_saturation(i) > (1.d0-1.d-7)) then
      liquid_saturation(i) = 1.d0-1.d-7
    endif
    call this%CapillaryPressure(liquid_saturation(i), &
                                capillary_pressure(i),dpc_dsatl(i),option)
    ! calculate numerical derivative dpc_dsatl_numerical
    pert = liquid_saturation(i) * perturbation
    if (liquid_saturation(i) > 0.5d0) then
      pert = -1.d0 * pert
    endif
    liquid_saturation_pert = liquid_saturation(i) + pert
    call this%CapillaryPressure(liquid_saturation_pert, &
                                capillary_pressure_pert,dummy_real,option)
    dpc_dsatl_numerical(i) = (capillary_pressure_pert - &
         & capillary_pressure(i))/pert 
  enddo
  count = num_values

  write(string,*) cc_name
  string = trim(cc_name) // '_pc_from_sat.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "capillary pressure", "dpc/dsat", &
              &dpc_dsat_numerical"'
  do i = 1, count
    write(86,'(4es14.6)') liquid_saturation(i), capillary_pressure(i), &
                          dpc_dsatl(i), dpc_dsatl_numerical(i)
  enddo
  close(86)

end subroutine SFBaseTest

! ************************************************************************** !

subroutine SFBaseUnitTest(this,cc_name,input_filename,option)

  use Option_module
  use Material_Aux_class

  implicit none
  
  class(sat_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  character(len=MAXWORDLENGTH) :: input_filename
  type(option_type), intent(inout) :: option

  character(len=MAXWORDLENGTH), pointer :: name(:)      ! [-]
  character(len=MAXWORDLENGTH), pointer :: temp_name(:) ! [-]
  PetscReal, pointer :: capillary_pressure(:)           ! [Pa]
  PetscReal, pointer :: temp_capillary_pressure(:)      ! [Pa]
  PetscReal, pointer :: corr_liq_saturation(:)          ! [-]
  PetscReal, pointer :: temp_corr_liq_saturation(:)     ! [-]
  PetscReal, pointer :: liq_saturation(:)               ! [-]
  PetscReal, pointer :: temp_liq_saturation(:)          ! [-]
  PetscReal, pointer :: corr_capillary_pressure(:)      ! [Pa]
  PetscReal, pointer :: temp_corr_capillary_pressure(:) ! [Pa]
  PetscReal :: liq_saturation_calcd, capillary_press_calcd
  PetscReal, parameter :: tolerance = 1.d-8
  PetscReal :: dum1
  PetscReal :: diff
  PetscInt :: i, j, k
  character(len=MAXWORDLENGTH) :: pass_fail
  character(len=MAXWORDLENGTH) :: filename_out, id
  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscInt :: values(8)
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time

  allocate(temp_name(150))
  allocate(temp_capillary_pressure(150))
  allocate(temp_corr_liq_saturation(150))
  allocate(temp_liq_saturation(150))
  allocate(temp_corr_capillary_pressure(150))

  i = 1
  open(action='read', file=trim(input_filename), iostat=rc_in, &
       newunit=fu_in)
  read(fu_in, *) ! skip header line
  do
    read (fu_in, *, iostat=rc_in) temp_name(i), &
                                  temp_capillary_pressure(i), &
                                  temp_corr_liq_saturation(i), &
                                  temp_liq_saturation(i), &
                                  temp_corr_capillary_pressure(i)
    if (rc_in /= 0) exit 
    i = i + 1 
    if (i > 150) exit
  enddo

  allocate(name(i-1))
  allocate(capillary_pressure(i-1))
  allocate(corr_liq_saturation(i-1))
  allocate(liq_saturation(i-1))
  allocate(corr_capillary_pressure(i-1))
  name(:) = temp_name(1:i-1)
  capillary_pressure(:) = temp_capillary_pressure(1:i-1)
  corr_liq_saturation(:) = temp_corr_liq_saturation(1:i-1)
  liq_saturation(:) = temp_liq_saturation(1:i-1)
  corr_capillary_pressure(:) = temp_corr_capillary_pressure(1:i-1)

  filename_out = 'sf_' // trim(cc_name) //'.out'

  open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
  call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
  write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                  time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
  write(fu_out,*)
  write(fu_out,'(a)') 'NOTE: The input file provided was:'
  write(fu_out,'(a,a)') '      ', trim(input_filename)
  write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                               tolerance, '.'
  write(fu_out,*)

  i = 0
  do k=1,size(name)
    write(fu_out,'(a,I3,a)') '||-----------TEST-#',k,'-----------------------&
                             &--------------||'
    write(fu_out,'(a)') '  [in]  characteristic curves name:'
    write(fu_out,'(a)') name(k)
    write(fu_out,'(a)') '  [in]  capillary pressure [Pa]:'
    write(fu_out,'(d17.10)') capillary_pressure(k)
    call this%Saturation(capillary_pressure(k),liq_saturation_calcd,dum1,option)
    write(fu_out,'(a)') '  [out]  liquid saturation [-]:'
    write(fu_out,'(d17.10)') liq_saturation_calcd
    write(fu_out,'(a)') '  [correct]  liquid saturation [-]:'
    write(fu_out,'(d17.10)') corr_liq_saturation(k)
    diff = abs(corr_liq_saturation(k)-liq_saturation_calcd)
    if (diff > (tolerance*corr_liq_saturation(k))) then
      pass_fail = 'FAIL!'
      i = i + 1
    else
      pass_fail = 'pass'
    endif
    write(fu_out,'(a)') trim(pass_fail)

    write(fu_out,'(a)') '  [in]  liquid saturation [-]:'
    write(fu_out,'(d17.10)') liq_saturation(k)
    call this%CapillaryPressure(liq_saturation(k),capillary_press_calcd, &
                                dum1,option)
    write(fu_out,'(a)') '  [out]  capillary pressure [Pa]:'
    write(fu_out,'(d17.10)') capillary_press_calcd
    write(fu_out,'(a)') '  [correct]  capillary pressure [Pa]:'
    write(fu_out,'(d17.10)') corr_capillary_pressure(k)
    diff = abs(corr_capillary_pressure(k)-capillary_press_calcd)
    if (diff > (tolerance*corr_capillary_pressure(k))) then
      pass_fail = 'FAIL!'
      i = i + 1
    else
      pass_fail = 'pass'
    endif
    write(fu_out,'(a)') trim(pass_fail)

    write(fu_out,*)
  enddo

  write(fu_out,'(a)') 'TEST SUMMARY:'
  if (i == 0) then
    write(fu_out,'(a)') ' All tests passed!'
  else
    write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
  endif

  close(fu_out)
  close(fu_in)

  deallocate(temp_name)
  deallocate(temp_capillary_pressure)
  deallocate(temp_corr_liq_saturation)
  deallocate(temp_liq_saturation)
  deallocate(temp_corr_capillary_pressure)
  deallocate(name)
  deallocate(capillary_pressure)
  deallocate(corr_liq_saturation)
  deallocate(liq_saturation)
  deallocate(corr_capillary_pressure)

end subroutine SFBaseUnitTest

! ************************************************************************** !
! ************************************************************************** !

subroutine RPF_Base_RelPerm(this,liquid_saturation,relative_permeability, &
                            dkr_sat,option)
  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'RPF_Base_RelPerm must be extended.'
  call PrintErrMsg(option)
  
end subroutine RPF_Base_RelPerm

! ************************************************************************** !

subroutine RPF_Base_Test(this,cc_name,phase,option)

  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  character(len=MAXWORDLENGTH) :: phase
  type(option_type), intent(inout) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, parameter :: num_values = 101
  PetscReal :: perturbation
  PetscReal :: liquid_saturation(num_values)
  PetscReal :: liquid_saturation_pert(num_values)
  PetscReal :: kr(num_values)
  PetscReal :: kr_pert(num_values)
  PetscReal :: dkr_dsat(num_values)
  PetscReal :: dkr_dsat_numerical(num_values)
  PetscReal :: dummy_real(num_values)
  
  perturbation = 1.d-6

  do i = 1, num_values
    liquid_saturation(i) = dble(i-1)*0.01d0
    call this%RelativePermeability(liquid_saturation(i),kr(i),dkr_dsat(i), &
                                   option)
    ! calculate numerical derivative dkr_dsat_numerical
    liquid_saturation_pert(i) = liquid_saturation(i) &
                                + liquid_saturation(i)*perturbation
    call this%RelativePermeability(liquid_saturation_pert(i),kr_pert(i), &
                                   dummy_real(i),option)
    if( i>1 ) then
      dkr_dsat_numerical(i) = (kr_pert(i) - kr(i))/ &
                              (liquid_saturation(i)*perturbation)
    else
! Trap case of i=0 as liquid_saturation is 0 and will otherwise divide by zero
      dkr_dsat_numerical(i) = 0.0
    endif
  enddo

  write(string,*) cc_name
  string = trim(cc_name) // '_' // trim(phase) // '_rel_perm.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "' // trim(phase) // ' relative permeability", "' &
              // trim(phase) // ' dkr/dsat", "' // trim(phase) // &
              ' dkr/dsat_numerical"'
  do i = 1, size(liquid_saturation)
    write(86,'(4es14.6)') liquid_saturation(i), kr(i), dkr_dsat(i), &
                          dkr_dsat_numerical(i)
  enddo
  close(86)

end subroutine RPF_Base_Test

! ************************************************************************** !

subroutine RPF_Base_UnitTest(this,cc_name,phase,input_filename,option)

  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  character(len=MAXWORDLENGTH) :: phase
  character(len=MAXWORDLENGTH) :: input_filename
  type(option_type), intent(inout) :: option

  character(len=MAXWORDLENGTH), pointer :: name(:)      ! [-]
  character(len=MAXWORDLENGTH), pointer :: temp_name(:) ! [-]
  PetscReal, pointer :: saturation(:)                   ! [-]
  PetscReal, pointer :: temp_saturation(:)              ! [-]
  PetscReal, pointer :: rel_perm(:)                     ! [m2]
  PetscReal, pointer :: corr_rel_perm(:)                ! [m2]
  PetscReal, pointer :: temp_corr_rel_perm(:)           ! [m2]          
  PetscReal, parameter :: tolerance = 1.d-8
  PetscReal :: dum1
  PetscReal :: diff
  PetscInt :: i, j, k
  character(len=MAXWORDLENGTH) :: pass_fail
  character(len=MAXWORDLENGTH) :: filename_out, id
  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscBool :: input_file_given
  PetscInt :: values(8)
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time

  allocate(temp_name(99))
  allocate(temp_saturation(99))
  allocate(temp_corr_rel_perm(99))

  i = 1
  open(action='read', file=trim(input_filename), iostat=rc_in, &
       newunit=fu_in)
  read(fu_in, *) ! skip header line
  do
    read (fu_in, *, iostat=rc_in) temp_name(i), &
                                  temp_saturation(i), &
                                  temp_corr_rel_perm(i)
    if (rc_in /= 0) exit 
    i = i + 1 
    if (i > 99) exit
  enddo

  allocate(name(i-1))
  allocate(saturation(i-1))
  allocate(corr_rel_perm(i-1))
  allocate(rel_perm(i-1))
  name(:) = temp_name(1:i-1)
  saturation(:) = temp_saturation(1:i-1)
  corr_rel_perm(:) = temp_corr_rel_perm(1:i-1)

  filename_out = 'rpf_' // trim(phase) // '_' // trim(cc_name) // '.out'

  open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
  call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
  write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                  time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
  write(fu_out,*)
  write(fu_out,'(a)') 'NOTE: The input file provided was:'
  write(fu_out,'(a,a)') '      ', trim(input_filename)
  write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                               tolerance, '.'
  write(fu_out,*)

  i = 0
  do k=1,size(name)
    write(fu_out,'(a,I2,a)') '||-----------TEST-#',k,'-----------------------&
                             &--------------||'
    write(fu_out,'(a)') '  [in]  characteristic curves name:'
    write(fu_out,'(a)') name(k)
    write(fu_out,'(a,a,a)') '  [in]  saturation (', trim(phase), ') [Pa]:'
    write(fu_out,'(d17.10)') saturation(k)

    call this%RelativePermeability(saturation(k),rel_perm(k),dum1,option)
    write(fu_out,'(a,a,a)') '  [out]  relative permeability (', trim(phase), &
                            ') [m2]:'
    write(fu_out,'(d17.10)') rel_perm(k)
    write(fu_out,'(a,a,a)') '  [correct]  relative permeability (', &
                            trim(phase), ') [m2]:'
    write(fu_out,'(d17.10)') corr_rel_perm(k)
    diff = abs(corr_rel_perm(k)-rel_perm(k))
    if (diff > (tolerance*corr_rel_perm(k))) then
      pass_fail = 'FAIL!'
      i = i + 1
    else
      pass_fail = 'pass'
    endif
    write(fu_out,'(a)') trim(pass_fail)

    write(fu_out,*)
  enddo
  write(fu_out,'(a)') 'TEST SUMMARY:'
  if (i == 0) then
    write(fu_out,'(a)') ' All tests passed!'
  else
    write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
  endif

  close(fu_out)
  close(fu_in)

  deallocate(name)
  deallocate(saturation)
  deallocate(corr_rel_perm)
  deallocate(rel_perm)
  deallocate(temp_name)
  deallocate(temp_saturation)
  deallocate(temp_corr_rel_perm)

end subroutine RPF_Base_UnitTest

! ************************************************************************** !

subroutine PolynomialDestroy(poly)
  !
  ! Destroys a polynomial smoother
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  !

  implicit none

  type(polynomial_type), pointer :: poly

  if (.not.associated(poly)) return

  deallocate(poly)
  nullify(poly)

end subroutine PolynomialDestroy
  
! ************************************************************************** !

subroutine SaturationFunctionDestroy(sf)
  !
  ! Destroys a saturuation function
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  !

  implicit none

  class(sat_func_base_type), pointer :: sf

  if (.not.associated(sf)) return
  
  call PolynomialDestroy(sf%sat_poly)
  call PolynomialDestroy(sf%sat_poly)
  deallocate(sf)
  nullify(sf)

end subroutine SaturationFunctionDestroy

! ************************************************************************** !

subroutine PermeabilityFunctionDestroy(rpf)
  ! 
  ! Destroys a saturuation function
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(rel_perm_func_base_type), pointer :: rpf
  
  if (.not.associated(rpf)) return
  
  call PolynomialDestroy(rpf%poly)
  deallocate(rpf)
  nullify(rpf)

end subroutine PermeabilityFunctionDestroy

subroutine SFBaseSurfaceTension(this,T,sigma)
  
  !Surface tension of water equation from Revised Release on Surface
  !Tension of Ordinary Water Substance, June 2014. Valid from -25C to
  !373 C
  
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: sigma
  
  PetscReal, parameter :: Tc = 647.096d0
  PetscReal, parameter :: B = 235.8d0
  PetscReal, parameter :: b_2 = -0.625d0
  PetscReal, parameter :: mu = 1.256d0
  PetscReal, parameter :: sigma_base = 0.073d0
  PetscReal :: Temp
  PetscReal :: tao
  
  Temp=T+273.15d0
  
  if (T <= 373.d0) then
    tao = 1.d0-Temp/Tc
    sigma = B*(tao**mu)*(1+b_2*tao)
    sigma = sigma * 1.d-3
  else
    sigma = 0.d0
  endif
  sigma= sigma/sigma_base

  !TOUGH3 way (not pressure-dependent)
  !if (Temp >= 101) sigma = 0
  
end subroutine SFBaseSurfaceTension
  
end module Characteristic_Curves_Base_module
