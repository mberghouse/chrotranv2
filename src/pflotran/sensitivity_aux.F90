module Sensitivity_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  
  implicit none

  private
  
  !Output format
  PetscInt, parameter, public :: SENSITIVITY_OUTPUT_HDF5          = 0
  PetscInt, parameter, public :: SENSITIVITY_OUTPUT_MATLAB        = 1
  PetscInt, parameter, public :: SENSITIVITY_OUTPUT_BINARY        = 2
  PetscInt, parameter, public :: SENSITIVITY_OUTPUT_ASCII         = 3
  
  !output variable
  PetscInt, parameter, public :: SENSITIVITY_PRESSURE             = 0
  PetscInt, parameter, public :: SENSITIVITY_PERMEABILITY         = 1
  PetscInt, parameter, public :: SENSITIVITY_POROSITY             = 2
  
  !Output times
  PetscInt, parameter, public :: SENSITIVITY_SYNC_OUTPUT          = 0
  PetscInt, parameter, public :: SENSITIVITY_EVERY_X_TIMESTEP     = 1
  PetscInt, parameter, public :: SENSITIVITY_LAST                 = 2
  
  !Miscenalleous
  PetscInt, parameter, public :: SENSITIVITY_MAX_FACE_PER_CELL_OUTPUT = 60
  
  
  type, public :: sensitivity_output_variable_type
    character(len=MAXWORDLENGTH) :: name   ! string that appears in hdf5 file
    character(len=MAXWORDLENGTH) :: units
    PetscInt :: ivar
    type(sensitivity_output_variable_type), pointer :: next
  end type sensitivity_output_variable_type
  
  type, public :: sensitivity_output_variable_list_type
    type(sensitivity_output_variable_type), pointer :: first
    type(sensitivity_output_variable_type), pointer :: last
    !type(sensitivity_output_variable_type), pointer :: array(:)
    PetscInt :: nvars
  end type sensitivity_output_variable_list_type
  
  type, public :: sensitivity_output_option_type
    PetscInt :: plot_number
    PetscBool :: first_plot_flag
    PetscInt :: format !hdf5, binary, matlab
    PetscInt :: output_time_option !sync, x timestep, last
    PetscInt :: output_every_timestep !when x timestep option
    PetscInt :: timestep_per_hdf5_file !x timestep per output file
    PetscReal :: time !simulation time
    PetscBool :: plot_flag
    type(sensitivity_output_variable_list_type), pointer :: output_variables
  end type sensitivity_output_option_type
  
  public :: SensitivityOutputVariableCreate, &
            SensitivityOutputOptionInit, &
            SensitivityOutputOptionIsTimeToOutput, &
            SensitivityAddOutputVariableToList, &
            SensitivityOutputOptionDestroy

contains

! ************************************************************************** !

function  SensitivityOutputVariableCreate()
  ! 
  ! Create an sensitivity output variable
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/12/2021
  ! 
  
  implicit none
  
  type(sensitivity_output_variable_type), pointer :: &
                                    SensitivityOutputVariableCreate
  type(sensitivity_output_variable_type), pointer :: variable
  
  allocate(variable)
  variable%name = ""
  variable%units = ""
  variable%ivar = UNINITIALIZED_INTEGER
  nullify(variable%next)
  
  SensitivityOutputVariableCreate => variable
  
end function SensitivityOutputVariableCreate

! ************************************************************************** !

subroutine SensitivityOutputOptionInit(sensitivity_output_option)
  ! 
  ! Initialize an output sensitivity option instance
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  ! 
  
  implicit none
  
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  
  sensitivity_output_option%plot_number = 0
  sensitivity_output_option%first_plot_flag = PETSC_TRUE
  sensitivity_output_option%format = SENSITIVITY_OUTPUT_HDF5
  sensitivity_output_option%output_time_option = SENSITIVITY_SYNC_OUTPUT
  sensitivity_output_option%output_every_timestep = 0
  sensitivity_output_option%timestep_per_hdf5_file = 0 !no limit
  sensitivity_output_option%time = 0.d0
  sensitivity_output_option%plot_flag = PETSC_FALSE
  allocate(sensitivity_output_option%output_variables)
  call SensitivityOutputVaraibleListInit( &
                                 sensitivity_output_option%output_variables)
  
end subroutine SensitivityOutputOptionInit

! ************************************************************************** !

subroutine SensitivityOutputVaraibleListInit(list)
  ! 
  ! Initialize an output sensitivity option output variable list
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  ! 
  
  implicit none
  
  type(sensitivity_output_variable_list_type), pointer :: list
  
  nullify(list%first)
  nullify(list%last)
  !nullify(list%array)
  list%nvars = 0 
  
end subroutine SensitivityOutputVaraibleListInit

! ************************************************************************** !

subroutine SensitivityAddOutputVariableToList(sensitivity_output_option, &
                                              new_variable)
  ! 
  ! Add an output sensitivity variable to the output variable list
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  ! 
  
  implicit none
  
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  type(sensitivity_output_variable_type), pointer :: new_variable
  nullify(new_variable%next)
  
  if (.not. associated(sensitivity_output_option%output_variables%first)) then
    sensitivity_output_option%output_variables%first => new_variable
  else
    sensitivity_output_option%output_variables%last%next => new_variable
  endif
  sensitivity_output_option%output_variables%last => new_variable
  
  sensitivity_output_option%output_variables%nvars = &
                           sensitivity_output_option%output_variables%nvars+1
  
end subroutine SensitivityAddOutputVariableToList

! ************************************************************************** !

subroutine SensitivityOutputOptionIsTimeToOutput(sensitivity_output_option,&
                                                 timestep_flag,last_flag,&
                                                 sync_flag)
  ! 
  ! Determine if it's the time to output the sensitivity
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  ! 
  
  implicit none
  
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  PetscBool :: timestep_flag, last_flag, sync_flag
  
  select case(sensitivity_output_option%output_time_option)
    case(SENSITIVITY_SYNC_OUTPUT)
      if (sync_flag) sensitivity_output_option%plot_flag = PETSC_TRUE
    case(SENSITIVITY_EVERY_X_TIMESTEP)
      if (timestep_flag) sensitivity_output_option%plot_flag = PETSC_TRUE
    case(SENSITIVITY_LAST)
      if (last_flag) sensitivity_output_option%plot_flag = PETSC_TRUE
  end select
  
end subroutine SensitivityOutputOptionIsTimeToOutput

! ************************************************************************** !

subroutine SensitivityDestroyOutputVariableList(var_list)
  ! 
  ! Destroy an sensitivity output variable list object
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/12/2021
  ! 
  
  implicit none
  
  type(sensitivity_output_variable_list_type), pointer :: var_list
  type(sensitivity_output_variable_type), pointer :: var, prev_var
  
  if (.not.associated(var_list)) return
  
  var => var_list%first
  do 
    if (.not.associated(var)) exit
    prev_var => var
    var => var%next
    nullify(prev_var%next)
  enddo
  
  var_list%nvars = 0
  nullify(var_list%first)
  nullify(var_list%last)
  
  deallocate(var_list)
  nullify(var_list)

end subroutine SensitivityDestroyOutputVariableList

! ************************************************************************** !

subroutine SensitivityOutputOptionDestroy(sensitivity_output_option)
  ! 
  ! Destroy an output sensitivity option instance
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/12/2021
  ! 
  
  implicit none
  
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  
  call SensitivityDestroyOutputVariableList( &
                                 sensitivity_output_option%output_variables)
  
end subroutine SensitivityOutputOptionDestroy

! ************************************************************************** !
 
end module
