
module Output_Sensitivity_module

#include "petsc/finclude/petscsys.h"
  use petscsys
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Logging_module
  use hdf5
  use PFLOTRAN_Constants_module
  
  implicit none

  private
  
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscInt, parameter :: OUTPUT_HDF5=0, OUTPUT_MATLAB=1, OUTPUT_BINARY=2
  PetscInt, parameter :: SYNC_OUTPUT=0, EVERY_X_TIMESTEP=1, LAST=2
  
  !output variable
  PetscInt, parameter :: PRESSURE         = 0
  PetscInt, parameter :: PERMEABILITY     = 1
  PetscInt, parameter :: POROSITY         = 2
  
  type, public :: output_sensitivity_variable_list_type
    type(output_sensitivity_variable_type), pointer :: first
    type(output_sensitivity_variable_type), pointer :: last
    !type(output_sensitivity_variable_type), pointer :: array(:)
    PetscInt :: nvars
  end type output_sensitivity_variable_list_type
  
  type, public :: output_sensitivity_variable_type
    character(len=MAXWORDLENGTH) :: name   ! string that appears in hdf5 file
    character(len=MAXWORDLENGTH) :: units
    PetscInt :: ivar
    type(output_sensitivity_variable_type), pointer :: next
  end type output_sensitivity_variable_type
  
  type, public :: output_sensitivity_option_type
    PetscInt :: plot_number
    PetscInt :: format !hdf5, binary, matlab
    PetscInt :: output_time_option !sync, x timestep, last
    PetscInt :: output_every_timestep !when x timestep option
    PetscInt :: timestep_per_hdf5_file !x timestep per output file
    PetscReal :: time !simulation time
    PetscBool :: plot_flag
    type(output_sensitivity_variable_list_type), pointer :: output_variables
  end type output_sensitivity_option_type
  
  public :: OutputSensitivity, &
            OutputSensitivityOptionInit, &
            OutputSensitivityOptionIsTimeToOutput, &
            OutputSensitivityOptionSetOutputTimeOption, &
            OutputSensitivityOptionSetOutputFormat, &
            OutputSensitivityAddVariable

contains

! ************************************************************************** !

subroutine OutputSensitivityOptionInit(output_sensitivity_option)
  ! 
  ! Initialize an output sensitivity option instance
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  ! 
  
  implicit none
  
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  
  output_sensitivity_option%plot_number = 0
  output_sensitivity_option%format = OUTPUT_HDF5
  output_sensitivity_option%output_time_option = SYNC_OUTPUT
  output_sensitivity_option%output_every_timestep = 0
  output_sensitivity_option%timestep_per_hdf5_file = 0 !no limit
  output_sensitivity_option%time = 0.d0
  output_sensitivity_option%plot_flag = PETSC_FALSE
  allocate(output_sensitivity_option%output_variables)
  call OutputSensitivityVaraibleListInit( &
                                 output_sensitivity_option%output_variables)
  
end subroutine OutputSensitivityOptionInit

! ************************************************************************** !

subroutine OutputSensitivityVaraibleListInit(list)
  ! 
  ! Initialize an output sensitivity option output variable list
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  ! 
  
  implicit none
  
  type(output_sensitivity_variable_list_type), pointer :: list
  
  nullify(list%first)
  nullify(list%last)
  !nullify(list%array)
  list%nvars = 0 
  
end subroutine OutputSensitivityVaraibleListInit

! ************************************************************************** !

subroutine OutputSensitivityAddVariable(output_sensitivity_option,word)
  ! 
  ! Add an output sensitivity variable to the output variable list
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  ! 
  
  implicit none
  
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  character(len=MAXWORDLENGTH) :: word
  
  type(output_sensitivity_variable_type), pointer :: variable
  allocate(variable)
  
  select case(trim(word))
    case('PRESSURE')
      variable%name = 'Pressure'
      variable%units = '1/Pa' ! TODO (moise)
      variable%ivar = PRESSURE
    case('PERMEABILITY')
      variable%name = 'Permeability'
      variable%units = '1/m/s' ! TODO (moise)
      variable%ivar = PERMEABILITY
    case('POROSITY')
      variable%name = 'Porosity'
      variable%units = '' ! TODO (moise)
      variable%ivar = POROSITY
  end select
  
  if (.not. associated(output_sensitivity_option%output_variables%first)) then
    output_sensitivity_option%output_variables%first => variable
  else
    output_sensitivity_option%output_variables%last%next => variable
  endif
  output_sensitivity_option%output_variables%last => variable
  
  output_sensitivity_option%output_variables%nvars = &
                           output_sensitivity_option%output_variables%nvars+1
  
end subroutine OutputSensitivityAddVariable

! ************************************************************************** !

subroutine OutputSensitivityOptionIsTimeToOutput(output_sensitivity_option,&
                                                 timestep_flag,last_flag,&
                                                 sync_flag)
  ! 
  ! Determine if it's the time to output the sensitivity
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  ! 
  
  implicit none
  
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  PetscBool :: timestep_flag, last_flag, sync_flag
  
  select case(output_sensitivity_option%output_time_option)
    case(SYNC_OUTPUT)
      if (sync_flag) output_sensitivity_option%plot_flag = PETSC_TRUE
    case(EVERY_X_TIMESTEP)
      if (timestep_flag) output_sensitivity_option%plot_flag = PETSC_TRUE
    case(LAST)
      if (last_flag) output_sensitivity_option%plot_flag = PETSC_TRUE
  end select
  
end subroutine OutputSensitivityOptionIsTimeToOutput

! ************************************************************************** !

subroutine OutputSensitivityOptionSetOutputTimeOption( &
                                                output_sensitivity_option,&
                                                option,word,input)
  ! 
  ! Set the time when to output the sensitivity
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  ! 
  
  use Input_Aux_module
  use Option_module
  
  implicit none
  
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  type(input_type), pointer :: input
  
  select case(trim(word))
    case('SYNC_WITH_SNAPSHOT_FILE')
      output_sensitivity_option%output_time_option = SYNC_OUTPUT
    case('PERIODIC_TIMESTEP')
      output_sensitivity_option%output_time_option = EVERY_X_TIMESTEP
      call InputReadInt(input,option,&
                        output_sensitivity_option%output_every_timestep)
    case('LAST_TIMESTEP')
      output_sensitivity_option%output_time_option = LAST
    case default
      call InputKeywordUnrecognized(input,word,'SENSITIVITY_RICHARDS,&
                                    &OUTPUT_TIME',option)
  end select

end subroutine OutputSensitivityOptionSetOutputTimeOption

! ************************************************************************** !

subroutine OutputSensitivityOptionSetOutputFormat(output_sensitivity_option,&
                                                option,word)
  ! 
  ! Set the time when to output the sensitivity
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  ! 
  
  use Input_Aux_module
  use Option_module
  
  implicit none
  
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  character(len=MAXWORDLENGTH) :: word
  
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  
  select case(trim(word))
    case('HDF5')
      output_sensitivity_option%format = OUTPUT_HDF5
    case('MATLAB')
      output_sensitivity_option%format = OUTPUT_MATLAB
    case('BINARY')
      output_sensitivity_option%format = OUTPUT_BINARY
    case default
      call InputKeywordUnrecognized(input,word,'SENSITIVITY_RICHARDS,&
                                    &FORMAT',option)
  end select

end subroutine OutputSensitivityOptionSetOutputFormat

! ************************************************************************** !

subroutine OutputSensitivity(J,option,output_sensitivity_option,variable_name)
  ! 
  ! Output the J matrix in the desired format
  ! 
  ! Author: Moise Rousseau
  ! Date: 20 dec 2020
  ! 
  
  use Option_module
  use Output_Aux_module
  use HDF5_module
  
  implicit none
  
  
  Mat :: J
  type(option_type), pointer :: option
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  
  PetscBool :: first
  integer(HID_T) :: file_id
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: variable_name
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  call OutputSensitivityCreateFilename(output_sensitivity_option,option,&
                                       variable_name,filename)
  
  select case(output_sensitivity_option%format)
    case(OUTPUT_HDF5)
      call OutputSensitivityOpenHDF5(option,output_sensitivity_option,&
                                     filename,file_id,first)
      if (first) then
        !call OutputHDF5Provenance(option, output_option, file_id)
        call OutputSensitivityWriteMatrixIJ(J,option,file_id)
      endif
      call OutputSensitivityWriteMatrixData(J,option,file_id)
      call OutputHDF5CloseFile(file_id)
    case(OUTPUT_MATLAB)
      option%io_buffer = '--> writing sensitivity to matlab ascii file: ' // &
                                                     trim(filename)
      call PrintMsg(option)
      call PetscViewerASCIIOpen(option%mycomm,filename, &
                                 viewer,ierr);CHKERRQ(ierr)
      call MatView(J,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerPopFormat(viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    case(OUTPUT_BINARY)
      option%io_buffer = '--> writing sensitivity to matlab binary file: ' // &
                                                     trim(filename)
      call PrintMsg(option)
      call PetscViewerBinaryOpen(option%mycomm,filename, &
                                 FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
      call MatView(J,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerPopFormat(viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  end select
  
end subroutine OutputSensitivity

! ************************************************************************** !

subroutine OutputSensitivityOpenHDF5(option,output_sensitivity_option, &
                                     filename,file_id,first)
  !
  ! Determine the propper hdf5 sensitivity output file name and open it.
  !
  ! Return the file handle and 'first' flag indicating if this is the
  ! first time the file has been opened.
  !
  
  use Option_module
  use Output_Aux_module
  use hdf5
  use String_module

  implicit none

  type(option_type), pointer :: option
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  character(len=MAXSTRINGLENGTH) :: filename
  PetscBool, intent(out) :: first
  integer(HID_T), intent(out) :: file_id
  
  PetscErrorCode :: ierr
  integer(HID_T) :: prop_id
  
  PetscMPIInt :: hdf5_err

  !create filename

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  if (output_sensitivity_option%plot_number > 0) then
    first = PETSC_FALSE
  else
    first = PETSC_TRUE
  endif
  if (.not. first) then
    call h5eset_auto_f(OFF,hdf5_err)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf5_err,prop_id)
    if (hdf5_err /= 0) first = PETSC_TRUE
    call h5eset_auto_f(ON,hdf5_err)
  endif
  if (first) then 
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf5_err, &
                      H5P_DEFAULT_F,prop_id)
  endif
  call h5pclose_f(prop_id,hdf5_err)

  if (first) then
    option%io_buffer = '--> creating sensitivity hdf5 output file: ' // &
                                                     trim(filename)
  else
    option%io_buffer = '--> appending to sensitivity hdf5 output file: ' // &
                                                     trim(filename)
  endif
  call PrintMsg(option)

end subroutine OutputSensitivityOpenHDF5

! ************************************************************************** !

subroutine OutputSensitivityCreateFilename(output_sensitivity_option, &
                                           option,variable_name,filename)
  ! 
  ! Create sensitivity output filename
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  
  use Option_module
  use String_module

  implicit none
  
  type(output_sensitivity_option_type), pointer :: output_sensitivity_option
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: variable_name
  character(len=MAXSTRINGLENGTH), intent(out) :: filename
  
  PetscInt :: file_number
  character(len=MAXSTRINGLENGTH) :: string, string2
  
  if (output_sensitivity_option%format == OUTPUT_HDF5) then
    if (output_sensitivity_option%timestep_per_hdf5_file > 0) &
      file_number = floor(real(output_sensitivity_option%plot_number) / &
                         output_sensitivity_option%timestep_per_hdf5_file)
  else
    file_number = output_sensitivity_option%plot_number
  endif
  if (file_number < 10) then
    write(string,'("00",i1)') file_number
  else if (file_number < 100) then
    write(string,'("0",i2)') file_number  
  else if (file_number < 1000) then
    write(string,'(i3)') file_number  
  else if (file_number < 10000) then
    write(string,'(i4)') file_number
  else if (file_number < 100000) then
    write(string,'(i5)') file_number
  else
    option%io_buffer = 'Plot number exceeds current maximum of 10^5.'
    call PrintErrMsgToDev(option,'ask for a higher maximum')
  endif 
  string = adjustl(string)

  if (output_sensitivity_option%format == OUTPUT_HDF5) then
    if (output_sensitivity_option%timestep_per_hdf5_file > 0) then
      filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                  '-sensitivity-' // trim(string) // '.h5'
    else
      filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                  '-sensitivity.h5'
    endif
  else
    string2 = trim(variable_name)
    call StringToLower(string2)
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                '-sensitivity-' // string2 // '-' // trim(string)
    select case(output_sensitivity_option%format)
      case(OUTPUT_MATLAB)
        filename = trim(filename) // '.mat'
      case(OUTPUT_BINARY)
        filename = trim(filename) // '.bin'
    end select
  endif
 
end subroutine OutputSensitivityCreateFilename

! ************************************************************************** !

subroutine OutputHDF5CloseFile(file_id)

  use hdf5

  implicit none

  integer(HID_T), intent(in) :: file_id
  integer :: hdf5_err

  call h5fclose_f(file_id, hdf5_err)

end subroutine OutputHDF5CloseFile

! ************************************************************************** !

subroutine OutputSensitivityWriteMatrixIJ(J, option, file_id)
  ! 
  ! Write matrix IJ data
  ! 
  ! Author: Moise Rousseau
  ! Date: 20 dec 2020
  ! 

  use Option_module
  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec
  use String_module
  
  Mat :: J
  type(option_type), pointer :: option
  integer(HID_T) :: file_id
  
  Vec :: global_i, global_j, datas
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: info(MAT_INFO_SIZE)
  PetscInt :: mat_non_zeros_global, count
  PetscInt :: nrows, irow, ncols, icol
  PetscInt, allocatable :: cols(:)
  PetscReal, allocatable :: values(:)
  PetscInt :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  integer(HID_T) :: grp_id
  PetscMPIInt :: hdf5_err
  
  !create mpi output vec
  !MAT_GLOBAL=2
  call MatGetInfo(J, 2, info, ierr);CHKERRQ(ierr)
  mat_non_zeros_global = info(MAT_INFO_NZ_ALLOCATED)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, mat_non_zeros_global, &
                    global_i,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, mat_non_zeros_global, &
                    global_j,ierr);CHKERRQ(ierr)
  
  !fill them
  call MatGetSize(J,nrows,ncols,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_i,vec_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(global_j,vec_ptr2,ierr);CHKERRQ(ierr)
  vec_ptr(:) = UNINITIALIZED_DOUBLE
  vec_ptr2(:) = UNINITIALIZED_DOUBLE
  count = 1
  do irow = 1, nrows
    call MatGetRow(J,irow,ncols,cols,values,ierr);CHKERRQ(ierr)
    do icol = 1, ncols
      vec_ptr(count) = irow
      vec_ptr2(count) = cols(icol)
      count = count + 1
    enddo
  enddo
  call VecRestoreArrayF90(global_i,vec_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(global_j,vec_ptr2,ierr);CHKERRQ(ierr)
  
  !merge all data in a single vector 
  ! not needed but verify this TODO
  
  !output it
  string = "Mat indice I"
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  call HDF5WriteDataSetFromVec(string,option,global_i,grp_id, &
                               H5T_NATIVE_DOUBLE)
  call h5gclose_f(grp_id,hdf5_err)
  string = "Mat indice J"
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  call HDF5WriteDataSetFromVec(string,option,global_j,grp_id, &
                               H5T_NATIVE_DOUBLE)
  call h5gclose_f(grp_id,hdf5_err)
  
  call VecDestroy(global_i,ierr);CHKERRQ(ierr)
  call VecDestroy(global_j,ierr);CHKERRQ(ierr)
                    
end subroutine OutputSensitivityWriteMatrixIJ

! ************************************************************************** !

subroutine OutputSensitivityWriteMatrixData(J,option,file_id)
  ! 
  ! Write matrix data
  ! 
  ! Author: Moise Rousseau
  ! Date: 20 dec 2020
  ! 
  
  use Option_module
  use Output_Aux_module
  use hdf5
  use String_module
  !use HDF5_module
  
  implicit none
  
  Mat :: J
  type(option_type), pointer :: option
  !type(output_option_type), pointer :: output_option
  integer(HID_T) :: file_id
  
  character(len=MAXSTRINGLENGTH) :: string
  integer(HID_T) :: grp_id
  PetscMPIInt :: hdf5_err
  
  ! create a group for the data set TODO review this
!  write(string,'(''Time'',es13.5,x,a1)') &
!        option%time/output_option%tconv,output_option%tunit
!  if (len_trim(output_option%plot_name) > 2) then
!    string = trim(string) // ' ' // output_option%plot_name
!  endif
  !string = trim(string3) // ' ' // trim(string)

  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)

  ! write group attributes
  !call OutputHDF5WriteSnapShotAtts(grp_id,option) TODO

end subroutine OutputSensitivityWriteMatrixData

! ************************************************************************** !
 
end module
