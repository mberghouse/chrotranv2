
module Sensitivity_Output_module

#include "petsc/finclude/petscsys.h"
  use petscsys
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use hdf5
  use PFLOTRAN_Constants_module
  use Sensitivity_Aux_module
  
  implicit none

  private
  
  PetscMPIInt, parameter :: ON=1, OFF=0
  
  public :: OutputSensitivity

contains

! ************************************************************************** !

subroutine OutputSensitivity(J,option,output_option, &
                             sensitivity_output_option,variable,mode)
  ! 
  ! Output the J matrix in the desired format
  ! 
  ! Author: Moise Rousseau
  ! Date: 20 dec 2020
  ! 
  
  use Option_module
  use Output_Aux_module
  use HDF5_module
  use Output_HDF5_module
  
  implicit none
  
  Mat :: J
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  character(len=MAXWORDLENGTH) :: mode
  
  PetscBool :: first
  integer(HID_T) :: file_id
  character(len=MAXSTRINGLENGTH) :: filename
  type(sensitivity_output_variable_type) :: variable
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  call OutputSensitivityCreateFilename(sensitivity_output_option,option,&
                                       variable%name,mode,filename)
  
  select case(sensitivity_output_option%format)
    case(SENSITIVITY_OUTPUT_HDF5)
      call OutputSensitivityOpenHDF5(option,sensitivity_output_option,&
                                     filename,file_id,first)
      if (first) then
        call OutputSensitivityWriteMatrixIJ(J,option,file_id)
      endif
        !call OutputHDF5Provenance(option, output_option, file_id)
      call OutputSensitivityWriteMatrixData(J,option,output_option, &
                                            sensitivity_output_option, &
                                            variable,file_id)
        !call OutputHDF5Provenance(option, output_option, file_id)
      call OutputHDF5CloseFile(option,file_id)
    case(SENSITIVITY_OUTPUT_ASCII)
      option%io_buffer = '--> writing sensitivity to ascii file: ' // &
                                                     trim(filename)
      call PrintMsg(option)
      call PetscViewerASCIIOpen(option%mycomm,filename, &
                                 viewer,ierr);CHKERRQ(ierr)
      call MatView(J,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    case(SENSITIVITY_OUTPUT_MATLAB)
      option%io_buffer = '--> writing sensitivity to matlab ascii file: ' // &
                                                     trim(filename)
      call PrintMsg(option)
      call PetscViewerASCIIOpen(option%mycomm,filename, &
                                 viewer,ierr);CHKERRQ(ierr)
      call PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB, &
                                 ierr);CHKERRQ(ierr)
      call MatView(J,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerPopFormat(viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    case(SENSITIVITY_OUTPUT_BINARY)
      option%io_buffer = '--> writing sensitivity to matlab binary file: ' // &
                                                     trim(filename)
      call PrintMsg(option)
      call PetscViewerBinaryOpen(option%mycomm,filename, &
                                 FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
      call MatView(J,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerPopFormat(viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    case default
  end select
  
end subroutine OutputSensitivity

! ************************************************************************** !

subroutine OutputSensitivityOpenHDF5(option,sensitivity_output_option, &
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
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  character(len=MAXSTRINGLENGTH) :: filename
  PetscBool, intent(out) :: first
  integer(HID_T), intent(out) :: file_id
  
  PetscErrorCode :: ierr
  integer(HID_T) :: prop_id
  PetscMPIInt :: hdf5_err

  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
  if (sensitivity_output_option%first_plot_flag) then
    first = PETSC_TRUE
    sensitivity_output_option%first_plot_flag = PETSC_FALSE
  else
    first = PETSC_FALSE
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

subroutine OutputSensitivityCreateFilename(sensitivity_output_option, &
                                           option,variable_name,mode,filename)
  ! 
  ! Create sensitivity output filename
  ! 
  ! Author: Moise Rousseau
  ! Date: 01/04/2021
  
  use Option_module
  use String_module

  implicit none
  
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: variable_name
  character(len=MAXWORDLENGTH) :: mode
  character(len=MAXSTRINGLENGTH), intent(out) :: filename
  
  PetscInt :: file_number
  character(len=MAXSTRINGLENGTH) :: string, string2
  
  if (sensitivity_output_option%format == SENSITIVITY_OUTPUT_HDF5 .and. &
      sensitivity_output_option%timestep_per_hdf5_file > 0) then
      file_number = floor(real(sensitivity_output_option%plot_number) / &
                         sensitivity_output_option%timestep_per_hdf5_file)
  else
    file_number = sensitivity_output_option%plot_number
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
    !write(option%io_buffer, '(i10)') file_number
    !call PrintMsg(option)
    option%io_buffer = 'Plot number exceeds current maximum of 10^5.'
    call PrintErrMsgToDev(option,'ask for a higher maximum')
  endif 
  string = adjustl(string)

  if (sensitivity_output_option%format == SENSITIVITY_OUTPUT_HDF5) then
    if (sensitivity_output_option%timestep_per_hdf5_file > 0) then
      filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                  '-sensitivity-' // trim(mode) // '-' // trim(string) // '.h5'
    else
      filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                  '-sensitivity-flow.h5'
    endif
  else
    string2 = trim(variable_name)
    call StringToLower(string2)
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
                '-sensitivity-' // trim(mode) // '-'  // trim(string2) // &
                 '-' // trim(string)
    select case(sensitivity_output_option%format)
      case(SENSITIVITY_OUTPUT_MATLAB)
        filename = trim(filename) // '.mat'
      case(SENSITIVITY_OUTPUT_BINARY)
        filename = trim(filename) // '.bin'
      case(SENSITIVITY_OUTPUT_ASCII)
        filename = trim(filename) // '.txt'
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

subroutine OutputSensitivityWriteMatrixIJ(J,option,file_id)
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
  PetscInt :: local_size
  PetscInt :: cols(SENSITIVITY_MAX_FACE_PER_CELL_OUTPUT)
  PetscInt :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  integer(HID_T) :: grp_id
  PetscMPIInt :: hdf5_err
  
  !get matrix size and non zeros
  !MAT_GLOBAL=2
  call MatGetInfo(J, 2, info, ierr);CHKERRQ(ierr)
  mat_non_zeros_global = info(MAT_INFO_NZ_ALLOCATED)
  call MatGetSize(J,nrows,ncols,ierr);CHKERRQ(ierr)
  
  !create mpi output vec
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, mat_non_zeros_global, &
                    global_i,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, mat_non_zeros_global, &
                    global_j,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_i,local_size,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(global_j,local_size,ierr);CHKERRQ(ierr)
  
  !fill them
  if (option%mycommsize == 1) then
		call VecGetArrayF90(global_i,vec_ptr,ierr);CHKERRQ(ierr)
		call VecGetArrayF90(global_j,vec_ptr2,ierr);CHKERRQ(ierr)
		vec_ptr(:) = UNINITIALIZED_DOUBLE
		vec_ptr2(:) = UNINITIALIZED_DOUBLE
		count = 1
		do irow = 0, nrows-1
		  call MatGetRow(J,irow,ncols,cols,PETSC_NULL_SCALAR,ierr);CHKERRQ(ierr)
		  do icol = 1, ncols
		    vec_ptr(count) = dble(irow)
		    vec_ptr2(count) = dble(cols(icol))
		    count = count + 1
		  enddo
		  call MatRestoreRow(J,irow,ncols,cols,PETSC_NULL_SCALAR, &
		                     ierr);CHKERRQ(ierr)
		enddo
		call VecRestoreArrayF90(global_i,vec_ptr,ierr);CHKERRQ(ierr)
		call VecRestoreArrayF90(global_j,vec_ptr2,ierr);CHKERRQ(ierr)
  else
    ! TODO (moise) check output_common.F90 line 477 or test it in parallel
  endif
  
  !merge all data in a single vector 
  ! not needed but verify this TODO
  
  !output it
  string = "Mat Structure"
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  string = "Rows"
  call HDF5WriteDataSetFromVec(string,option,global_i,grp_id, &
                               H5T_NATIVE_INTEGER)
  string = "Columns"
  call HDF5WriteDataSetFromVec(string,option,global_j,grp_id, &
                               H5T_NATIVE_INTEGER)
  call h5gclose_f(grp_id,hdf5_err)
  
  call VecDestroy(global_i,ierr);CHKERRQ(ierr)
  call VecDestroy(global_j,ierr);CHKERRQ(ierr)
                    
end subroutine OutputSensitivityWriteMatrixIJ

! ************************************************************************** !

subroutine OutputSensitivityWriteMatrixData(J,option,output_option, &
                                            sensitivity_output_option, &
                                            variable,file_id)
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
  use HDF5_module, only : HDF5WriteDataSetFromVec
  
  implicit none
  
  Mat :: J
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  type(sensitivity_output_variable_type) :: variable
  integer(HID_T) :: file_id
  
  Vec :: datas
  PetscReal, pointer :: vec_ptr(:)
  PetscReal :: info(MAT_INFO_SIZE)
  PetscInt :: mat_non_zeros_global, count
  PetscInt :: nrows, irow, ncols, icol
  PetscInt :: local_size
  PetscReal :: values(SENSITIVITY_MAX_FACE_PER_CELL_OUTPUT)
  PetscInt :: ierr
  character(len=MAXSTRINGLENGTH) :: string
  integer(HID_T) :: grp_id
  PetscMPIInt :: hdf5_err
  
  !get matrix size and non zeros
  !MAT_GLOBAL=2
  call MatGetInfo(J, 2, info, ierr);CHKERRQ(ierr)
  mat_non_zeros_global = info(MAT_INFO_NZ_ALLOCATED)
  call MatGetSize(J,nrows,ncols,ierr);CHKERRQ(ierr)
  
  !create mpi output vec
  call VecCreateMPI(option%mycomm,PETSC_DECIDE, mat_non_zeros_global, &
                    datas,ierr);CHKERRQ(ierr)
  call VecGetLocalSize(datas,local_size,ierr);CHKERRQ(ierr)
  
  !fill them
  if (option%mycommsize == 1) then
		call VecGetArrayF90(datas,vec_ptr,ierr);CHKERRQ(ierr)
		vec_ptr(:) = UNINITIALIZED_DOUBLE
		count = 1
		do irow = 0, nrows-1
		  call MatGetRow(J,irow,ncols,PETSC_NULL_INTEGER,values,ierr);CHKERRQ(ierr)
		  do icol = 1, ncols
		    vec_ptr(count) = dble(values(icol))
		    count = count + 1
		  enddo
		  call MatRestoreRow(J,irow,ncols,PETSC_NULL_INTEGER,values, &
		                     ierr);CHKERRQ(ierr)
		enddo
		call VecRestoreArrayF90(datas,vec_ptr,ierr);CHKERRQ(ierr)
  else
    ! TODO (moise) check output_common.F90 line 477 or test it in parallel
  endif
  
  !open or create the corresponding time group
  write(string,'(''Time:'',es13.5,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)
  !output dataset
  string = trim(variable%name) // " [" // trim(variable%units) // ']'
  call HDF5WriteDataSetFromVec(string,option,datas,grp_id, &
                               H5T_NATIVE_DOUBLE)
  call h5gclose_f(grp_id,hdf5_err)
  
  call VecDestroy(datas,ierr);CHKERRQ(ierr)

end subroutine OutputSensitivityWriteMatrixData

! ************************************************************************** !
 
end module
