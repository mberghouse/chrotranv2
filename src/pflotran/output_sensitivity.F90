
module Sensitivity_Output_module

#include "petsc/finclude/petscsys.h"
  use petscsys
#include "petsc/finclude/petscsnes.h"
  use petscsnes
#include "petsc/finclude/petscmat.h"
  use petscmat
  use hdf5
  use PFLOTRAN_Constants_module
  use Sensitivity_Aux_module
  
  implicit none

  private
  
  PetscMPIInt, parameter :: ON=1, OFF=0
  
  public :: OutputSensitivity

contains

! ************************************************************************** !

subroutine OutputSensitivity(J,grid,option,output_option, &
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
  use Grid_module
  
  implicit none
  
  Mat :: J
  type(grid_type), pointer :: grid
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
        call OutputSensitivityWriteMatrixIJ(J,grid,option,file_id)
        !call OutputHDF5Provenance(option, output_option, file_id)
      endif
      call OutputSensitivityWriteMatrixData(J,grid,option,output_option, &
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

subroutine OutputSensitivityWriteMatrixIJ(J,grid,option,file_id)
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
  use Grid_module
  use Connection_module
  
  Mat :: J
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(connection_set_type), pointer :: cur_connection_set
  integer(HID_T) :: file_id
  
  !get indices
  PetscInt, allocatable :: local_i(:), local_j(:)
  PetscReal :: info(MAT_INFO_SIZE)
  PetscInt :: mat_non_zeros_local, mat_non_zeros_global
  PetscInt :: local_id, ghosted_id, natural_id, count
  PetscInt :: nlmax, iconn, id_up, id_dn
  PetscInt :: ierr
  !save indices
  character(len=MAXSTRINGLENGTH) :: string
  integer(HID_T) :: file_space_id,memory_space_id, data_set_id, prop_id
  PetscMPIInt :: rank_mpi
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  integer(HID_T) :: grp_id
  PetscMPIInt :: hdf5_err
  PetscInt :: istart
  
  cur_connection_set => grid%internal_connection_set_list%first
  
  !get matrix size for output vector
  !MAT_LOCAL = 1
  call MatGetInfo(J,1,info,ierr);CHKERRQ(ierr)
  mat_non_zeros_local = info(MAT_INFO_NZ_ALLOCATED)
  !MAT_SUM_GLOBAL = 3
  call MatGetInfo(J, 3, info, ierr);CHKERRQ(ierr)
  mat_non_zeros_global = info(MAT_INFO_NZ_ALLOCATED)
  
  !create output array
  allocate(local_i(mat_non_zeros_local))
  allocate(local_j(mat_non_zeros_local))
  local_i = UNINITIALIZED_INTEGER
  local_j = UNINITIALIZED_INTEGER
  if (associated(grid%unstructured_grid)) then 
    nlmax = grid%unstructured_grid%nlmax
  else
    nlmax = grid%structured_grid%nlmax
  endif
  
  count = 1
  !diag term
  do local_id = 1, nlmax
    ghosted_id = local_id
    natural_id = grid%nG2A(ghosted_id)
    local_i(count) = natural_id
    local_j(count) = natural_id
    count = count + 1
  enddo
  
  !neighbors term
  cur_connection_set => grid%internal_connection_set_list%first
  do
    if (.not. associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      id_up = cur_connection_set%id_up(iconn)
      id_dn = cur_connection_set%id_dn(iconn)
      !get values
      if (id_up <= nlmax) then ! local
        natural_id = grid%nG2A(id_up)
        local_i(count) = natural_id
        natural_id = grid%nG2A(id_dn)
        local_j(count) = natural_id
        count = count + 1
      endif
      if (id_dn <= nlmax) then ! local
        natural_id = grid%nG2A(id_dn)
        local_i(count) = natural_id
        natural_id = grid%nG2A(id_up)
        local_j(count) = natural_id
        count = count + 1
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
  
	
  !output it
  !create group
  string = "Mat Structure"
  call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)

  ! -- ROWS -- !
  ! Ask for space and organize it
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims(1) = mat_non_zeros_local
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
  ! file space which is a 2D block
  dims(1) = mat_non_zeros_global
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  
  string = "Row Indices" // CHAR(0)
  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(grp_id,string,data_set_id,hdf5_err)
  if (hdf5_err < 0) then
    call h5eset_auto_f(ON,hdf5_err)
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(grp_id,string,H5T_NATIVE_INTEGER,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5eset_auto_f(ON,hdf5_err)
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif
  call h5pclose_f(prop_id,hdf5_err)
  
  !geh: cannot use dims(1) in MPI_Allreduce as it causes errors on 
  !     Juqueen
  istart = 0
  call MPI_Exscan(mat_non_zeros_local, istart, ONE_INTEGER_MPI, MPIU_INTEGER, &
                  MPI_SUM, option%mycomm, ierr);CHKERRQ(ierr)
  start(1) = istart
  length(1) = mat_non_zeros_local
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif
  !call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_INTEGER,local_i,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  !call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr) 
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  ! -- COLUMNS -- !
  ! Ask for space and organize it
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims(1) = mat_non_zeros_local
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
  ! file space which is a 2D block
  dims(1) = mat_non_zeros_global
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  
  string = "Column Indices" // CHAR(0)
  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(grp_id,string,data_set_id,hdf5_err)
  if (hdf5_err < 0) then
    call h5eset_auto_f(ON,hdf5_err)
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(grp_id,string,H5T_NATIVE_INTEGER,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5eset_auto_f(ON,hdf5_err)
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif
  call h5pclose_f(prop_id,hdf5_err)
  
  !geh: cannot use dims(1) in MPI_Allreduce as it causes errors on 
  !     Juqueen
  istart = 0
  call MPI_Exscan(mat_non_zeros_local, istart, ONE_INTEGER_MPI, MPIU_INTEGER, &
                  MPI_SUM, option%mycomm, ierr);CHKERRQ(ierr)
  start(1) = istart
  length(1) = mat_non_zeros_local
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif
  !call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_INTEGER,local_j,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  !call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr) 
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call h5gclose_f(grp_id,hdf5_err)
  deallocate(local_i)
  deallocate(local_j)
                    
end subroutine OutputSensitivityWriteMatrixIJ

! ************************************************************************** !

subroutine OutputSensitivityWriteMatrixData(J,grid,option,output_option, &
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
  use Connection_module
  use HDF5_module, only : HDF5WriteDataSetFromVec
  use Grid_module
  
  implicit none
  
  Mat :: J
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(sensitivity_output_option_type), pointer :: sensitivity_output_option
  type(sensitivity_output_variable_type) :: variable
  integer(HID_T) :: file_id
  
  !compute sensitivity
  PetscReal, allocatable :: datas(:)
  PetscReal :: info(MAT_INFO_SIZE)
  PetscInt :: mat_non_zeros_local, mat_non_zeros_global
  PetscInt :: local_id, count
  PetscInt :: temp_id_in(1), temp_id_out(1), neighbor_id(1)
  PetscInt :: nlmax, iconn, id_up, id_dn
  PetscInt :: ierr
  ISLocalToGlobalMapping :: rmapping,cmapping
  type(connection_set_type), pointer :: cur_connection_set
  !save sensitivity
  character(len=MAXSTRINGLENGTH) :: string
  integer(HID_T) :: file_space_id,memory_space_id, data_set_id, prop_id
  PetscMPIInt :: rank_mpi
  integer(HSIZE_T) :: dims(3)
  integer(HSIZE_T) :: start(3), length(3), stride(3)
  integer(HID_T) :: grp_id
  PetscMPIInt :: hdf5_err
  PetscInt :: istart
  
  !get matrix size for output vector
  !MAT_LOCAL = 1
  call MatGetInfo(J, 1, info, ierr);CHKERRQ(ierr)
  mat_non_zeros_local = info(MAT_INFO_NZ_ALLOCATED)
  !MAT_SUM_GLOBAL = 3
  call MatGetInfo(J, 3, info, ierr);CHKERRQ(ierr)
  mat_non_zeros_global = info(MAT_INFO_NZ_ALLOCATED)
  call MatGetLocalToGlobalMapping(J,rmapping,cmapping,ierr);CHKERRQ(ierr)

  allocate(datas(mat_non_zeros_local))
  datas = UNINITIALIZED_DOUBLE
  if (associated(grid%unstructured_grid)) then 
    nlmax = grid%unstructured_grid%nlmax
  else
    nlmax = grid%structured_grid%nlmax
  endif
  
  count = 1
  !diag term
  do local_id = 1, nlmax
    temp_id_in(1) = local_id - 1
    call ISLocalToGlobalMappingApply(rmapping,1,temp_id_in, &
                                     temp_id_out,ierr);CHKERRQ(ierr)
    call MatGetValues(J,1,temp_id_out(1),1,temp_id_out(1), &
                      datas(count),ierr);CHKERRQ(ierr)
    count = count + 1
  enddo
  
  !neighbors term
  cur_connection_set => grid%internal_connection_set_list%first
  do
    if (.not. associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      id_up = cur_connection_set%id_up(iconn)
      id_dn = cur_connection_set%id_dn(iconn)
      !convert ghosted id to mat local index
      temp_id_in(1) = id_up - 1
      call ISLocalToGlobalMappingApply(rmapping,1,temp_id_in, &
                                       temp_id_out,ierr);CHKERRQ(ierr)
      temp_id_in(1) = id_dn - 1
      call ISLocalToGlobalMappingApply(cmapping,1,temp_id_in, &
                                       neighbor_id,ierr);CHKERRQ(ierr)
      !get values
      if (id_up <= nlmax) then ! local
        call MatGetValues(J,1,temp_id_out(1),1,neighbor_id(1), &
                                 datas(count),ierr);CHKERRQ(ierr)
        count = count + 1
      endif
      if (id_dn <= nlmax) then ! local
        call MatGetValues(J,1,neighbor_id(1),1,temp_id_out(1), &
                                 datas(count),ierr);CHKERRQ(ierr)
        count = count + 1
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
  
  ! -- Output it ---
  !open or create the corresponding time group
  write(string,'(''Time:'',es13.5,x,a1)') &
        option%time/output_option%tconv,output_option%tunit
  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)

  ! Ask for space and organize it
  ! memory space which is a 1D vector
  rank_mpi = 1
  dims(1) = mat_non_zeros_local
  call h5screate_simple_f(rank_mpi,dims,memory_space_id,hdf5_err,dims)
  ! file space which is a 2D block
  dims(1) = mat_non_zeros_global
  call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,hdf5_err)
  
  string = trim(variable%name) // " [" // trim(variable%units) // ']'
  call h5eset_auto_f(OFF,hdf5_err)
  call h5dopen_f(grp_id,string,data_set_id,hdf5_err)
  if (hdf5_err < 0) then
    call h5eset_auto_f(ON,hdf5_err)
    ! if the dataset does not exist, create it
    call h5screate_simple_f(rank_mpi,dims,file_space_id,hdf5_err,dims)
    call h5dcreate_f(grp_id,string,H5T_NATIVE_DOUBLE,file_space_id, &
                     data_set_id,hdf5_err,prop_id)
  else
    call h5eset_auto_f(ON,hdf5_err)
    call h5dget_space_f(data_set_id,file_space_id,hdf5_err)
  endif
  call h5pclose_f(prop_id,hdf5_err)
  
  !geh: cannot use dims(1) in MPI_Allreduce as it causes errors on 
  !     Juqueen
  istart = 0
  call MPI_Exscan(mat_non_zeros_local, istart, ONE_INTEGER_MPI, MPIU_INTEGER, &
                  MPI_SUM, option%mycomm, ierr);CHKERRQ(ierr)
  start(1) = istart
  length(1) = mat_non_zeros_local
  stride = 1
  call h5sselect_hyperslab_f(file_space_id,H5S_SELECT_SET_F,start,length, &
                             hdf5_err,stride,stride)
  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
    call h5pset_dxpl_mpio_f(prop_id,H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
#endif
  !call PetscLogEventBegin(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr)
  call h5dwrite_f(data_set_id,H5T_NATIVE_DOUBLE,datas,dims, &
                  hdf5_err,memory_space_id,file_space_id,prop_id)
  !call PetscLogEventEnd(logging%event_h5dwrite_f,ierr);CHKERRQ(ierr) 
  
  call h5pclose_f(prop_id,hdf5_err)
  call h5dclose_f(data_set_id,hdf5_err)
  call h5sclose_f(file_space_id,hdf5_err)

  call h5gclose_f(grp_id,hdf5_err)
  deallocate(datas)

end subroutine OutputSensitivityWriteMatrixData

! ************************************************************************** !
 
end module
