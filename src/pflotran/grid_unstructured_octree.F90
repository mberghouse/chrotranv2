module Grid_Unstructured_Octree_module

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Geometry_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Explicit_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: UGridOctreeReadVTK, &
            UGridOctreeToExplicit, &
            UGridOctCompCellCentroidVol, &
            UGridOctCompFaceCentroidArea, &
            UGridOctCompFaceCenAreaPlane            
!            UGridOctreeReadHDF5, &
!            UGridOctreeConstruct

contains

! ************************************************************************** !

subroutine UGridOctreeReadVTK(unstructured_grid,filename,option)
  !
  ! Reads an octree unstructured grid vtk file in serial
  !
  ! Author: Bryan He
  ! Date: 07/11/22
  !

  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(unstructured_octree_type), pointer :: octree_grid
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string, hint
  character(len=MAXWORDLENGTH) :: word, card
  PetscInt :: fileid, icell, iconn, irank, remainder, temp_int

  PetscInt :: num_cells, num_connections, num_elems
  PetscInt :: num_cells_local, num_cells_local_save
  PetscInt :: num_connections_local, num_connections_local_save
  PetscInt :: num_elems_local, num_elems_local_save
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscErrorCode :: ierr
  PetscReal, allocatable :: temp_real_array(:,:)
  PetscInt, allocatable :: temp_int_array(:,:)
  PetscInt :: ivertex,ivertex_offset,num_vertices, num_grid_vertices
  PetscInt :: ioff,ip,num_grid_vertices_offset,num_grid_vertices_conn
  PetscInt :: num_line,extre_num_grid_vertices,line,i_offset,i
  PetscReal :: temp_x,temp_y,temp_z,centroids(3),volume
  PetscInt, allocatable :: octree_pts_offset(:),octree_pts_conn(:)
  PetscReal, allocatable:: eightpoints(:,:),fourpoints(:,:),face_centroids(:,:)
  PetscReal, allocatable:: face_area(:),face_coeff(:,:)
  PetscReal, allocatable:: octree_vertex_coordinates(:,:)
  PetscInt, allocatable :: octree_cell_ids(:)
  PetscReal, allocatable:: octree_cell_volumes(:),octree_cell_centroids(:,:)
  PetscReal :: plane(5)
  PetscInt ::  num_face_total,iface,jface
  PetscInt, allocatable :: face_cell(:)
  PetscReal :: centroid_0(3),centroid_1(3)
  PetscReal :: coeff_0(5),coeff_1(5)
  PetscReal :: tol = 0.000001
  PetscInt, allocatable :: conn_id_tmp(:,:),conn_face_id_tmp(:)
  PetscInt, allocatable :: octree_conn_id(:,:)
  PetscReal, allocatable:: octree_conn_face_centroids(:,:),octree_conn_face_area(:)
  PetscInt :: temp_conn_id(2),fdir,i0,i1
  PetscInt, allocatable :: face_corners(:,:)
  PetscInt :: cell0,cell1,num_conn,face_replicate
  PetscReal :: x0,y0,xx_array(4),yy_array(4)
  PetscInt :: ip_tmp,i_tmp
  PetscReal :: p_tmp(3)
  PetscBool :: inbox0,inbox1
  PetscInt :: istart,iend,num_to_send
  PetscInt, allocatable :: my_icell_offset(:),icell_offset(:)
  PetscInt, allocatable :: my_iconn_offset(:),iconn_offset(:)
  octree_grid => unstructured_grid%octree_grid

! -----------------------------------------------------------------
!  print *, 'in UGridOctreeReadVTK '
  call OptionSetBlocking(option,PETSC_FALSE)
  if (OptionIsIORank(option)) then
    fileid = 86
    input => InputCreate(fileid,filename,option)
 
    call InputReadPflotranString(input,option)
    call InputReadCard(input,option,card,PETSC_FALSE)
    call InputReadPflotranString(input,option)
    call InputReadCard(input,option,card,PETSC_FALSE)
    call InputReadPflotranString(input,option)
    call InputReadCard(input,option,card,PETSC_FALSE)
    call InputReadPflotranString(input,option)
    call InputReadCard(input,option,card,PETSC_FALSE)
    word = 'POINTS'
    if (StringCompare(word,card)) then
      card = 'Octree Unstruct. Grid POINTS'
      call InputReadInt(input,option,num_grid_vertices)
!      print *, 'number of points:',num_grid_vertices
      call InputErrorMsg(input,option,'number of points',card)
      call InputReadWord(input,option,word,PETSC_TRUE) ! for the word 'float'
      allocate(octree_vertex_coordinates(3,num_grid_vertices))
      octree_vertex_coordinates(:,:) = 0.0
      hint = 'Octree Unstructured Grid points'
      num_line = num_grid_vertices/3
      extre_num_grid_vertices = MOD(num_grid_vertices,3)
!      print *, 'num_line,extre_num_grid_vertices:',num_line,extre_num_grid_vertices
      line = 0
      do while (line<num_line)
        call InputReadPflotranString(input,option)
        do ip = 1+line*3, 3+line*3
          call InputReadDouble(input,option,temp_x)
          call InputErrorMsg(input,option,'point x coordinate',hint)
          call InputReadDouble(input,option,temp_y)
          call InputErrorMsg(input,option,'point y coordinate',hint)
          call InputReadDouble(input,option,temp_z)
          call InputErrorMsg(input,option,'point z coordinate',hint)
          octree_vertex_coordinates(1,ip) = temp_x
          octree_vertex_coordinates(2,ip) = temp_y
          octree_vertex_coordinates(3,ip) = temp_z
!          print *, '# of point:,',ip,' . x,y,z:',temp_x,temp_y,temp_z
        enddo
        line = line + 1
      enddo
      if (extre_num_grid_vertices>0) then
        call InputReadPflotranString(input,option)
        do ip = 1+num_line*3, extre_num_grid_vertices+num_line*3
           call InputReadDouble(input,option,temp_x)
          call InputErrorMsg(input,option,'point x coordinate',hint)
          call InputReadDouble(input,option,temp_y)
          call InputErrorMsg(input,option,'point y coordinate',hint)
          call InputReadDouble(input,option,temp_z)
          call InputErrorMsg(input,option,'point z coordinate',hint)
          octree_vertex_coordinates(1,ip) = temp_x
          octree_vertex_coordinates(2,ip) = temp_y
          octree_vertex_coordinates(3,ip) = temp_z
!          print *, '# of point:,',ip,' . x,y,z:',temp_x,temp_y,temp_z
        enddo
      endif
    endif
!    print *, 'Done reading points.'
    call InputReadPflotranString(input,option)
    call InputReadCard(input,option,card,PETSC_FALSE)
    word = 'CELLS'
    call InputErrorMsg(input,option,word,card)
    call InputReadInt(input,option,temp_int)
    call InputErrorMsg(input,option,'number of cells',hint)
    num_cells = temp_int - 1
!    print *, 'number of cells:',num_cells
    num_grid_vertices_offset = temp_int
    allocate(octree_pts_offset(num_grid_vertices_offset))
    octree_pts_offset(1:num_grid_vertices_offset) = 0
    num_line = num_grid_vertices_offset/9
    extre_num_grid_vertices = MOD(num_grid_vertices_offset,9)
!    print *, 'num_line,extre_num_grid_vertices:',num_line,extre_num_grid_vertices
    call InputReadInt(input,option,temp_int)
    num_grid_vertices_conn = temp_int
    call InputErrorMsg(input,option,'number of connectivity',hint)
    call InputReadPflotranString(input,option)
    call InputReadCard(input,option,card,PETSC_FALSE)
    word = 'OFFSETS'
    if (StringCompare(word,card)) then
      line = 0
      do while (line<num_line)
        call InputReadPflotranString(input,option)
        do ioff = 1+line*9, 9+line*9
          call InputReadInt(input,option,temp_int)
          call InputErrorMsg(input,option,'offset for points',hint) 
          octree_pts_offset(ioff) = temp_int
!          print *, '# of offset: ',ioff,' . Offset: ',temp_int
        enddo
        line = line + 1
      enddo
      if (extre_num_grid_vertices>0) then
        call InputReadPflotranString(input,option)
        do ioff = 1+num_line*9, extre_num_grid_vertices+num_line*9
          call InputReadInt(input,option,temp_int)
          call InputErrorMsg(input,option,'offset for points',hint)
          octree_pts_offset(ioff) = temp_int
!          print *, '# of offset: ',ioff,' . Offset: ',temp_int
        enddo
      endif
    endif
!    print *, 'Done reading cell offsets.'
    allocate(octree_pts_conn(num_grid_vertices_conn))
    octree_pts_conn(1:num_grid_vertices_conn) = 0
    call InputReadPflotranString(input,option)
    call InputReadCard(input,option,card,PETSC_FALSE)
    word = 'CONNECTIVITY'
    num_line = num_grid_vertices_conn/9
    extre_num_grid_vertices = MOD(num_grid_vertices_conn,9)
!    print *, 'num_line,extre_num_grid_vertices:',num_line,extre_num_grid_vertices
    if (StringCompare(word,card)) then
      line = 0
      do while (line<num_line)
        call InputReadPflotranString(input,option)
        do iconn = 1+line*9, 9+line*9
          call InputReadInt(input,option,temp_int)
          call InputErrorMsg(input,option,'connectivity for points',hint)
          octree_pts_conn(iconn) = temp_int
!          print *, '# of connectivity: ',iconn,' . Connectivity: ',temp_int
        enddo
        line = line + 1
      enddo
      if (extre_num_grid_vertices>0) then
        call InputReadPflotranString(input,option)
        do iconn = 1+num_line*9, extre_num_grid_vertices+num_line*9
          call InputReadInt(input,option,temp_int)
          call InputErrorMsg(input,option,'connectivity for points',hint)
          octree_pts_conn(iconn) = temp_int
!          print *, '# of connectivity: ',iconn,' . Connectivity: ',temp_int
        enddo
      endif
    endif
!    print *, 'Done reading cell connectivity.'

! Now find 1. cell volume and centroids
!          2. face area, centroids and face orientation 
    allocate(octree_cell_ids(num_cells))
    allocate(octree_cell_volumes(num_cells))
    allocate(octree_cell_centroids(3,num_cells))
    octree_cell_ids(:) = 0
    octree_cell_volumes(:) = 0.0
    octree_cell_centroids(:,:) = 0.0
    allocate(eightpoints(3,8))
    allocate(fourpoints(3,4))
    num_face_total = num_cells*6 ! all the faces include replicated ones
    allocate(face_centroids(3,num_face_total))
    allocate(face_area(num_face_total))
    allocate(face_coeff(5,num_face_total))
    allocate(face_cell(num_face_total))
    allocate(face_corners(4,num_face_total))
    face_centroids = 0.0
    face_area = 0.0
    face_coeff = 0.0 
    face_cell = 0
    face_corners = 0
    do icell = 1,num_cells
      octree_cell_ids(icell) = icell
      eightpoints = 0.0
      fourpoints = 0.0
      i_offset = octree_pts_offset(icell)
!       print *, 'icell,i:',icell,i_offset
      do i = 1,8
        i_tmp = i+i_offset
        eightpoints(1:3,i) = octree_vertex_coordinates(1:3,i_tmp)
      enddo
      call UGridOctCompCellCentroidVol(eightpoints ,&
                    volume,centroids)
      octree_cell_volumes(icell) = volume
      octree_cell_centroids(1:3,icell) = centroids(1:3)
!       print *, 'icell:',icell
!       print *, 'volume:',octree_grid%cell_volumes(icell)
!       print *, 'centroids:',octree_grid%cell_centroids(icell)%x,&
!                             octree_grid%cell_centroids(icell)%y,&
!                             octree_grid%cell_centroids(icell)%z
!  west face  
      iface = 1+(icell-1)*6
      fourpoints(:,1) = eightpoints(:,1)
      fourpoints(:,2) = eightpoints(:,3)
      fourpoints(:,3) = eightpoints(:,5)
      fourpoints(:,4) = eightpoints(:,7)     
      call UGridOctCompFaceCenAreaPlane(fourpoints,face_area(iface),&
                            face_centroids(1:3,iface),face_coeff(1:5,iface))
      face_cell(iface) = -icell ! 1 based cell ID, negative sign means the cell
                                 ! is at the downwind side of the face
      face_corners(1:4,iface) = (/1+i_offset,3+i_offset,5+i_offset,7+i_offset/)
!       print *,'iface,area,centroid,face_coeff:',iface,face_area(iface), &
!                             face_centroids(1:3,iface),face_coeff(1:5,iface)
!       print *, 'iface,point index:',iface,face_corners(1:4,iface)
!  east face
      iface = 2+(icell-1)*6
      fourpoints(:,1) = eightpoints(:,2)
      fourpoints(:,2) = eightpoints(:,4)
      fourpoints(:,3) = eightpoints(:,6)
      fourpoints(:,4) = eightpoints(:,8)
      call UGridOctCompFaceCenAreaPlane(fourpoints,face_area(iface),&
                           face_centroids(1:3,iface),face_coeff(1:5,iface))
      face_cell(iface) = icell ! 1 based cell ID, positive sign means the cell
                                 ! is at the upwind side of the face
      face_corners(1:4,iface) = (/2+i_offset,4+i_offset,6+i_offset,8+i_offset/)
!       print *, 'iface,point index:',iface,face_corners(1:4,iface)
!       print *,'iface,area,centroid,face_coeff:',iface,face_area(iface), &
!                             face_centroids(1:3,iface),face_coeff(1:5,iface)
!  south face
      iface = 3+(icell-1)*6
      fourpoints(:,1) = eightpoints(:,1)
      fourpoints(:,2) = eightpoints(:,2)
      fourpoints(:,3) = eightpoints(:,5)
      fourpoints(:,4) = eightpoints(:,6)
      call UGridOctCompFaceCenAreaPlane(fourpoints,face_area(iface),&
                           face_centroids(1:3,iface),face_coeff(1:5,iface))
      face_cell(iface) = -icell ! 1 based cell ID, negative sign means the cell
                                ! is at the downwind side of the face 
      face_corners(1:4,iface) = (/1+i_offset,2+i_offset,5+i_offset,6+i_offset/)
!      print *, 'iface,point index:',iface,face_corners(1:4,iface)
!       print *,'iface,area,centroid,face_coeff:',iface,face_area(iface), &
!                             face_centroids(1:3,iface),face_coeff(1:5,iface)
!  north face
      iface = 4+(icell-1)*6
      fourpoints(:,1) = eightpoints(:,3)
      fourpoints(:,2) = eightpoints(:,4)
      fourpoints(:,3) = eightpoints(:,7)
      fourpoints(:,4) = eightpoints(:,8)
      call UGridOctCompFaceCenAreaPlane(fourpoints,face_area(iface),&
                           face_centroids(1:3,iface),face_coeff(1:5,iface))
      face_cell(iface) = icell ! 1 based cell ID, positive sign means the cell
                                ! is at the upwind side of the face
      face_corners(1:4,iface) = (/3+i_offset,4+i_offset,7+i_offset,8+i_offset/)
!       print *, 'iface,point index:',iface,face_corners(1:4,iface)
!       print *,'iface,area,centroid,face_coeff:',iface,face_area(iface), &
!                             face_centroids(1:3,iface),face_coeff(1:5,iface)
!  botom face
      iface = 5+(icell-1)*6
      fourpoints(:,1) = eightpoints(:,1)
      fourpoints(:,2) = eightpoints(:,2)
      fourpoints(:,3) = eightpoints(:,3)
      fourpoints(:,4) = eightpoints(:,4)
      call UGridOctCompFaceCenAreaPlane(fourpoints,face_area(iface),&
                           face_centroids(1:3,iface),face_coeff(1:5,iface))
      face_cell(iface) = -icell ! 1 based cell ID, negative sign means the cell
                                ! is at the downwind side of the face
      face_corners(1:4,iface) = (/1+i_offset,2+i_offset,3+i_offset,4+i_offset/)
!       print *, 'iface,point index:',iface,face_corners(1:4,iface)
!       print *,'iface,area,centroid,face_coeff:',iface,face_area(iface), &
!                             face_centroids(1:3,iface),face_coeff(1:5,iface)
!  top face
      iface = 6+(icell-1)*6
      fourpoints(:,1) = eightpoints(:,5)
      fourpoints(:,2) = eightpoints(:,6)
      fourpoints(:,3) = eightpoints(:,7)
      fourpoints(:,4) = eightpoints(:,8)
      call UGridOctCompFaceCenAreaPlane(fourpoints,face_area(iface),&
                           face_centroids(1:3,iface),face_coeff(1:5,iface))
      face_cell(iface) = icell ! 1 based cell ID, positive sign means the cell
                                ! is at the upwind side of the face
      face_corners(1:4,iface) = (/5+i_offset,6+i_offset,7+i_offset,8+i_offset/)
!      print *, 'iface,point index:',iface,face_corners(1:4,iface)
!          print *,'iface,area,centroid,face_coeff:',iface,face_area(iface), &
!                             face_centroids(1:3,iface),face_coeff(1:5,iface)
    enddo
    deallocate(eightpoints)
    deallocate(fourpoints)

! Next find overalpping faces
    allocate(conn_id_tmp(2,num_face_total)) ! use maximum face number to allocate the array
    allocate(conn_face_id_tmp(num_face_total))
    conn_id_tmp = 0
    conn_face_id_tmp = 0
    num_conn = 0

    do iface = 1,num_face_total
      centroid_0 = face_centroids(1:3,iface)
      coeff_0 = face_coeff(1:5,iface)
!      do jface = 1,num_face_total
      do jface = iface+1,num_face_total ! reduce the number of iterations
        centroid_1 = face_centroids(1:3,jface)
        coeff_1 = face_coeff(1:5,jface)
        if (all(abs(coeff_0(1:4)-coeff_1(1:4))<tol)) then
          ! find the faces have the same plane coefficients
          fdir = int(coeff_0(5)) 
          if (fdir == 1) then ! x face
            i0 = 2
            i1 = 3
          elseif (fdir == 2) then ! y face
            i0 = 1
            i1 = 3
          elseif (fdir == 3) then ! z face
            i0 = 1
            i1 = 2
          endif
            ! if centroid of i face is in j face
          x0 = face_centroids(i0,iface)
          y0 = face_centroids(i1,iface) 
          do i = 1,4
            ip_tmp = face_corners(i,jface)
            p_tmp(1:3) = octree_vertex_coordinates(1:3,ip_tmp)
            xx_array(i) = p_tmp(i0)
            yy_array(i) = p_tmp(i1) 
          enddo
          inbox0 = GeometryPointInRectangle1(x0,y0,xx_array,yy_array)
            ! if centroid of j face is in i face
          x0 = face_centroids(i0,jface)
          y0 = face_centroids(i1,jface)
          do i = 1,4
            ip_tmp = face_corners(i,iface)
            p_tmp(1:3) = octree_vertex_coordinates(1:3,ip_tmp)
            xx_array(i) = p_tmp(i0)
            yy_array(i) = p_tmp(i1)
          enddo
          inbox1 = GeometryPointInRectangle1(x0,y0,xx_array,yy_array)
          if (inbox0==.True. .or. inbox1==.True.) then
            ! check if centroids are in each other's plane,
            cell0 = face_cell(iface)
            cell1 = face_cell(jface)
            temp_conn_id(1) = cell0
            temp_conn_id(2) = -cell1
            num_conn = num_conn + 1
            conn_id_tmp(1:2,num_conn) = temp_conn_id(1:2)
            if (inbox0 ==.True.) then
              conn_face_id_tmp(num_conn) = iface
            else
              conn_face_id_tmp(num_conn) = jface
            endif
!              print *, '# of conn:',num_conn
!              print *, 'up and downwind cells:',conn_id_tmp(1:2,num_conn)
!              print *, '# of face:',conn_face_id_tmp(num_conn)
          endif 
        endif
      enddo
    enddo
    deallocate(face_coeff)
    deallocate(face_cell)
    deallocate(face_corners)
!    print *, 'Total number of connectivity:',num_conn
    allocate(octree_conn_id(2,num_conn))
    allocate(octree_conn_face_centroids(3,num_conn))
    allocate(octree_conn_face_area(num_conn))
    octree_conn_id(:,:) = conn_id_tmp(:,1:num_conn)
    do iconn = 1, num_conn
      iface = conn_face_id_tmp(iconn)
      octree_conn_face_centroids(1:3,iconn) = face_centroids(1:3,iface)
      octree_conn_face_area(iconn) = face_area(iface)
    enddo
    deallocate(face_centroids)
    deallocate(face_area)
    deallocate(conn_id_tmp)
    deallocate(conn_face_id_tmp)
!    print *, 'Total number of cells:',num_cells
!    do icell = 1,num_cells
!       print *, icell,octree_cell_centroids(:,icell),octree_cell_volumes(icell)
!    enddo
!    print *, 'Total number of connectivity:',num_conn
!    do iconn = 1,num_conn 
!       print *,octree_conn_id(:,iconn),octree_conn_face_centroids(:,iconn),octree_conn_face_area(iconn)
!    enddo
  endif
  call OptionSetBlocking(option,PETSC_TRUE)
  call OptionCheckNonBlockingError(option)

! Now boardcast number of cells to all procs
  call MPI_Bcast(num_cells,ONE_INTEGER_MPI,MPI_INTEGER, &
                 option%driver%io_rank,option%mycomm,ierr)
  octree_grid%num_cells_global = num_cells
  octree_grid%num_vertices = num_cells*8
! divide cells across ranks
  num_cells_local = num_cells/option%comm%mycommsize
  num_cells_local_save = num_cells_local
  remainder = num_cells - &
              num_cells_local*option%comm%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                 num_cells_local + 1
  octree_grid%num_vertices_local = num_cells_local*8
!  print *,'myrank:',option%myrank,'num_cells_local:',num_cells_local
  ! Find istart as offset
  istart = 0
  call MPI_Exscan(num_cells_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
  allocate(my_icell_offset(option%comm%mycommsize))
  allocate(icell_offset(option%comm%mycommsize))
  my_icell_offset = 0
  icell_offset = 0
  my_icell_offset(option%myrank+1) = istart
  call MPI_ALLREDUCE(my_icell_offset,icell_offset,option%comm%mycommsize, &
                     MPI_INT,MPI_SUM,option%mycomm,ierr)
  deallocate(my_icell_offset)
!  print *, 'myrank:',option%myrank,'cell offset:',icell_offset
  allocate(octree_grid%cell_ids(num_cells_local))
  octree_grid%cell_ids = 0
  allocate(octree_grid%cell_volumes(num_cells_local))
  octree_grid%cell_volumes = 0
  allocate(octree_grid%cell_centroids(num_cells_local))
  do icell = 1, num_cells_local
    octree_grid%cell_centroids(icell)%x = 0.d0
    octree_grid%cell_centroids(icell)%y = 0.d0
    octree_grid%cell_centroids(icell)%z = 0.d0
  enddo
  allocate(octree_grid%vertex_coordinates(num_cells_local*8))
  do ivertex = 1, num_cells_local*8
    octree_grid%vertex_coordinates(ivertex)%x = 0.d0
    octree_grid%vertex_coordinates(ivertex)%y = 0.d0
    octree_grid%vertex_coordinates(ivertex)%z = 0.d0
  enddo
  ! Now pass cells from io rank to all others 
  call OptionSetBlocking(option,PETSC_FALSE)
  if (OptionIsIORank(option)) then
!    allocate(temp_real_array(5,num_cells_local_save+1))
    allocate(temp_real_array(29,num_cells_local_save+1))
    do irank = 0, option%comm%mycommsize-1   
      temp_real_array = UNINITIALIZED_DOUBLE
      num_to_send = num_cells_local_save
      istart = icell_offset(irank+1)
!      print *, 'irank,istart:',irank,istart
      if (irank < remainder) num_to_send = num_to_send + 1
      do icell = 1, num_to_send
        temp_real_array(1,icell) = dble(octree_cell_ids(icell+istart))
        temp_real_array(2:4,icell) = octree_cell_centroids(1:3,icell+istart)
        temp_real_array(5,icell) = octree_cell_volumes(icell+istart)         
        do ivertex = 1, 8
          ivertex_offset = (ivertex-1)*3
          temp_real_array(6+ivertex_offset:8+ivertex_offset,icell) = & 
                        octree_vertex_coordinates(1:3,(istart+icell-1)*8+ivertex)
!          print *, 'irank,icell,ivertex:',irank,icell,ivertex
!          print *, 'array index:',6+ivertex_offset,8+ivertex_offset
!          print *, 'vertex index:',(istart+icell-1)*8+ivertex
!          print *, 'x,y,z:',octree_vertex_coordinates(1:3,(istart+icell-1)*8+ivertex)
        enddo
      enddo        
        ! if the cells reside on io_rank
      if (irank == option%driver%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_cells_local
        string = trim(adjustl(string)) // ' cells stored on p0'
        print *, trim(string)
#endif
        do icell = 1, num_cells_local
          octree_grid%cell_ids(icell) = int(temp_real_array(1,icell))
          octree_grid%cell_centroids(icell)%x = temp_real_array(2,icell)
          octree_grid%cell_centroids(icell)%y = temp_real_array(3,icell)
          octree_grid%cell_centroids(icell)%z = temp_real_array(4,icell)
          octree_grid%cell_volumes(icell) = temp_real_array(5,icell)
          do ivertex = 1, 8
            ivertex_offset = (ivertex-1)*3
            octree_grid%vertex_coordinates((icell-1)*8+ivertex)%x = &
                                 temp_real_array(6+ivertex_offset,icell)
            octree_grid%vertex_coordinates((icell-1)*8+ivertex)%y = &
                                 temp_real_array(7+ivertex_offset,icell)
            octree_grid%vertex_coordinates((icell-1)*8+ivertex)%z = &
                                 temp_real_array(8+ivertex_offset,icell)
          enddo          
        enddo
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_send
        write(word,*) irank
        string = trim(adjustl(string)) // ' cells sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif        
        int_mpi = num_to_send*29
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                    num_to_send,option%mycomm,ierr)
      endif
    enddo
  else
     ! other ranks post the recv
#if UGRID_DEBUG
    write(string,*) num_cells_local
    write(word,*) option%myrank
    string = trim(adjustl(string)) // ' cells received from p0 at p' // &
              trim(adjustl(word))
    print *, trim(string)
#endif
    allocate(temp_real_array(29,num_cells_local))
    int_mpi = num_cells_local*29
    call MPI_Recv(temp_real_array,int_mpi, &
                 MPI_DOUBLE_PRECISION,option%driver%io_rank, &
                 MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
    do icell = 1, num_cells_local
      octree_grid%cell_ids(icell) = int(temp_real_array(1,icell))
      octree_grid%cell_centroids(icell)%x = temp_real_array(2,icell)
      octree_grid%cell_centroids(icell)%y = temp_real_array(3,icell)
      octree_grid%cell_centroids(icell)%z = temp_real_array(4,icell)
      octree_grid%cell_volumes(icell) = temp_real_array(5,icell)
      do ivertex = 1, 8
        ivertex_offset = (ivertex-1)*3
        octree_grid%vertex_coordinates((icell-1)*8+ivertex)%x = &
                                         temp_real_array(6+ivertex_offset,icell)
        octree_grid%vertex_coordinates((icell-1)*8+ivertex)%y = &
                                         temp_real_array(7+ivertex_offset,icell)
        octree_grid%vertex_coordinates((icell-1)*8+ivertex)%z = &
                                         temp_real_array(8+ivertex_offset,icell)
      enddo
    enddo
  endif
  deallocate(icell_offset)
  if (OptionIsIORank(option)) then
    deallocate(octree_cell_ids)
    deallocate(octree_cell_centroids) 
    deallocate(octree_cell_volumes)
    deallocate(octree_vertex_coordinates)
  endif
!  do irank = 0,option%comm%mycommsize-1
!  if (option%myrank==irank) then
!  do i = 1,num_cells_local
!     print *, 'myrank:',option%myrank,' icell:',i,'ID:',octree_grid%cell_ids(i)
!     print *, '       centroids:',octree_grid%cell_centroids(i)
!     print *, '       volume:',octree_grid%cell_volumes(i)   
!     print *, '       vertex coord:'
!     do ivertex = 1,8
!       print *, 'iv=',ivertex,'x,y,z:',octree_grid%vertex_coordinates((i-1)*8+ivertex)%x,&
!                                       octree_grid%vertex_coordinates((i-1)*8+ivertex)%y,&
!                                       octree_grid%vertex_coordinates((i-1)*8+ivertex)%z
!     enddo
!  enddo
!  endif
!  enddo
  call OptionSetBlocking(option,PETSC_TRUE)
  call OptionCheckNonBlockingError(option)
  deallocate(temp_real_array)

! Now do the same for the connections
  call MPI_Bcast(num_conn,ONE_INTEGER_MPI,MPI_INTEGER, &
                 option%driver%io_rank,option%mycomm,ierr)
  num_connections = num_conn
  octree_grid%num_connections_global = num_connections
   ! divide cells across ranks
  num_connections_local = num_connections/option%comm%mycommsize
  num_connections_local_save = num_connections_local
  remainder = num_connections - &
              num_connections_local*option%comm%mycommsize
  if (option%myrank < remainder) num_connections_local = &
                                 num_connections_local + 1

  allocate(octree_grid%connections(2,num_connections_local))
  octree_grid%connections = 0
  allocate(octree_grid%face_areas(num_connections_local))
  octree_grid%face_areas = 0
  allocate(octree_grid%face_centroids(num_connections_local))
  do iconn = 1, num_connections_local
    octree_grid%face_centroids(iconn)%x = 0.d0
    octree_grid%face_centroids(iconn)%y = 0.d0
    octree_grid%face_centroids(iconn)%z = 0.d0
  enddo
  
  ! Find istart
  istart = 0
  call MPI_Exscan(num_connections_local, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)
!  print *,'myrank:',option%myrank,'istart:',istart
  allocate(my_iconn_offset(option%comm%mycommsize))
  allocate(iconn_offset(option%comm%mycommsize))
  my_iconn_offset = 0
  iconn_offset = 0
  my_iconn_offset(option%myrank+1) = istart
  call MPI_ALLREDUCE(my_iconn_offset,iconn_offset,option%comm%mycommsize, &
                     MPI_INT,MPI_SUM,option%mycomm,ierr)
  deallocate(my_iconn_offset)
  ! Now pass connections from io rank to all others
  call OptionSetBlocking(option,PETSC_FALSE)
  if (OptionIsIORank(option)) then
    allocate(temp_real_array(6,num_connections_local_save+1))
    do irank = 0, option%comm%mycommsize-1
      temp_real_array = UNINITIALIZED_DOUBLE
      num_to_send = num_connections_local_save
      istart = iconn_offset(irank+1)
      if (irank < remainder) num_to_send = num_to_send + 1
      do iconn = 1, num_to_send
        temp_real_array(1:2,iconn) = dble(octree_conn_id(1:2,iconn+istart))
        temp_real_array(3:5,iconn) = octree_conn_face_centroids(1:3,iconn+istart)
        temp_real_array(6,iconn) = octree_conn_face_area(iconn+istart)
      enddo
      ! if the cells reside on io_rank
      if (irank == option%driver%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_connections_local
        string = trim(adjustl(string)) // ' connections stored on p0'
        print *, trim(string)
#endif
        do iconn = 1, num_connections_local
          octree_grid%connections(1,iconn) = int(temp_real_array(1,iconn))
          octree_grid%connections(2,iconn) = int(temp_real_array(2,iconn))
          octree_grid%face_centroids(iconn)%x = temp_real_array(3,iconn)
          octree_grid%face_centroids(iconn)%y = temp_real_array(4,iconn)
          octree_grid%face_centroids(iconn)%z = temp_real_array(5,iconn)
          octree_grid%face_areas(iconn) = temp_real_array(6,iconn)
        enddo
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_send
        write(word,*) irank
        string = trim(adjustl(string)) // ' connections sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_send*6
        call MPI_Send(temp_real_array,int_mpi,MPI_DOUBLE_PRECISION,irank, &
                      num_to_send,option%mycomm,ierr)
      endif
    enddo
  else        
    ! other ranks post the recv
#if UGRID_DEBUG
    write(string,*) num_connections_local
    write(word,*) option%myrank
    string = trim(adjustl(string)) // ' connections received from p0 at p' // &
              trim(adjustl(word))
    print *, trim(string)
#endif
    allocate(temp_real_array(6,num_connections_local))
    int_mpi = num_connections_local*6
    call MPI_Recv(temp_real_array,int_mpi, &
                  MPI_DOUBLE_PRECISION,option%driver%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
    do iconn = 1, num_connections_local
      octree_grid%connections(1,iconn) = int(temp_real_array(1,iconn))
      octree_grid%connections(2,iconn) = int(temp_real_array(2,iconn))
      octree_grid%face_centroids(iconn)%x = temp_real_array(3,iconn)
      octree_grid%face_centroids(iconn)%y = temp_real_array(4,iconn)
      octree_grid%face_centroids(iconn)%z = temp_real_array(5,iconn)
      octree_grid%face_areas(iconn) = temp_real_array(6,iconn)
    enddo
  endif
!  do irank = 0,option%comm%mycommsize-1
!  if (option%myrank==irank) then
!  do i = 1,num_connections_local
!     print *, 'myrank:',option%myrank,' icell:',i,'ID:',octree_grid%connections(:,i)
!     print *, '       centroids:',octree_grid%face_centroids(i)
!     print *, '       area:',octree_grid%face_areas(i)
!  enddo
!  endif
!  enddo
  deallocate(iconn_offset)
  deallocate(temp_real_array)
  if (OptionIsIORank(option)) then
    deallocate(octree_conn_id)
    deallocate(octree_conn_face_centroids)
    deallocate(octree_conn_face_area)
    call InputDestroy(input)
  endif
  call OptionSetBlocking(option,PETSC_TRUE)
  call OptionCheckNonBlockingError(option)
!  print *, 'rank',option%myrank,' done reading'
!  stop 
end subroutine UGridOctreeReadVTK

subroutine UGridOctreeToExplicit(unstructured_grid,option)
  !
  ! Currently convert all the octree grids to 
  ! explicit grids 
  !
  ! Author: Bryan He
  ! Date: 08/01/22
  !

  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(grid_unstructured_type) :: unstructured_grid
  type(unstructured_explicit_type), pointer :: explicit_grid
  type(unstructured_octree_type), pointer :: octree_grid
  type(option_type) :: option

  PetscInt :: icell, iconn, remainder
  PetscInt :: num_cells, num_connections
  PetscInt :: num_cells_local, num_cells_local_save
  PetscInt :: num_connections_local, num_connections_local_save
  PetscInt :: ivertex, num_vertices, num_grid_vertices
  PetscInt :: ivertex_id
  call OptionSetBlocking(option,PETSC_FALSE)
!  print *, 'rank',option%myrank,' in UGridOctreeToExplicit'
  explicit_grid => unstructured_grid%explicit_grid
  octree_grid => unstructured_grid%octree_grid

! copy all the octree cells infor into explicit grid type data
  explicit_grid%num_cells_global = octree_grid%num_cells_global
  num_cells = explicit_grid%num_cells_global
  explicit_grid%num_vertices = octree_grid%num_vertices
  explicit_grid%num_vertices_local = octree_grid%num_vertices_local 
  ! divide cells across ranks
  num_cells_local = num_cells/option%comm%mycommsize
  num_cells_local_save = num_cells_local
  remainder = num_cells - &
              num_cells_local*option%comm%mycommsize
  if (option%myrank < remainder) num_cells_local = &
                                 num_cells_local + 1
  allocate(explicit_grid%cell_ids(num_cells_local))
  allocate(explicit_grid%cell_volumes(num_cells_local))
  allocate(explicit_grid%cell_centroids(num_cells_local))
  allocate(explicit_grid%vertex_coordinates(num_cells_local*8))
  do icell = 1, num_cells_local
    explicit_grid%cell_ids(icell) = octree_grid%cell_ids(icell)
    explicit_grid%cell_volumes(icell) = octree_grid%cell_volumes(icell)
    explicit_grid%cell_centroids(icell)%x = octree_grid%cell_centroids(icell)%x
    explicit_grid%cell_centroids(icell)%y = octree_grid%cell_centroids(icell)%y
    explicit_grid%cell_centroids(icell)%z = octree_grid%cell_centroids(icell)%z
    do ivertex = 1, 8
      ivertex_id = (icell-1)*8+ivertex
      explicit_grid%vertex_coordinates(ivertex_id)%x = &
                    octree_grid%vertex_coordinates(ivertex_id)%x
      explicit_grid%vertex_coordinates(ivertex_id)%y = &
                    octree_grid%vertex_coordinates(ivertex_id)%y
      explicit_grid%vertex_coordinates(ivertex_id)%z = &
                    octree_grid%vertex_coordinates(ivertex_id)%z
    enddo
  enddo
!  print *, 'rank',option%myrank,' in UGridOctreeToExplicit done cells'
! copy all the octree connections infor into explicit type data
  num_connections = octree_grid%num_connections_global
   ! divide conns across ranks
  num_connections_local = num_connections/option%comm%mycommsize
  num_connections_local_save = num_connections_local
  remainder = num_connections - &
              num_connections_local*option%comm%mycommsize
  if (option%myrank < remainder) num_connections_local = &
                                 num_connections_local + 1

  allocate(explicit_grid%connections(2,num_connections_local))
  allocate(explicit_grid%face_areas(num_connections_local))
  allocate(explicit_grid%face_centroids(num_connections_local))
  do iconn = 1, num_connections_local
    explicit_grid%connections(1:2,iconn) = octree_grid%connections(1:2,iconn)
    explicit_grid%face_areas(iconn) = octree_grid%face_areas(iconn)
    explicit_grid%face_centroids(iconn)%x = octree_grid%face_centroids(iconn)%x
    explicit_grid%face_centroids(iconn)%y = octree_grid%face_centroids(iconn)%y
    explicit_grid%face_centroids(iconn)%z = octree_grid%face_centroids(iconn)%z
  enddo
!  print *, 'rank',option%myrank,' in UGridOctreeToExplicit done conns'
  call OptionSetBlocking(option,PETSC_TRUE)
  call OptionCheckNonBlockingError(option)

  call UGridOctreeDestroy(octree_grid)
end subroutine UGridOctreeToExplicit

subroutine UGridOctCompCellCentroidVol(coordinates,volume,centroid)
  ! Determine the cell centroid and volume
  ! from eight corner points
  ! 
  ! Author: Bryan He
  ! Date: 07/18/2022
  !
  implicit none

  PetscReal :: coordinates(:,:)
  PetscReal :: x_min,y_min,z_min,x_max,y_max,z_max
  PetscInt :: i, num_coordinates
  PetscReal :: xx_array(8), yy_array(8),zz_array(8)
  PetscReal :: volume,centroid(3)
  num_coordinates = size(coordinates,2)
  if (num_coordinates .NE. 8) then
    print *, 'Error: Exactly eight points are needed to compute the cell centoids and volume.'
    stop
  endif

  do i = 1,num_coordinates
    xx_array(i) = coordinates(1,i)
    yy_array(i) = coordinates(2,i)
    zz_array(i) = coordinates(3,i)
  enddo

  x_min = minval(xx_array)
  y_min = minval(yy_array)
  z_min = minval(zz_array)
  x_max = maxval(xx_array)
  y_max = maxval(yy_array)
  z_max = maxval(zz_array)

  volume = (x_max-x_min)*(y_max-y_min)*(z_max-z_min)
 
  centroid(1) = 0.5*(x_min+x_max)
  centroid(2) = 0.5*(y_min+y_max)
  centroid(3) = 0.5*(z_min+z_max)

  
end subroutine UGridOctCompCellCentroidVol

subroutine UGridOctCompFaceCentroidArea(coordinates,area,centroid)
  ! Determine the face centroid and area
  ! from four corner points
  ! 
  ! Author: Bryan He
  ! Date: 07/18/2022
  !
  implicit none

  PetscReal :: coordinates(:,:)
  PetscReal :: x_min,y_min,z_min,x_max,y_max,z_max
  PetscInt :: i, num_coordinates
  PetscReal :: xx_array(4), yy_array(4),zz_array(4)
  PetscReal :: area,centroid(3),dx,dy,dz
  num_coordinates = size(coordinates,2)
  if (num_coordinates .NE. 4) then
    print *, 'Error: Exactly four points are needed to compute the face centoids and area.'
    stop
  endif

  do i = 1,num_coordinates
    xx_array(i) = coordinates(1,i)
    yy_array(i) = coordinates(2,i)
    zz_array(i) = coordinates(3,i)
  enddo

  x_min = minval(xx_array)
  y_min = minval(yy_array)
  z_min = minval(zz_array)
  x_max = maxval(xx_array)
  y_max = maxval(yy_array)
  z_max = maxval(zz_array)
  
  dx = x_max-x_min
  dy = y_max-y_min
  dz = z_max-z_min
  if (dx==0) then
    area=dy*dz
  else if (dy==0) then
    area=dx*dz
  else if (dz==0) then
    area=dx*dy
  endif
  centroid(1) = 0.5*(x_min+x_max)
  centroid(2) = 0.5*(y_min+y_max)
  centroid(3) = 0.5*(z_min+z_max)

  
end subroutine UGridOctCompFaceCentroidArea

subroutine GeometryComputePlaneWithPoints1s(plane,point1,point2,point3)
  !
  ! Calculates the plane defined by three points
  ! 
  ! Modified by Bryan He to (1) include the face orientation information
  !                             plane(5): 1-X face
  !                                       2-Y face
  !                                       3-Z face
  !                         (2) normalize the coefficients by the 1st non-zero
  !                             coefficient
  ! Date: 07/19/2022
  ! Author: Glenn Hammond
  ! Date: 10/30/09, 02/01/17
  !

  implicit none

!  type(plane_type) :: plane
  PetscReal :: point1(3),point2(3),point3(3)
  PetscReal :: x1,y1,z1
  PetscReal :: x2,y2,z2
  PetscReal :: x3,y3,z3
  PetscReal :: plane(5)

  PetscReal :: x12, y12, z12
  PetscReal :: x13, y13, z13
  PetscInt :: i
  PetscReal :: coeff_scale
  x1 = point1(1)
  y1 = point1(2)
  z1 = point1(3)
  x2 = point2(1)
  y2 = point2(2)
  z2 = point2(3)
  x3 = point3(1)
  y3 = point3(2)
  z3 = point3(3)
  x12 = x2-x1
  y12 = y2-y1
  z12 = z2-z1
  x13 = x3-x1
  y13 = y3-y1
  z13 = z3-z1
  plane(1) = y12*z13-z12*y13
  plane(2) = z12*x13-x12*z13
  plane(3) = x12*y13-y12*x13
  plane(4) = -1.d0*(plane(1)*x1+plane(2)*y1+plane(3)*z1)
  do i = 1,3
    if (plane(i) .NE. 0.0) then
      coeff_scale = plane(i)
      plane(5) = float(i)
      exit
    endif
  enddo
  plane(1:4) = plane(1:4)/coeff_scale
   
end subroutine GeometryComputePlaneWithPoints1s

subroutine UGridOctCompFaceCenAreaPlane(fourpoints,area,centroid,plane)
  ! Find the face area and centroids from 4 points
  ! Compute the plane coefficient from 3 points
  !
  ! Author: Bryan He 
  ! Date: 07/19/2022
  !

  implicit none

  PetscReal :: fourpoints(3,4)
  PetscReal :: plane(5)
  PetscReal :: area,centroid(3)  

  call UGridOctCompFaceCentroidArea(fourpoints,area,centroid) 
  
  call GeometryComputePlaneWithPoints1s(plane,fourpoints(:,1), &
                                 fourpoints(:,2),fourpoints(:,3))

end subroutine UGridOctCompFaceCenAreaPlane

! ************************************************************************** !

function GeometryPointInRectangle1(x,y,x_array,y_array)
  !
  ! Determines whether a point in xy space is within
  ! a box
  ! Modified by Bryan He to be exclusive of the boundary values
  ! Data: 07/20/2022
  ! Author: Glenn Hammond
  ! Date: 01/12/12
  !

  implicit none

  PetscReal :: x, y
  PetscReal :: x_array(:), y_array(:)

  PetscBool :: GeometryPointInRectangle1

  PetscReal :: x1, y1
  PetscReal :: x2, y2

  GeometryPointInRectangle1 = PETSC_FALSE

  ! only using first 2 values in each array

  x1 = minval(x_array(:))
  y1 = minval(y_array(:))
  x2 = maxval(x_array(:))
  y2 = maxval(y_array(:))


  if (x > x1 .and. x < x2 .and. y > y1 .and. y < y2) &
    GeometryPointInRectangle1 = PETSC_TRUE

end function GeometryPointInRectangle1

end module Grid_Unstructured_Octree_module
