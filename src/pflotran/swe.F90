module SWE_module

#include "petsc/finclude/petscts.h"
  use petscts

  implicit none

  private

  public :: SWESetup, &
            SWEUpdateAuxVars, &
            SWERHSFunction

contains

! ************************************************************************** !

subroutine SWESetup(surface_realization)
  !
  ! Sets up the variable list for output and observation.
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use Realization_Surface_class
  use Output_Aux_module
  use Patch_module
  use Option_module
  use Grid_module
  use SWE_Aux_module

  implicit none

  class (realization_surface_type) :: surface_realization
  type(output_variable_list_type), pointer :: list

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(swe_auxvar_type), pointer :: swe_auxvars(:)
  PetscInt :: ghosted_id

  list => surface_realization%output_option%output_snap_variable_list
  call SWESetPlotVariables(list)

  list => surface_realization%output_option%output_obs_variable_list
  call SWESetPlotVariables(list)

  option => surface_realization%option
  patch => surface_realization%patch
  patch%surf_aux%SWE => SWEAuxCreate(option)
  grid => patch%grid

  allocate(swe_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call SWEAuxVarInit(swe_auxvars(ghosted_id))
  enddo
  patch%surf_aux%SWE%auxvars => swe_auxvars
  patch%surf_aux%SWE%num_aux = grid%ngmax

end subroutine SWESetup

! ************************************************************************** !

subroutine SWESetPlotVariables(list)
  !
  ! Adds variables to be printed to list
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use Output_Aux_module
  use Variables_module
  use PFLOTRAN_Constants_module

  implicit none

  type(output_variable_list_type), pointer :: list

  character(len=MAXWORDLENGTH) :: name, units

  if (associated(list%first)) then
    return
  endif

  name = 'H'
  units = 'm'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_LIQUID_HEAD)

  name = 'hu'
  units = 'm^2/s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_FLOW_X_MOMENTUM)

  name = 'hv'
  units = 'm^2/s'
  call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                               SURFACE_FLOW_Y_MOMENTUM)

end subroutine SWESetPlotVariables

! ************************************************************************** !

subroutine SWEUpdateAuxVars(surface_realization)
  !
  ! Sets up the variable list for output and observation.
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use Realization_Surface_class
  use Patch_module
  use Option_module
  use Field_Surface_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Surface_Global_Aux_module
  use SWE_Aux_module

  implicit none

  class (realization_surface_type) :: surface_realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_surface_type), pointer :: field_surface
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(swe_auxvar_type), pointer :: swe_auxvars(:)
  type(swe_auxvar_type), pointer :: swe_auxvars_bc(:)
  type(swe_auxvar_type), pointer :: swe_auxvars_ss(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_bc(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars_ss(:)

  PetscInt :: ghosted_id, istart, iend
  PetscReal, pointer :: xx_loc_p(:)
  PetscErrorCode :: ierr

  option => surface_realization%option
  patch => surface_realization%patch
  grid => patch%grid
  field_surface => surface_realization%field_surface

  swe_auxvars => patch%surf_aux%SWE%auxvars
  swe_auxvars_bc => patch%surf_aux%SWE%auxvars_bc
  swe_auxvars_ss => patch%surf_aux%SWE%auxvars_ss

  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars
  surf_global_auxvars_bc => patch%surf_aux%SurfaceGlobal%auxvars_bc
  surf_global_auxvars_ss => patch%surf_aux%SurfaceGlobal%auxvars_ss

  call VecGetArrayF90(field_surface%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !Ignore inactive cells with inactive materials
    if (associated(patch%imat)) then
      if (patch%imat(ghosted_id) <= 0) cycle
    endif

    iend = ghosted_id*option%nflowdof
    istart = iend - option%nflowdof + 1

    call SWEAuxVarCompute(xx_loc_p(istart:iend), swe_auxvars(ghosted_id), &
                          surf_global_auxvars(ghosted_id), option)

  enddo


  call VecRestoreArrayF90(field_surface%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

end subroutine SWEUpdateAuxVars

! ************************************************************************** !

subroutine ComputeRoeFlux(swe_auxvar_up,surf_global_auxvar_up, &
                          swe_auxvar_dn,surf_global_auxvar_dn, &
                          sn,cn,option,flux,amax)
  !
  ! Compute Rimenann flux based on Roe Solver
  !
  ! Author: Gautam Bisht
  ! Date: 05/08/23
  !
  use Realization_Surface_class
  use Option_module
  use Field_Surface_module
  use Field_Surface_module
  use Discretization_module
  use Patch_module
  use Grid_module
  use Connection_module
  use SWE_Aux_module
  use Surface_Global_Aux_module

  implicit none

  type(swe_auxvar_type) ::  swe_auxvar_up, swe_auxvar_dn
  type(surface_global_auxvar_type) :: surf_global_auxvar_up, surf_global_auxvar_dn
  PetscReal :: sn, cn
  type(option_type) :: option
  PetscReal :: flux(option%nflowdof)
  PetscReal :: amax

  PetscReal :: hl,ul,vl
  PetscReal :: hr,ur,vr
  PetscReal :: g
  PetscReal :: duml,dumr
  PetscReal :: cl,cr
  PetscReal :: hhat,uhat,vhat,chat
  PetscReal :: dh,du,dv
  PetscReal :: uperp, dupar, duperp
  PetscReal :: dW(3),R(3,3),A(3,3),FL(3),FR(3)
  PetscReal :: uperpl,uperpr
  PetscReal :: al1,al3,ar1,ar3
  PetscReal :: da1,da3
  PetscReal :: a1,a2,a3

  ! `left = up` and `right = dn`

  hl = surf_global_auxvar_up%h
  ul = swe_auxvar_up%u
  vl = swe_auxvar_up%v

  hr = surf_global_auxvar_dn%h
  ur = swe_auxvar_dn%u
  vr = swe_auxvar_dn%v

  g = abs(option%gravity(3))

  duml = hl**0.5d0
  dumr = hr**0.5d0
  cl = (g*hl)**0.5d0
  cr = (g*hr)**0.5d0

  hhat = duml*dumr
  uhat = (duml*ul + dumr*ur)/(duml+dumr)
  vhat = (duml*vl + dumr*vr)/(duml+dumr)
  chat = (0.5d0*g*(hl+hr))**0.5d0

  uperp = uhat*cn + vhat*sn

  dh = hr - hl
  du = ur - ul
  dv = vr - vl

  dupar = -du*sn + dv*cn
  duperp=  du*cn + dv*sn

  dW(1) = 0.5*(dh - hhat*duperp/chat)
  dW(2) = hhat*dupar
  dW(3) = 0.5*(dh + hhat*duperp/chat)

  uperpl = ul*cn + vl*sn
  uperpr = ur*cn + vr*sn

  al1 = uperpl - cl
  al3 = uperpl + cl
  ar1 = uperpr - cr
  ar3 = uperpr + cr

  R(1,1) = 1.d0
  R(1,2) = 0.d0
  R(1,3) = 1.d0
  R(2,1) = uhat - chat * cn
  R(2,2) = -sn
  R(2,3) = uhat + chat * cn
  R(3,1) = vhat - chat * sn
  R(3,2) = cn
  R(3,3) = vhat + chat * sn

  da1 = max(0.d0, 2.d0 * (ar1 - al1))
  da3 = max(0.d0, 2.d0 * (ar3 - al3))

  a1  = abs(uperp - chat)
  a2  = abs(uperp)
  a3  = abs(uperp + chat)

  ! Critical flow fix
  if (a1 < da1) then
    a1 = 0.5d0 * (a1 * a1 / da1 + da1)
  endif
  if (a3 < da3) then
    a3 = 0.5d0 * (a3 * a3 / da3 + da3)
  endif

  ! Compute interface flux
  A(:,:) = 0.d0

  A(1,1) = a1
  A(2,2) = a2
  A(3,3) = a3

  FL(1) = uperpl * hl;
  FL(2) = ul * uperpl * hl + 0.5d0 * g * hl * hl * cn;
  FL(3) = vl * uperpl * hl + 0.5d0 * g * hl * hl * sn;

  FR(1) = uperpr * hr;
  FR(2) = ur * uperpr * hr + 0.5d0 * g * hr * hr * cn;
  FR(3) = vr * uperpr * hr + 0.5d0 * g * hr * hr * sn;

  ! fflux = 0.5*(FL + FR - matmul(R,matmul(A,dW))
  flux(1) = 0.5d0 * (FL(1) + FR(1) - R(1,1) * A(1,1) * dW(1) - R(1,2) * A(2,2) * dW(2) - R(1,3) * A(3,3) * dW(3));
  flux(2) = 0.5d0 * (FL(2) + FR(2) - R(2,1) * A(1,1) * dW(1) - R(2,2) * A(2,2) * dW(2) - R(2,3) * A(3,3) * dW(3));
  flux(3) = 0.5d0 * (FL(3) + FR(3) - R(3,1) * A(1,1) * dW(1) - R(3,2) * A(2,2) * dW(2) - R(3,3) * A(3,3) * dW(3));

  amax = chat + abs(uperp)

end subroutine ComputeRoeFlux

! ************************************************************************** !

subroutine SWERHSFunctionInternalConn(f,surface_realization,ierr)
  !
  ! Sets up the variable list for output and observation.
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use Realization_Surface_class
  use Option_module
  use Field_Surface_module
  use Field_Surface_module
  use Discretization_module
  use Patch_module
  use Grid_module
  use Connection_module
  use SWE_Aux_module
  use Surface_Global_Aux_module

  implicit none

  Vec :: f
  class (realization_surface_type) :: surface_realization
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(discretization_type), pointer :: discretization
  type(field_surface_type), pointer :: field_surface
  type(option_type), pointer :: option
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(swe_auxvar_type), pointer :: swe_auxvars(:)
  type(surface_global_auxvar_type), pointer :: surf_global_auxvars(:)

  PetscInt :: iconn, sum_connection
  PetscInt :: local_id_up, local_id_dn, face_id
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: vertex_id_up, vertex_id_dn
  PetscReal, pointer :: f_p(:)
  PetscReal :: x_up, y_up, x_dn, y_dn
  PetscReal :: dx, dy, ds
  PetscReal :: sn, cn
  PetscInt, pointer :: face_to_vertex(:,:)
  PetscReal :: flux(surface_realization%option%nflowdof)
  PetscReal :: amax
  PetscReal, pointer :: area_p(:)
  PetscReal :: edge_len,area_up,area_dn
  PetscReal :: cnum, max_courant_num
  PetscInt :: idof, istart

  field_surface => surface_realization%field_surface
  discretization => surface_realization%discretization
  option => surface_realization%option
  patch => surface_realization%patch
  grid => discretization%grid

  swe_auxvars => patch%surf_aux%SWE%auxvars
  surf_global_auxvars => patch%surf_aux%SurfaceGlobal%auxvars

  call VecGetArrayF90(f,f_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field_surface%area,area_p,ierr);CHKERRQ(ierr)

  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  face_to_vertex => grid%unstructured_grid%face_to_vertex
  sum_connection = 0
  max_courant_num = 0.d0

  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      face_id = cur_connection_set%face_id(iconn)

      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)

      vertex_id_up = face_to_vertex(1,face_id)
      vertex_id_dn = face_to_vertex(2,face_id)

      ! If 'h' up and down of the edge are both below the given threshold,
      ! assume both cells are dry and skip them
      if ((surf_global_auxvars(local_id_up)%h < tiny_h .and. &
           surf_global_auxvars(local_id_dn)%h < tiny_h)) cycle



      !
      !                   vertex_dn
      !                       o
      !                      /
      !                     /
      !                    /
      !   ghosted_up -----/-----> ghosted_dn
      !                  /
      !                 /
      !                /
      !               o
      !            vertex_up

      x_up = grid%unstructured_grid%vertices(vertex_id_up)%x
      y_up = grid%unstructured_grid%vertices(vertex_id_up)%y

      x_dn = grid%unstructured_grid%vertices(vertex_id_dn)%x
      y_dn = grid%unstructured_grid%vertices(vertex_id_dn)%y

      dx = x_dn - x_up
      dy = y_dn - y_up
      ds = (dx**2.d0 + dy**2.d0)**0.5d0

      sn = -dx/ds
      cn = dy/ds

      edge_len = cur_connection_set%area(iconn)
      area_up = area_p(ghosted_id_up)
      area_dn = area_p(ghosted_id_dn)

      call ComputeRoeFlux(swe_auxvars(local_id_up),surf_global_auxvars(local_id_up), &
                          swe_auxvars(local_id_dn),surf_global_auxvars(local_id_dn), &
                          sn,cn,option,flux,amax)

      cnum = amax * edge_len / min(area_up, area_dn) * option%flow_dt;
      if (cnum > max_courant_num) then
        max_courant_num = cnum
      endif

      do idof = 1, option%nflowdof

        if (local_id_up > 0) then
          istart = (local_id_up-1)*option%nflowdof + idof
          f_p(istart) = f_p(istart) - flux(idof) * edge_len / area_up;
        endif

        if (local_id_dn > 0) then
          istart = (local_id_dn-1)*option%nflowdof + idof
          f_p(istart) = f_p(istart) + flux(idof) * edge_len / area_dn;
        endif
      enddo
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  call VecRestoreArrayF90(f,f_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field_surface%area,area_p,ierr);CHKERRQ(ierr)

end subroutine SWERHSFunctionInternalConn

! ************************************************************************** !

subroutine SWERHSFunction(ts,time,x,f,surface_realization,ierr)
  !
  ! Sets up the variable list for output and observation.
  !
  ! Author: Gautam Bisht
  ! Date: 04/01/23
  !
  use PFLOTRAN_Constants_module
  use Realization_Surface_class
  use Option_module
  use Field_Surface_module
  use Field_Surface_module
  use Discretization_module

  implicit none

  TS :: ts
  PetscReal :: time
  Vec :: x
  Vec :: f
  class (realization_surface_type) :: surface_realization
  PetscErrorCode :: ierr

  type(discretization_type), pointer :: discretization
  type(field_surface_type), pointer :: field_surface
  type(option_type), pointer :: option

  field_surface => surface_realization%field_surface
  discretization => surface_realization%discretization
  option => surface_realization%option

  write(*,*)'In SWERHSFunction'
  call VecZeroEntries(f,ierr);CHKERRQ(ierr)

  call DiscretizationGlobalToLocal(discretization,x,field_surface%flow_xx_loc,NFLOWDOF)

  call SWEUpdateAuxVars(surface_realization)

  call SWERHSFunctionInternalConn(f,surface_realization,ierr);CHKERRQ(ierr)

  write(*,*)'stopping in SWERHSFunction'
  call exit(0)

end subroutine SWERHSFunction

end module SWE_module
