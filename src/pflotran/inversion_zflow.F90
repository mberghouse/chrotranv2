module Inversion_ZFlow_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Inversion_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_subsurface_type) :: inversion_zflow_type
    PetscInt :: start_iteration          ! Starting iteration number
    PetscInt :: miniter,maxiter          ! min/max CGLS iterations

    PetscReal :: beta                    ! regularization parameter
    PetscReal :: beta_red_factor         ! beta reduction factor
    PetscReal :: minperm,maxperm         ! min/max permeability
    PetscReal :: target_chi2             ! target CHI^2 norm
    PetscReal :: current_chi2

    ! Cost/objective functions
    PetscReal :: min_phi_red             ! min change in cost function
    PetscReal :: phi_total_0,phi_total
    PetscReal :: phi_data_0,phi_data
    PetscReal :: phi_model_0,phi_model

    ! arrays for CGLS algorithm
    PetscReal, pointer :: b(:)           ! vector for CGLS RHS
    PetscReal, pointer :: p(:)           ! vector of dim -> num of inv cells
    PetscReal, pointer :: q(:)           ! product of Jacobian with p = Jp
    PetscReal, pointer :: r(:)           ! vector of dim -> num of measur
    PetscReal, pointer :: s(:)           ! product Jacobian transpose with r
    PetscReal, pointer :: del_perm(:)    ! permeability update vector

  contains
    procedure, public :: Init => InversionZFlowInit
    procedure, public :: ReadBlock => InversionZFlowReadBlock
    procedure, public :: Initialize => InversionZFlowInitialize
    procedure, public :: Step => InversionZFlowStep
    procedure, public :: UpdateParameters => InversionZFlowUpdateParameters
    procedure, public :: CalculateUpdate => InversionZFlowCalculateUpdate
    procedure, public :: CheckConvergence => InversionZFlowCheckConvergence
    procedure, public :: EvaluateCostFunction => InvZFlowEvaluateCostFunction
    procedure, public :: Finalize => InversionZFlowFinalize
    procedure, public :: Strip => InversionZFlowStrip
  end type inversion_zflow_type

  public :: InversionZFlowCreate, &
            InversionZFlowFinalize, &
            InversionZFlowStrip, &
            InversionZFlowDestroy

contains

! ************************************************************************** !

function InversionZFlowCreate(driver)
  !
  ! Allocates and initializes a new zflow inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !

  use Driver_module

  class(driver_type), pointer :: driver

  class(inversion_zflow_type), pointer :: InversionZFlowCreate

  allocate(InversionZFlowCreate)
  call InversionZFlowCreate%Init(driver)

end function InversionZFlowCreate

! ************************************************************************** !

subroutine InversionZFlowInit(this,driver)
  !
  ! Initializes a new zflow inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !
  use Variables_module, only : PERMEABILITY
  use Driver_module

  implicit none

  class(inversion_zflow_type) :: this
  class(driver_type), pointer :: driver

  call InversionSubsurfaceInit(this,driver)
  ! override default set in InversionSubsurfaceInit
  this%iqoi = PERMEABILITY

  ! Default inversion parameters
  this%miniter = 10
  this%maxiter = 50

  this%beta = 100.d0
  this%beta_red_factor = 0.5d0
  this%minperm = 1d-25
  this%maxperm = 1d-02
  this%target_chi2 = 1.d0
  this%min_phi_red = 0.2d0

  this%start_iteration = 1
  this%maximum_iteration = 20

  this%current_chi2 = UNINITIALIZED_DOUBLE
  this%phi_total_0 = UNINITIALIZED_DOUBLE
  this%phi_data_0 = UNINITIALIZED_DOUBLE
  this%phi_model_0 = UNINITIALIZED_DOUBLE
  this%phi_total = UNINITIALIZED_DOUBLE
  this%phi_data = UNINITIALIZED_DOUBLE
  this%phi_model = UNINITIALIZED_DOUBLE

  nullify(this%b)
  nullify(this%p)
  nullify(this%q)
  nullify(this%r)
  nullify(this%s)
  nullify(this%del_perm)

end subroutine InversionZFlowInit

! ************************************************************************** !

subroutine InversionZFlowAllocateWorkArrays(this)
  !
  ! Initialize inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/13/21
  !

  use Grid_module

  implicit none

  class(inversion_zflow_type) :: this

  type(grid_type), pointer :: grid

  PetscInt :: num_measurement
  PetscInt :: num_constraints

  grid => this%realization%patch%grid

  num_measurement = size(this%measurement)
  num_constraints = 0

  allocate(this%b(num_measurement + num_constraints))
  allocate(this%p(grid%nlmax))
  allocate(this%q(num_measurement + num_constraints))
  allocate(this%r(num_measurement + num_constraints))
  allocate(this%s(grid%nlmax))
  allocate(this%del_perm(grid%nlmax))

  this%b = 0.d0
  this%p = 0.d0
  this%q = 0.d0
  this%r = 0.d0
  this%s = 0.d0
  this%del_perm = 0.d0

end subroutine InversionZFlowAllocateWorkArrays

! ************************************************************************** !

subroutine InversionZFlowDeallocateWorkArrays(this)
  !
  ! Initialize inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/13/21
  !

  use Utility_module, only : DeallocateArray

  implicit none

  class(inversion_zflow_type) :: this

  call DeallocateArray(this%b)
  call DeallocateArray(this%p)
  call DeallocateArray(this%q)
  call DeallocateArray(this%r)
  call DeallocateArray(this%s)
  call DeallocateArray(this%del_perm)

end subroutine InversionZFlowDeallocateWorkArrays

! ************************************************************************** !

! ************************************************************************** !

subroutine InversionZFlowReadBlock(this,input,option)
  !
  ! Reads input file parameters associated an ZFlow inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !

  use Input_Aux_module
  use Option_module
  use String_module

  class(inversion_zflow_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'ZFlow Inversion'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_TRUE
    call InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))

      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine InversionZFlowReadBlock

! ************************************************************************** !

subroutine InversionZFlowInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !
  use Discretization_module
  use Material_module
  use Option_module
  use Variables_module, only : PERMEABILITY

  implicit none

  class(inversion_zflow_type) :: this

  PetscBool :: exists
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  ! theck to ensure that quantity of interest exists
  exists = PETSC_FALSE
  select case(this%iqoi)
    case(PERMEABILITY)
      if (this%realization%option%iflowmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'PERMEABILITY'
    case default
  end select
  if (.not.exists) then
    this%realization%option%io_buffer = 'Inversion for ' // trim(word) // &
      &' cannot be performed with the specified process models.'
    call PrintErrMsg(this%realization%option)
  endif

  call InversionBaseInitialize(this)

end subroutine InversionZFlowInitialize

! ************************************************************************** !

subroutine InversionZFlowSetup(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hamond
  ! Date: 10/11/21
  !
  use Discretization_module
  use Material_module
  use PM_ZFlow_class
  use Inversion_Aux_module

  class(inversion_zflow_type) :: this

  PetscErrorCode :: ierr

  call InversionSubsurfaceSetup(this)

  call DiscretizationGlobalToLocal(this%realization%discretization, &
                                   this%quantity_of_interest, &
                                   this%realization%field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                               this%realization%field%work_loc, &
                               this%iqoi,ZERO_INTEGER)

  ! this code must follow InversionSubsurfaceSetup() as Jsensitivity is
  ! is created within it.
  select type(pm => &
                this%forward_simulation%flow_process_model_coupler%pm_list)
    class is(pm_zflow_type)
      pm%inversion_aux => InversionAuxCreate()
      pm%inversion_aux%Jsensitivity = this%Jsensitivity
      pm%inversion_aux%imeasurement => this%imeasurement
      pm%inversion_aux%measurement => this%measurement
    class default
      call this%driver%PrintErrMsg('Unsupported process model in &
                                   &InversionZFlowSetup.')
  end select

  ! Build Wm matrix
  !call InversionZFlowBuildWm(this)

end subroutine InversionZFlowSetup

! ************************************************************************** !

subroutine InversionZFlowStep(this)
  !
  ! Execute a simulation
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21

  use Option_module
  use Factory_Forward_module

  class(inversion_zflow_type) :: this

  type(option_type), pointer :: option

  option => OptionCreate()
  write(option%group_prefix,'(i6)') this%iteration
  option%group_prefix = 'Run' // trim(adjustl(option%group_prefix))
  call OptionSetDriver(option,this%driver)
  call FactoryForwardInitialize(this%forward_simulation, &
                                this%forward_simulation_filename,option)
  this%realization => this%forward_simulation%realization
  call this%Initialize()
  call this%forward_simulation%InitializeRun()
  call InversionZFlowSetup(this)
  if (option%status == PROCEED) then
    call this%forward_simulation%ExecuteRun()
  endif
  call this%CheckConvergence()
print*,"PHI DATA: ", this%phi_data
  call this%CalculateUpdate()
  call this%UpdateParameters()
  call this%forward_simulation%FinalizeRun()
  call this%forward_simulation%Strip()
  deallocate(this%forward_simulation)
  nullify(this%forward_simulation)

  this%converg_flag = PETSC_FALSE
  if (this%iteration > this%maximum_iteration) this%converg_flag = PETSC_TRUE

end subroutine InversionZFlowStep

! ************************************************************************** !

subroutine InversionZFlowCheckConvergence(this)
  !
  ! Check Inversion convergence
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/15/21

  implicit none

  class(inversion_zflow_type) :: this

  this%converg_flag = PETSC_FALSE
  call this%EvaluateCostFunction()
  if ((this%current_chi2 <= this%target_chi2) .or. &
      (this%iteration > this%maximum_iteration)) this%converg_flag = PETSC_TRUE

end subroutine InversionZFlowCheckConvergence

! ************************************************************************** !

subroutine InvZFlowEvaluateCostFunction(this)
  !
  ! Evaluates cost functions for inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/15/21

  use Variables_module, only : LIQUID_PRESSURE
  use Realization_Base_class
  use Option_module

  implicit none

  class(inversion_zflow_type) :: this

  type(option_type), pointer :: option

  PetscInt :: icell
  PetscInt :: idata,num_measurement
  PetscReal :: wd,tempreal
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  num_measurement = size(this%measurement)

  call RealizationGetVariable(this%realization, &
                              this%realization%field%work, &
                              LIQUID_PRESSURE,ZERO_INTEGER)

  call VecGetArrayF90(this%realization%field%work,vec_ptr,ierr);CHKERRQ(ierr)

  ! Data part
  this%phi_data = 0.d0
  do idata=1,num_measurement

    wd = 0.05 * this%measurement(idata)
    wd = 1/wd

    icell = this%imeasurement(idata)
    tempreal = wd * (this%measurement(idata) - vec_ptr(icell))
    this%phi_data = this%phi_data + tempreal * tempreal
!print*,icell,vec_ptr(icell)
  enddo

  this%current_chi2 = this%phi_data

  call VecRestoreArrayF90(this%realization%field%work,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine InvZFlowEvaluateCostFunction

! ************************************************************************** !

subroutine InversionZFlowUpdateParameters(this)
  !
  ! Updates input parameters
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/14/21
  !

  use Material_module
  use Discretization_module
  use Field_module

  class(inversion_zflow_type) :: this

  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization

  PetscInt :: local_id,ghosted_id
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  field => this%realization%field
  discretization => this%realization%discretization

  call DiscretizationGlobalToLocal(discretization, &
                                   this%quantity_of_interest, &
                                   field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                               field%work_loc,this%iqoi,ZERO_INTEGER)

end subroutine InversionZFlowUpdateParameters

! ************************************************************************** !

subroutine InversionZFlowCalculateUpdate(this)
  !
  ! Calculates updated model parameters
  ! using m_new = m_old + del_m
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/13/21
  !

  use Patch_module
  use Grid_module

  implicit none

  class(inversion_zflow_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid

  PetscInt :: local_id
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  patch => this%realization%patch
  grid => patch%grid

  if (this%quantity_of_interest /= PETSC_NULL_VEC) then

    call InversionZFlowAllocateWorkArrays(this)

    ! get inversion%del_perm
    call InversionZFlowCGLSSolve(this)

    ! Get updated conductivity as m_new = m_old + del_m (where m = log(perm))
    call VecGetArrayF90(this%quantity_of_interest,vec_ptr,ierr);CHKERRQ(ierr)
    do local_id=1,grid%nlmax
     vec_ptr(local_id) = exp(log(vec_ptr(local_id)) + this%del_perm(local_id))
!print*,local_id,exp(this%del_perm(local_id))
     if (vec_ptr(local_id) > this%maxperm) vec_ptr(local_id) = this%maxperm
     if (vec_ptr(local_id) < this%minperm) vec_ptr(local_id) = this%minperm
    enddo
    call VecRestoreArrayF90(this%quantity_of_interest,vec_ptr, &
                                                          ierr);CHKERRQ(ierr)
    call InversionZFlowDeallocateWorkArrays(this)

  endif

end subroutine InversionZFlowCalculateUpdate

! ************************************************************************** !

subroutine InversionZFlowCGLSSolve(this)
  !
  ! Implements CGLS solver for least sqaure equivalent
  !            of the ZFlow normal equations
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/14/21
  !

  use Option_module
  use Timer_class
  use String_module

  implicit none

  class(inversion_zflow_type) :: this

  type(option_type), pointer :: option
  class(timer_type), pointer :: timer

  PetscInt :: i,nm,ncons
  PetscReal :: alpha,gbeta,gamma,gamma1,delta1,delta2,delta
  PetscReal :: norms0,norms,normx,xmax
  PetscReal :: resNE,resNE_old
  PetscBool :: exit_info,indefinite
  PetscErrorCode :: ierr

  PetscReal, parameter :: delta_initer = 1e-23
  PetscReal, parameter :: initer_conv  = 1e-24

  option => this%realization%option

  this%del_perm = 0.0d0

  timer => TimerCreate()
  call timer%Start()

  if (OptionPrintToScreen(option)) then
    write(*,'(" --> Solving ZFlow normal equation using CGLS solver:")')
  endif

  nm = size(this%measurement)
  ncons = 0

  ! Get RHS vector this%b
  call InversionZFlowCGLSRhs(this)

  this%r = this%b

  ! get this%s = J^tr
  call InversionZFlowComputeMatVecProductJtr(this)
  this%p = this%s

  gamma = dot_product(this%s,this%s)
  call MPI_Allreduce(MPI_IN_PLACE,gamma,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

  norms0 = sqrt(gamma)
  xmax = 0.d0
  normx = 0.d0
  resNE = 0.d0
  exit_info = PETSC_FALSE
  indefinite = PETSC_FALSE

  do i=1,this%maxiter

    if (exit_info) exit

    ! get this%q = Jp
    call InversionZFlowComputeMatVecProductJp(this)

    delta1 = dot_product(this%q(1:nm),this%q(1:nm))
    delta2 = dot_product(this%q(nm+1:nm+ncons),this%q(nm+1:nm+ncons))
    call MPI_Allreduce(MPI_IN_PLACE,delta2,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
    delta = delta1 + delta2

    if (delta < 0) indefinite = PETSC_TRUE
    if (delta == 0) delta = epsilon(delta)

    alpha = gamma / delta

    this%del_perm = this%del_perm + alpha * this%p
    this%r = this%r - alpha * this%q

    ! get this%s = J^tr
    call InversionZFlowComputeMatVecProductJtr(this)

    gamma1 = gamma
    gamma = dot_product(this%s,this%s)
    call MPI_Allreduce(MPI_IN_PLACE,gamma,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

    norms = sqrt(gamma)
    gbeta = gamma / gamma1
    this%p = this%s + gbeta * this%p

!print*,norm2(this%r),gamma

    normx = dot_product(this%del_perm,this%del_perm)
    call MPI_Allreduce(MPI_IN_PLACE,normx,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
    normx = sqrt(normx)
    if (xmax < normx) xmax = normx
    if ( (norms <= norms0 * initer_conv) .or. (normx * initer_conv >= 1)) &
                               exit_info = PETSC_TRUE

    resNE_old = resNE
    resNE = norms / norms0

    if ( abs((resNE_old - resNe) /resNE_old) < delta_initer .and. &
        i > this%miniter) exit_info = PETSC_TRUE

  enddo

!print*,this%del_perm
!stop

  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime())) &
    // ' seconds and ' // trim(StringWrite(i)) // &
    ' iterations to solve ZFlow normal equation.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InversionZFlowCGLSSolve

! ************************************************************************** !

subroutine InversionZFlowCGLSRhs(this)
  !
  ! Builds RHS for least-square equation for CGLS solver
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/14/21
  !

  use Variables_module, only : LIQUID_PRESSURE
  use Realization_Base_class
  use Patch_module
  use Material_Aux_class
  use Option_module

  implicit none

  class(inversion_zflow_type) :: this

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch

  PetscInt :: icell
  PetscInt :: idata,num_measurement
  PetscReal :: wd
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  option => this%realization%option
  patch => this%realization%patch

  this%b = 0.0d0

  num_measurement = size(this%measurement)

  call RealizationGetVariable(this%realization, &
                              this%realization%field%work, &
                              LIQUID_PRESSURE,ZERO_INTEGER)
  call VecGetArrayF90(this%realization%field%work,vec_ptr,ierr);CHKERRQ(ierr)

  ! Data part
  do idata=1,num_measurement

    wd = 0.05 * this%measurement(idata)
    wd = 1/wd

    icell = this%imeasurement(idata)
    this%b(idata) = wd * (this%measurement(idata) - vec_ptr(icell))
  enddo

  call VecRestoreArrayF90(this%realization%field%work,vec_ptr,ierr);CHKERRQ(ierr)

  ! Model part

end subroutine InversionZFlowCGLSRhs

! ************************************************************************** !

subroutine InversionZFlowComputeMatVecProductJp(this)
  !
  ! Computes product of Jacobian J with a vector p = Jp
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/14/21
  !

  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module
  use Option_module

  implicit none

  class(inversion_zflow_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization

  PetscInt :: num_measurement
  PetscReal, pointer :: pvec_ptr(:)
  PetscReal, pointer :: qvec_ptr(:)
  PetscErrorCode :: ierr

  Vec :: p1
  Vec :: q1

  option => this%realization%option
  field => this%realization%field
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid

  this%q = 0.d0

  num_measurement = size(this%measurement)

  ! Data part
  call VecDuplicate(field%work,p1,ierr);CHKERRQ(ierr)
  call VecCreateMPI(this%driver%comm%mycomm,num_measurement,num_measurement, &
                    q1,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(p1,pvec_ptr,ierr);CHKERRQ(ierr)
  pvec_ptr = this%p
  call VecRestoreArrayF90(p1,pvec_ptr,ierr);CHKERRQ(ierr)

  ! q = Jp -> data part
  call MatMult(this%Jsensitivity,p1,q1,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(q1,qvec_ptr,ierr);CHKERRQ(ierr)
  this%q = qvec_ptr
  call VecRestoreArrayF90(q1,qvec_ptr,ierr);CHKERRQ(ierr)

  ! Model part
  ! Get local this%p to ghosted in pvec_ptr
  !call DiscretizationGlobalToLocal(discretization,p1, &
  !                                 field%work_loc,ONEDOF)
  !call VecGetArrayF90(field%work_loc,pvec_ptr,ierr);CHKERRQ(ierr)


  !call VecRestoreArrayF90(field%work_loc,pvec_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(p1,ierr);CHKERRQ(ierr)
  call VecDestroy(q1,ierr);CHKERRQ(ierr)

end subroutine InversionZFlowComputeMatVecProductJp

!************************************************************************** !

subroutine InversionZFlowComputeMatVecProductJtr(this)
  !
  ! Computes product of Jacobian J transpose with a vector r = J^t x r
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/14/21
  !

  use Patch_module
  use Grid_module
  use Field_module
  use Discretization_module

  implicit none

  class(inversion_zflow_type) :: this

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization

  PetscInt :: num_measurement
  PetscReal, pointer :: rvec_ptr(:)
  PetscReal, pointer :: svec_ptr(:)
  PetscErrorCode :: ierr

  Vec :: r1
  Vec :: s1

  field => this%realization%field
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid

  this%s = 0.0d0

  num_measurement = size(this%measurement)

  ! Data part
  call VecCreateMPI(this%driver%comm%mycomm,num_measurement,num_measurement, &
                    r1,ierr);CHKERRQ(ierr)
  call VecDuplicate(field%work,s1,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(r1,rvec_ptr,ierr);CHKERRQ(ierr)
  rvec_ptr = this%r
  call VecRestoreArrayF90(r1,rvec_ptr,ierr);CHKERRQ(ierr)

  ! s = J^T*r -> data part
  call MatMultTranspose(this%Jsensitivity,r1,s1,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(s1,svec_ptr,ierr);CHKERRQ(ierr)
  this%s = svec_ptr
  call VecRestoreArrayF90(s1,svec_ptr,ierr);CHKERRQ(ierr)

end subroutine InversionZFlowComputeMatVecProductJtr

! ************************************************************************** !

subroutine InversionZFlowFinalize(this)
  !
  ! Finalizes inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !

  implicit none

  class(inversion_zflow_type) :: this

  call InversionBaseStrip(this)

end subroutine InversionZFlowFinalize

! ************************************************************************** !

subroutine InversionZFlowStrip(this)
  !
  ! Deallocates members of inversion zflow
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !

  implicit none

  class(inversion_zflow_type) :: this

   call InversionSubsurfaceStrip(this)

end subroutine InversionZFlowStrip

! ************************************************************************** !

subroutine InversionZFlowDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 10/11/21
  !

  use Utility_module, only : DeallocateArray

  implicit none

  class(inversion_zflow_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call DeallocateArray(inversion%b)
  call DeallocateArray(inversion%p)
  call DeallocateArray(inversion%q)
  call DeallocateArray(inversion%r)
  call DeallocateArray(inversion%s)
  call DeallocateArray(inversion%del_perm)

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionZFlowDestroy

end module Inversion_ZFlow_class