module Inversion_Perturbation_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Inversion_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_subsurface_type) :: &
                                            inversion_perturbation_type
    Vec :: quantity_of_interest_base
    Vec :: base_measurement_vec
    PetscInt :: ndof
    PetscInt :: idof_pert
    PetscReal :: pert
    PetscReal :: perturbation_tolerance
    PetscInt, pointer :: select_cells(:)

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

    ! For Wm
    PetscInt :: num_constraints_local    ! Number of constraints
    PetscInt :: num_constraints_total    ! Total number of constraints
    PetscInt, pointer :: rblock(:,:)     ! array stores info about reg.
    PetscReal, pointer :: Wm(:)          ! Regularization matrix

  contains
    procedure, public :: Init => InversionPerturbationInit
    procedure, public :: ReadBlock => InversionPerturbationReadBlock
    procedure, public :: Initialize => InversionPerturbationInitialize
    procedure, public :: Step => InversionPerturbationStep
    procedure, public :: CheckConvergence => &
                           InversionPerturbationCheckConvergence
    procedure, public :: EvaluateCostFunction => &
                           InvPerturbationEvaluateCostFunction
    procedure, public :: ConnectToForwardRun => &
                           InvPerturbationConnectForwardRun
    procedure, public :: CalculateSensitivity => &
                           InvPerturbationCalculateSensitivity
    procedure, public :: WriteIterationInfo => &
                           InversionPerturbationWriteIterationInfo
    procedure, public :: Strip => InversionPerturbationStrip
  end type inversion_perturbation_type

  public :: InversionPerturbationCreate, &
            InversionPerturbationStrip, &
            InversionPerturbationDestroy

contains

! ************************************************************************** !

function InversionPerturbationCreate(driver)
  !
  ! Allocates and initializes a new perturbation inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Driver_module

  class(driver_type), pointer :: driver

  class(inversion_perturbation_type), pointer :: InversionPerturbationCreate

  allocate(InversionPerturbationCreate)
  call InversionPerturbationCreate%Init(driver)

end function InversionPerturbationCreate

! ************************************************************************** !

subroutine InversionPerturbationInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Driver_module
  use ZFlow_Aux_module, only : zflow_calc_adjoint

  class(inversion_perturbation_type) :: this
  class(driver_type), pointer :: driver

  call InversionSubsurfaceInit(this,driver)

  this%quantity_of_interest_base = PETSC_NULL_VEC
  this%base_measurement_vec = PETSC_NULL_VEC

  this%ndof = 0
  this%idof_pert = 0
  this%pert = 0.d0
  this%perturbation_tolerance = 1.d-6
  nullify(this%select_cells)

  ! Default inversion parameters
  this%miniter = 10
  this%maxiter = 50

  this%beta = 100.d0
  this%beta_red_factor = 0.5d0
  this%minperm = 1d-17
  this%maxperm = 1d-07
  this%target_chi2 = 1.d0
  this%min_phi_red = 0.2d0

  this%start_iteration = 1
  this%maximum_iteration = 20
  this%num_constraints_local = UNINITIALIZED_INTEGER
  this%num_constraints_total = UNINITIALIZED_INTEGER
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
  nullify(this%Wm)
  nullify(this%rblock)

  zflow_calc_adjoint = PETSC_FALSE

end subroutine InversionPerturbationInit

! ************************************************************************** !

subroutine InversionPerturbationReadBlock(this,input,option)

  use Input_Aux_module
  use Option_module
  use String_module
  use Utility_module

  class(inversion_perturbation_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'Perturbation Inversion'

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
      case('PERTURBATION_TOLERANCE')
        call InputReadDouble(input,option,this%perturbation_tolerance)
        call InputErrorMsg(input,option,keyword,error_string)
      case('SELECT_CELLS')
        call UtilityReadArray(this%select_cells,ZERO_INTEGER,error_string, &
                              input,option)
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine InversionPerturbationReadBlock

! ************************************************************************** !

subroutine InversionPerturbationInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Inversion_TS_Aux_module
  use String_module

  class(inversion_perturbation_type) :: this

  PetscErrorCode :: ierr

  call InversionSubsurfInitialize(this)
  call InvForwardAuxDestroyList(this%inversion_aux%inversion_forward_aux, &
                                PETSC_FALSE)
  this%inversion_aux%inversion_forward_aux%store_adjoint = PETSC_FALSE

  if (this%idof_pert == 0) then
    if (associated(this%select_cells)) then
      this%ndof = size(this%select_cells)
      if (this%ndof > this%realization%patch%grid%nmax) then
        call this%driver%PrintErrMsg('Number of SELECT_CELLS is larger than &
                                     &the problem size: '// &
                                     trim(StringWrite(this%ndof))//' '// &
                  trim(StringWrite(this%realization%patch%grid%nmax)))
      endif
    else
      this%ndof = this%realization%patch%grid%nmax
    endif
    call VecDuplicate(this%measurement_vec,this%base_measurement_vec, &
                      ierr);CHKERRQ(ierr)
  else
    if (this%idof_pert > this%realization%patch%grid%nmax) then
      call this%driver%PrintErrMsg('SELECT_CELLS ID is larger than &
                                   &the problem size: '// &
                          trim(StringWrite(this%idof_pert))//' '// &
                    trim(StringWrite(this%realization%patch%grid%nmax)))
    endif
  endif

  if (Uninitialized(this%iqoi(1))) then
    call this%driver%PrintErrMsg('Quantity of interest not specified in &
      &InversionPerturbationInitialize.')
  endif

end subroutine InversionPerturbationInitialize

! ************************************************************************** !

subroutine InversionPerturbationStep(this)
  !
  ! Execute a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21

  use Option_module
  use Factory_Forward_module
  use String_module

  class(inversion_perturbation_type) :: this

  type(option_type), pointer :: option
  PetscInt :: iteration

  ! this routine runs the modeling, checks convergence, and
  ! thereafter builds sensitivity matrix using perturbation
  call this%CalculateSensitivity()
  call this%OutputSensitivity('')

  nullify(this%realization)
  call this%forward_simulation%FinalizeRun()
  call this%forward_simulation%Strip()
  deallocate(this%forward_simulation)
  nullify(this%forward_simulation)

  this%converg_flag = PETSC_FALSE
  if (this%iteration > this%maximum_iteration) this%converg_flag = PETSC_TRUE

end subroutine InversionPerturbationStep

! ************************************************************************** !

subroutine InvPerturbationCalculateSensitivity(this)
  !
  ! Calculates sensitivity matrix Jsensitivity
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 03/22/22
  !

  use Option_module
  use Factory_Forward_module
  use String_module

  class(inversion_perturbation_type) :: this

  type(option_type), pointer :: option
  PetscInt :: iteration,ierr

  iteration = 0
  do
    if (associated(this%select_cells) .and. iteration > 0) then
      this%idof_pert = this%select_cells(iteration)
    else
      this%idof_pert = iteration
    endif
    option => OptionCreate()
    option%group_prefix = 'Run' // trim(StringWrite(this%iteration)) // &
                          '_' // StringWrite(iteration)
    call OptionSetDriver(option,this%driver)
    call OptionSetInversionOption(option,this%inversion_option)
    call FactoryForwardInitialize(this%forward_simulation, &
                                  this%forward_simulation_filename,option)
    this%realization => this%forward_simulation%realization
    call this%Initialize()
    call this%forward_simulation%InitializeRun()
    call this%ConnectToForwardRun()
    if (option%status == PROCEED) then
      call this%forward_simulation%ExecuteRun()
    endif
    call InversionPerturbationFillRow(this,iteration)
    if (iteration == 0) then
      ! check convergence and write inversion iteration info
      call this%CheckConvergence()
      call this%WriteIterationInfo()
      !if (this%converg_flag) return
    endif
    if (iteration < this%ndof) then
      nullify(this%realization)
      call this%forward_simulation%FinalizeRun()
      call this%forward_simulation%Strip()
      deallocate(this%forward_simulation)
      nullify(this%forward_simulation)
    endif
    iteration = iteration + 1
    if (iteration > this%ndof) exit
  enddo

end subroutine InvPerturbationCalculateSensitivity

! ************************************************************************** !

subroutine InvPerturbationConnectForwardRun(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Discretization_module
  use Material_module

  class(inversion_perturbation_type) :: this

  type(discretization_type), pointer :: discretization
  Vec :: work
  Vec :: natural_vec
  PetscReal, pointer :: vec_ptr(:)
  PetscReal :: rmin, rmax
  PetscErrorCode :: ierr

  call InvSubsurfConnectToForwardRun(this)

  ! on first pass, store and set thereafter
  if (this%quantity_of_interest_base == PETSC_NULL_VEC) then
    call VecDuplicate(this%quantity_of_interest, &
                      this%quantity_of_interest_base,ierr);CHKERRQ(ierr)
  endif
  if (this%idof_pert == 0) then
    call VecCopy(this%quantity_of_interest,this%quantity_of_interest_base, &
                                     ierr);CHKERRQ(ierr)
  else
    discretization => this%realization%discretization
    work = this%realization%field%work
    call DiscretizationCreateVector(discretization, &
                                    ONEDOF,natural_vec, &
                                    NATURAL,this%realization%option)
    call VecCopy(this%quantity_of_interest_base,this%quantity_of_interest, &
                 ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%quantity_of_interest,vec_ptr, &
                        ierr);CHKERRQ(ierr)
    call VecZeroEntries(natural_vec,ierr);CHKERRQ(ierr)
    if (this%driver%comm%myrank == 0) then
      call VecSetValue(natural_vec,this%idof_pert-1, &
                       this%perturbation_tolerance,INSERT_VALUES, &
                       ierr);CHKERRQ(ierr)
    endif
    call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
    call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                       work,ONEDOF)
    call VecPointwiseMult(work,work,this%quantity_of_interest,ierr);CHKERRQ(ierr)
    call VecMax(work,PETSC_NULL_INTEGER,rmax,ierr);CHKERRQ(ierr)
    call VecMin(work,PETSC_NULL_INTEGER,rmin,ierr);CHKERRQ(ierr)
    if (rmax > 0.d0) then
      this%pert = rmax
    else
      this%pert = rmin
    endif
    call VecAXPY(this%quantity_of_interest,1.d0,work,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(this%realization%discretization, &
                                     this%quantity_of_interest, &
                                     this%realization%field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi(1),this%iqoi(2))
    call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  endif

end subroutine InvPerturbationConnectForwardRun

! ************************************************************************** !

subroutine InversionPerturbationFillRow(this,iteration)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Debug_module
  use Variables_module, only : LIQUID_PRESSURE
  use Realization_Base_class
  use String_module

  class(inversion_perturbation_type) :: this
  PetscInt :: iteration

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscErrorCode :: ierr


  call VecGetArrayF90(this%measurement_vec,vec_ptr,ierr)
  do i = 1, size(this%measurements)
    vec_ptr(i) = this%measurements(i)%simulated_value
  enddo
  call VecRestoreArrayF90(this%measurement_vec,vec_ptr,ierr)

  if (this%idof_pert == 0) then
    call VecCopy(this%measurement_vec,this%base_measurement_vec, &
                 ierr);CHKERRQ(ierr)
    call MatZeroEntries(this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
  else
    call VecAXPY(this%measurement_vec,-1.d0,this%base_measurement_vec, &
                 ierr);CHKERRQ(ierr)
    call VecScale(this%measurement_vec,1.d0/this%pert,ierr);CHKERRQ(ierr)
  endif

  if (this%idof_pert == 0) return

  ! don't need to use the distributed vec, but why not
  call VecScatterBegin(this%scatter_measure_to_dist_measure, &
                       this%measurement_vec,this%dist_measurement_vec, &
                       INSERT_VALUES,SCATTER_FORWARD_LOCAL, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(this%scatter_measure_to_dist_measure, &
                     this%measurement_vec,this%dist_measurement_vec, &
                     INSERT_VALUES,SCATTER_FORWARD_LOCAL, &
                     ierr);CHKERRQ(ierr)
  call VecGetArrayF90(this%dist_measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)
  do i = 1, size(vec_ptr)
    call MatSetValue(this%inversion_aux%JsensitivityT,this%idof_pert-1, &
                     this%dist_measurement_offset+i-1,vec_ptr(i), &
                     INSERT_VALUES,ierr);CHKERRQ(ierr)
  enddo
  call VecRestoreArrayF90(this%dist_measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)

  if (iteration == this%ndof) then
    call MatAssemblyBegin(this%inversion_aux%JsensitivityT,MAT_FINAL_ASSEMBLY, &
                          ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(this%inversion_aux%JsensitivityT,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
    call InversionPerturbationResetToBase(this)
  endif

end subroutine InversionPerturbationFillRow

! ************************************************************************** !

subroutine InversionPerturbationResetToBase(this)
  !
  ! resets all vectors to the base model
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 03/22/22
  !
  use Discretization_module
  use Material_module

  class(inversion_perturbation_type) :: this

  PetscErrorCode :: ierr

  ! reset measurement vectors to the base model
  call VecCopy(this%base_measurement_vec,this%measurement_vec, &
               ierr);CHKERRQ(ierr)
  call VecScatterBegin(this%scatter_measure_to_dist_measure, &
                       this%measurement_vec,this%dist_measurement_vec, &
                       INSERT_VALUES,SCATTER_FORWARD_LOCAL, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(this%scatter_measure_to_dist_measure, &
                     this%measurement_vec,this%dist_measurement_vec, &
                     INSERT_VALUES,SCATTER_FORWARD_LOCAL, &
                     ierr);CHKERRQ(ierr)

  ! reset material property to the base model
  call VecCopy(this%quantity_of_interest_base,this%quantity_of_interest, &
               ierr);CHKERRQ(ierr)
  call DiscretizationGlobalToLocal(this%realization%discretization, &
                                   this%quantity_of_interest, &
                                   this%realization%field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                               this%realization%field%work_loc, &
                               this%iqoi(1),this%iqoi(2))

end subroutine InversionPerturbationResetToBase

! ************************************************************************** !

subroutine InversionPerturbationCheckConvergence(this)
  !
  ! Check Inversion convergence
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 03/22/22

  implicit none

  class(inversion_perturbation_type) :: this

  this%converg_flag = PETSC_FALSE
  call this%EvaluateCostFunction()
  if ((this%current_chi2 <= this%target_chi2) .or. &
      (this%iteration > this%maximum_iteration)) this%converg_flag = PETSC_TRUE

end subroutine InversionPerturbationCheckConvergence

! ************************************************************************** !

subroutine InvPerturbationEvaluateCostFunction(this)
  !
  ! Evaluates cost functions for inversion
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 03/22/22

  implicit none

  class(inversion_perturbation_type) :: this

  PetscInt :: idata,num_measurement
  PetscReal :: wd,tempreal

  num_measurement = size(this%measurements)

  ! Data part
  this%phi_data = 0.d0
  do idata=1,num_measurement

    wd = 0.05 * this%measurements(idata)%value
    wd = 1/wd

    tempreal = wd * (this%measurements(idata)%value - &
                     this%measurements(idata)%simulated_value)
    this%phi_data = this%phi_data + tempreal * tempreal

  enddo

  this%current_chi2 = this%phi_data / num_measurement

end subroutine InvPerturbationEvaluateCostFunction

! ************************************************************************** !

subroutine InversionPerturbationWriteIterationInfo(this)
  !
  ! Writes inversion run info
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 03/22/22
  !

  use String_module

  implicit none

  class(inversion_perturbation_type) :: this

  PetscInt :: fid
  PetscInt, parameter :: zeronum = 0
  character(len=MAXWORDLENGTH) :: string

  if (this%driver%PrintToScreen()) then
    write(*,*)
    write(*,98)
    if (this%iteration == this%start_iteration) then
      write(*,'(/,2x,a,i6.4,/)') StringColor("CONVERGENCE STATISTICS AT &
                              &STARTING ITERATION:",C_RED), this%start_iteration
    else
      write(*,'(/,2x,a,i6.4,/)') StringColor("CONVERGENCE STATISTICS AFTER &
                                 &ITERATION:",C_RED),this%iteration
    endif
    write(*,99)
    write(*,*) StringColor("  Phi_data   ",C_GREEN), &
               StringColor("   Phi_Model  ",C_BLUE), &
               StringColor("  Phi_Model/Beta",C_MAGENTA), &
               StringColor("    Phi_Total   ",C_CYAN)
    write(*,102) this%phi_data,this%phi_model,this%phi_model/this%beta, &
                 this%phi_total
    write(*,*)
    if (this%num_constraints_total >= 0) then
      write(*,103) this%num_constraints_total
    else
      write(*,103) zeronum
    endif
    write(*,104) this%current_chi2
    write(*,105) this%target_chi2
    write(*,106) sqrt(this%current_chi2)
    write(*,107) this%beta
    write(*,108) this%beta_red_factor
    write(*,109) 100.d0*(this%phi_total_0 - this%phi_total)/this%phi_total_0
    write(*,110) 100.d0*this%min_phi_red
    write(*,99)
    write(*,*)
    write(*,98)
  endif

  if (this%driver%PrintToFile()) then
    fid = this%driver%fid_out
    write(fid,*)
    write(fid,98)
    if (this%iteration == this%start_iteration) then
      write(fid,'(/,2x,a,i6.4,/)') "CONVERGENCE STATISTICS AT STARTING &
                                   &ITERATION:", this%start_iteration
    else
      write(fid,'(/,2x,a,i6.4,/)') "CONVERGENCE STATISTICS AFTER ITERATION:", &
                                    this%iteration
    endif
    write(fid,99)
    write(fid,101) "  Phi_data   ","   Phi_Model   "," Phi_Model/Beta", &
                  &"   Phi_Total   "
    write(fid,102) this%phi_data,this%phi_model,this%phi_model/this%beta, &
                   this%phi_total
    write(fid,*)
    if (this%num_constraints_total >= 0) then
      write(fid,103) this%num_constraints_total
    else
      write(fid,103) zeronum
    endif
    write(fid,104) this%current_chi2
    write(fid,105) this%target_chi2
    write(fid,106) sqrt(this%current_chi2)
    write(fid,107) this%beta
    write(fid,108) this%beta_red_factor
    write(fid,109) 100.d0*(this%phi_total_0 - this%phi_total)/this%phi_total_0
    write(fid,110) 100.d0*this%min_phi_red
    write(fid,99)
    write(fid,*)
    write(fid,98)
    flush(fid)
  endif

98 format(40('=+'))
99 format(80('~'))
101 format(4a15)
102 format(4g15.7)

103 format(4x,'Number of Constraint Eqs:      ',2x,i15.10)
104 format(4x,'Current Chi2:                  ',2x,f15.4)
105 format(4x,'Target Ch2:                    ',2x,f15.4)
106 format(4x,'RMS error:                     ',2x,f15.4)
107 format(4x,'Beta:                          ',2x,f15.4)
108 format(4x,'Beta reduction factor:         ',2x,f15.4)
109 format(4x,'Reduction in Phi_Total:        ',2x,f15.4," %")
110 format(4x,'Minimum reduction in Phi_Total ' /,8x, &
                 &'before Beta reduction:     ',2x,f15.4," %")

end subroutine InversionPerturbationWriteIterationInfo

! ************************************************************************** !

subroutine InversionPerturbationStrip(this)
  !
  ! Deallocates members of inversion perturbation
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Utility_module

  class(inversion_perturbation_type) :: this

  PetscErrorCode :: ierr

  call InversionSubsurfaceStrip(this)

  call DeallocateArray(this%select_cells)
  if (this%quantity_of_interest_base /= PETSC_NULL_VEC) then
    call VecDestroy(this%quantity_of_interest_base,ierr);CHKERRQ(ierr)
  endif
  if (this%base_measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%base_measurement_vec,ierr);CHKERRQ(ierr)
  endif

end subroutine InversionPerturbationStrip

! ************************************************************************** !

subroutine InversionPerturbationDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  class(inversion_perturbation_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionPerturbationDestroy

end module Inversion_Perturbation_class
