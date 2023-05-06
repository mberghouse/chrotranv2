module Factory_Surface_Read_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Simulation_Surface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: &
            FactorySurfaceReadFlowPM, &
            FactorySurfaceReadRequiredCards, &
            FactorySurfaceReadInput
contains

! ************************************************************************** !

subroutine FactorySurfaceReadFlowPM(input,option,pm)
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_SWE_class
  use Init_Common_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_base_type), pointer :: pm

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'SIMULATION,PROCESS_MODELS,SURFACE_FLOW'

  nullify(pm)
  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('MODE')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call InputErrorMsg(input,option,'mode',error_string)
        call StringToUpper(word)

        select case(word)
          case('SWE')
            pm => PMSWECreate()
          case default
            error_string = trim(error_string) // ',MODE'
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
        pm%option => option
      case('OPTIONS')
        if (.not.associated(pm)) then
          option%io_buffer = 'MODE keyword must be read first under ' // &
                             trim(error_string)
          call PrintErrMsg(option)
        endif
        call pm%ReadSimulationOptionsBlock(input)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (.not.associated(pm)) then
    option%io_buffer = 'A flow MODE (card) must be included in the &
      &SURFACE_FLOW block in ' // trim(error_string) // '.'
    call PrintErrMsg(option)
  endif

end subroutine FactorySurfaceReadFlowPM

! ************************************************************************** !

subroutine FactorySurfaceReadRequiredCards(simulation,input)
  !
  ! Reads required cards from input file
  !
  !
  ! Author: Gautam Bisht
  ! Date: 03/22/23
  !
  use Option_module
  use Discretization_module
  use Grid_module
  use Input_Aux_module
  use String_module
  use Patch_module
  use Realization_Surface_class
  use HDF5_Aux_module

  use Simulation_Surface_class
  use General_module
  use Reaction_module
  use Reaction_Aux_module
  use NW_Transport_Aux_module
  use Init_Common_module

  implicit none

  class(simulation_surface_type) :: simulation

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  class(realization_surface_type), pointer :: surface_realization
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  type(input_type), pointer :: input
  PetscBool :: found
  PetscBool :: qerr

  character(len = MAXSTRINGLENGTH) :: wname

  surface_realization => simulation%surface_realization
  patch => surface_realization%patch
  option => surface_realization%option
  discretization => surface_realization%discretization

  qerr  = PETSC_FALSE
  wname = '<missing>'
  found = PETSC_FALSE

  call InputPushBlock(input,'SURFACE',option)

  ! GRID information - GRID is a required card for every simulation
  string = "GRID"
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call InputPushBlock(input,'GRID',option)
  call DiscretizationSurfaceReadRequiredCards(discretization,input,option)

  select case(discretization%itype)
    case(STRUCTURED_GRID,UNSTRUCTURED_GRID)
      patch => PatchCreate()
      patch%grid => discretization%grid
      if (.not.associated(surface_realization%patch_list)) then
        surface_realization%patch_list => PatchCreateList()
      endif
      call PatchAddToList(patch,surface_realization%patch_list)
      surface_realization%patch => patch
  end select
  call InputPopBlock(input,option)

  call InputPopBlock(input,option)

end subroutine FactorySurfaceReadRequiredCards

! ************************************************************************** !

subroutine FactorySurfaceReadInput(simulation,input)
  !
  ! Reads pflow input file
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !

  use Option_module
  use Field_Surface_module
  use Grid_module
  use Grid_Unstructured_Aux_module
  use Grid_Structured_module
  use Solver_module
  use Material_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Dataset_module
  use Dataset_Common_HDF5_class
  use Realization_Base_class
  use Region_module
  use Condition_module
  use Coupler_module
  use Strata_module
  use Observation_module
  use Waypoint_module
  use Debug_module
  use Patch_module
  use Discretization_module
  use Input_Aux_module
  use String_module
  use Units_module
  use Regression_module
  use Output_Aux_module
  use Output_module
  use Output_Tecplot_module
  use Data_Mediator_Dataset_class
  use Utility_module
  use Checkpoint_module
  use Simulation_Subsurface_class
  use PMC_Base_class
  use PM_Base_class
  use Print_module
  use Timestepper_Base_class
  use Timestepper_TS_class
  use Time_Storage_module
  use Realization_Surface_class
  use Material_Surface_module
  use Realization_Subsurface_class

  implicit none

  class(simulation_surface_type) :: simulation

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXSTRINGLENGTH) :: string, temp_string
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: error_string

  character(len=1) :: backslash
  PetscReal :: temp_real, temp_real2
  PetscReal, pointer :: temp_real_array(:)
  PetscInt :: temp_int

  PetscBool :: vel_cent
  PetscBool :: vel_face
  PetscBool :: fluxes
  PetscBool :: mass_flowrate
  PetscBool :: energy_flowrate
  PetscBool :: aveg_mass_flowrate
  PetscBool :: aveg_energy_flowrate

  PetscInt :: flag1

  type(region_type), pointer :: region
  type(flow_condition_type), pointer :: flow_condition
  type(coupler_type), pointer :: coupler
  type(strata_type), pointer :: strata
  type(observation_type), pointer :: observation
  class(pmc_base_type), pointer :: master_pmc

  type(waypoint_type), pointer :: waypoint

  type(material_property_type), pointer :: material_property

  class(realization_surface_type), pointer :: realization
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_surface_type), pointer :: field_surface
  type(patch_type), pointer :: patch
  type(output_option_type), pointer :: output_option
  class(dataset_base_type), pointer :: dataset
  class(dataset_ascii_type), pointer :: dataset_ascii
  type(time_storage_type), pointer :: default_time_storage
  class(data_mediator_dataset_type), pointer :: flow_data_mediator
  class(data_mediator_dataset_type), pointer :: rt_data_mediator
  type(waypoint_list_type), pointer :: waypoint_list
  type(waypoint_list_type), pointer :: waypoint_list_time_card
  type(input_type), pointer :: input
  type(material_surface_property_type), pointer :: material_surface_property

  PetscReal :: dt_init
  PetscReal :: dt_min
  PetscReal :: units_conversion

  class(timestepper_base_type), pointer :: temp_timestepper

  PetscReal :: msfsalt, msfwatr, mlfsalt, mlfwatr

  class(pm_base_type), pointer :: pm_flow

  internal_units = 'not_assigned'
  nullify(default_time_storage)
  nullify(waypoint_list_time_card)

  realization => simulation%surface_realization
  output_option => simulation%output_option
  waypoint_list => simulation%waypoint_list_surface
  patch => realization%patch

  if (associated(patch)) grid => patch%grid

  option => realization%option
  field_surface => realization%field_surface

  master_pmc => simulation%surface_flow_process_model_coupler

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++

  call InputRewind(input)
  string = 'SURFACE'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)

  call InputPushBlock(input,'SURFACE',option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    card = trim(word)

    option%io_buffer = 'pflotran card:: ' // trim(card)
    call PrintMsg(option)

    select case(trim(card))

!....................
      case ('GRID')
        call DiscretizationRead(realization%discretization,input,option)

!....................
      case ('SURFACE_MATERIAL_PROPERTY')
        material_surface_property => MaterialSurfacePropertyCreate()

        call InputReadWord(input,option,material_surface_property%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','MATERIAL_PROPERTY')
        call MaterialSurfacePropertyRead(material_surface_property,input,option)
        call MaterialSurfacePropertyAddToList(material_surface_property, &
                                      realization%surf_material_properties)
        nullify(material_surface_property)


!.....................
      case ('TIME')
        dt_init = UNINITIALIZED_DOUBLE
        dt_min = UNINITIALIZED_DOUBLE
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'word','TIME')
          select case(trim(word))
            case('SCREEN_UNITS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Screen Units','TIME')
              internal_units = 'sec'
              temp_real2 = UnitsConvertToInternal(word,internal_units, &
                                                  'TIME,SCREEN_UNITS',option)
              output_option%tunit = trim(word)
              output_option%tconv = temp_real2
            case('STEADY_STATE')
              option%io_buffer = 'STEADY_STATE no longer supported under &
                &TIME card. Please enter under process model OPTIONS.'
              call PrintErrMsg(option)
            case('FINAL_TIME')
              ! cannot use InputReadAndConvertUnits here because we need to
              ! store the units if output_option%tunit is not set
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Final Time',card)
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'Final Time Units',card)
              internal_units = 'sec'
              temp_real2 = UnitsConvertToInternal(word,internal_units, &
                                                  'TIME,FINAL_TIME',option)
              if (len_trim(output_option%tunit) == 0) then
                output_option%tunit = trim(word)
                output_option%tconv = temp_real2
              endif
              waypoint => WaypointCreate()
              waypoint%final = PETSC_TRUE
              waypoint%time = temp_real*temp_real2
              waypoint%print_snap_output = PETSC_TRUE
              ! do not place final time in waypoint_list_time_card
              call WaypointInsertInList(waypoint,waypoint_list)
            case('INITIAL_TIMESTEP_SIZE')
              call InputReadDouble(input,option,dt_init)
              call InputErrorMsg(input,option,'INITIAL_TIMESTEP_SIZE',card)
              internal_units = 'sec'
              call InputReadAndConvertUnits(input,dt_init,internal_units, &
                                            'TIME,INITIAL_TIMESTEP_SIZE', &
                                            option)
            case('MINIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,dt_min)
              call InputErrorMsg(input,option,'MINIMUM_TIMESTEP_SIZE',card)
              internal_units = 'sec'
              call InputReadAndConvertUnits(input,dt_min,internal_units, &
                                            'TIME,MINIMUM_TIMESTEP_SIZE', &
                                            option)
            case('MAXIMUM_TIMESTEP_SIZE')
              call InputReadDouble(input,option,temp_real)
              call InputErrorMsg(input,option,'Minimum Timestep Size',card)
              internal_units = 'sec'
              call InputReadAndConvertUnits(input,temp_real,internal_units, &
                                            'TIME,MAXIMUM_TIMESTEP_SIZE', &
                                            option)
              waypoint => WaypointCreate()
              waypoint%dt_max = temp_real
              call InputReadCard(input,option,word)
              if (input%ierr == 0) then
                call StringToUpper(word)
                if (StringCompare(word,'AT',TWO_INTEGER)) then
                  call InputReadDouble(input,option,waypoint%time)
                  call InputErrorMsg(input,option,'MAXIMUM_TIMESTEP_SIZE &
                                                  &Update Time',card)
                  internal_units = 'sec'
                  call InputReadAndConvertUnits(input,waypoint%time, &
                                                internal_units, &
                                                'TIME,MAXIMUM_TIMESTEP_SIZE,&
                                                &Update Time',option)
                else
                  option%io_buffer = 'Keyword under "MAXIMUM_TIMESTEP_SIZE" &
                                     &after maximum timestep size should &
                                     &be "AT".'
                  call PrintErrMsg(option)
                endif
              else
                waypoint%time = 0.d0
              endif
              if (.not.associated(waypoint_list_time_card)) then
                waypoint_list_time_card => WaypointListCreate()
              endif
              call WaypointInsertInList(waypoint, &
                                        waypoint_list_time_card)
            case default
              call InputKeywordUnrecognized(input,word,'TIME',option)
          end select
        enddo
        call InputPopBlock(input,option)

        ! we store dt_init and dt_min in local variables so that they
        ! cannot overwrite what has previously been set in the respective
        ! timestepper object member variable
        if (Initialized(dt_init)) then
          if (Initialized(master_pmc%timestepper%dt_init)) then
            option%io_buffer = 'INITIAL_TIMESTEP_SIZE may be included &
              &under either the TIME or TIMESTEPPER ' // &
              trim(master_pmc%timestepper%name) // ' card, but not both.'
            call PrintErrMsg(option)
          endif
          if (associated(simulation%surface_flow_process_model_coupler)) then
            temp_timestepper => &
              simulation%surface_flow_process_model_coupler%timestepper
            if (associated(temp_timestepper)) then
              if (Uninitialized(temp_timestepper%dt_init)) then
                temp_timestepper%dt_init = dt_init
              endif
            endif
          endif
        endif
!....................
      case ('REGION')
        region => RegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','REGION')
        call PrintMsg(option,region%name)
        call RegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call RegionAddToList(region,realization%surf_region_list)
        nullify(region)

!....................
      case ('SURFACE_FLOW_CONDITION')
        flow_condition => FlowConditionCreate(option)
        call InputReadWord(input,option,flow_condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'SURF_FLOW_CONDITION','name')
        call PrintMsg(option,flow_condition%name)
        call FlowConditionRead(flow_condition,input,option)
        call FlowConditionAddToList(flow_condition, &
                                    realization%surf_flow_conditions)
        nullify(flow_condition)
!....................
      case ('SURFACE_INITIAL_CONDITION')
        coupler => CouplerCreate(INITIAL_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Initial Condition name')
        call CouplerRead(coupler,input,option)
        call RealizationAddCoupler(realization%patch,coupler)
        nullify(coupler)

        !....................
      case ('STRATIGRAPHY','STRATA')
        strata => StrataCreate()
        call StrataRead(strata,input,option)
        call RealizationAddStrata(realization%patch,strata)
        nullify(strata)

!....................
      case ('END_SURFACE')
        exit

      case default
        call InputKeywordUnrecognized(input,word, &
                                      'SurfaceReadInput()',option)
    end select

  enddo
  call InputPopBlock(input,option) ! SURFACE

  call PrintInitFlags(option%print_flags,option%driver%print_flags)

  ! must come after setup of timestepper steady above. otherwise, the
  ! destruction of the waypoint lists will fail with to pointer to the
  ! same list.
  if (associated(master_pmc%timestepper%local_waypoint_list) .and. &
      associated(waypoint_list_time_card)) then
    option%io_buffer = 'MAXIMUM_TIMESTEP_SIZE may be included under either &
      &the TIME or TIMESTEPPER ' // trim(master_pmc%timestepper%name) // &
      ' card, but not both.'
    call PrintErrMsg(option)
  endif
  if (associated(waypoint_list_time_card)) then
    call WaypointListMerge(simulation%waypoint_list_surface, &
                           waypoint_list_time_card,option)
    ! DO NOT destroy as both pointer point to the same list
    nullify(waypoint_list_time_card)
  else
    call WaypointListMerge(simulation%waypoint_list_surface, &
                           master_pmc%timestepper%local_waypoint_list,option)
  endif

end  subroutine FactorySurfaceReadInput

end module Factory_Surface_Read_module
