module Hydrate_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module

  implicit none
  
  private 

  PetscBool, public :: hydrate_print_state_transition = PETSC_TRUE
  PetscBool, public :: hydrate_analytical_derivatives = PETSC_FALSE
  PetscBool, public :: hydrate_immiscible = PETSC_FALSE
  PetscBool, public :: hydrate_central_diff_jacobian = PETSC_FALSE
  PetscBool, public :: hydrate_restrict_state_chng = PETSC_FALSE
  PetscReal, public :: window_epsilon = 1.d-4 !0.d0
  PetscReal, public :: hydrate_phase_chng_epsilon = 0.d0 !1.d-6
  PetscReal, public :: fmw_comp(3) = [FMWH2O,FMWAIR,FMWCO2]
  PetscReal, public :: hydrate_max_pressure_change = 5.d4
  PetscInt, public :: hydrate_max_it_before_damping = UNINITIALIZED_INTEGER
  PetscReal, public :: hydrate_damping_factor = 0.6d0
  PetscInt, public :: hydrate_debug_cell_id = UNINITIALIZED_INTEGER
  PetscBool, public :: hydrate_diffuse_xmol = PETSC_TRUE
  PetscBool, public :: hydrate_temp_dep_gas_air_diff = PETSC_TRUE
  PetscBool, public :: hydrate_harmonic_diff_density = PETSC_TRUE
  PetscInt, public :: hydrate_newton_iteration_number = 0
  PetscBool, public :: hydrate_allow_state_change = PETSC_TRUE
  PetscBool, public :: hydrate_force_iteration = PETSC_FALSE
  PetscBool, public :: hydrate_state_changed = PETSC_FALSE

  PetscBool, public :: hyd_chk_max_dpl_liq_state_only = PETSC_FALSE
  PetscBool, public :: hydrate_high_temp_ts_cut = PETSC_FALSE

  ! debugging
  PetscInt, public :: hydrate_ni_count
  PetscInt, public :: hydrate_ts_cut_count
  PetscInt, public :: hydrate_ts_count

  ! thermodynamic state of fluid ids
  !PetscInt, parameter, public :: NULL_STATE = 0
  !PetscInt, parameter, public :: L_STATE = 1
  !PetscInt, parameter, public :: G_STATE = 2
  !PetscInt, parameter, public :: GA_STATE = 3
  !PetscInt, parameter, public :: ANY_STATE = 4
  !PetscInt, parameter, public :: MULTI_STATE = 5
  
  PetscInt, parameter, public :: PREV_TS = 1
  PetscInt, parameter, public :: PREV_IT = 2

  PetscInt, parameter, public :: HYDRATE_LIQUID_PRESSURE_DOF = 1
  PetscInt, parameter, public :: HYDRATE_GAS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: HYDRATE_2PH_STATE_AIR_PRESSURE_DOF = 3
  PetscInt, parameter, public :: HYDRATE_G_STATE_AIR_PRESSURE_DOF = 2
  PetscInt, parameter, public :: HYDRATE_GAS_SATURATION_DOF = 2
  
  PetscInt, parameter, public :: HYDRATE_ENERGY_DOF = 3
  PetscInt, parameter, public :: HYDRATE_L_STATE_X_MOLE_DOF = 2
  PetscInt, parameter, public :: HYDRATE_L_STATE_X2_MOLE_DOF = 4
  PetscInt, parameter, public :: HYDRATE_HYD2_SATURATION_DOF = 4
  PetscInt, parameter, public :: HYDRATE_G_STATE_CO2_PRESSURE_DOF = 4  

  PetscInt, parameter, public :: HYDRATE_STATE_INDEX = 1
  PetscInt, parameter, public :: HYDRATE_LIQUID_EQUATION_INDEX = 1
  PetscInt, parameter, public :: HYDRATE_GAS_EQUATION_INDEX = 2
  PetscInt, parameter, public :: HYDRATE_ENERGY_EQUATION_INDEX = 3
  PetscInt, parameter, public :: HYDRATE_CO2_EQUATION_INDEX = 4
  
  PetscInt, parameter, public :: HYDRATE_LIQUID_PRESSURE_INDEX = 2
  PetscInt, parameter, public :: HYDRATE_GAS_PRESSURE_INDEX = 3
  PetscInt, parameter, public :: HYDRATE_AIR_PRESSURE_INDEX = 4
  PetscInt, parameter, public :: HYDRATE_LIQ_MOLE_FRACTION_INDEX = 5
  PetscInt, parameter, public :: HYDRATE_HYD_MOLE_FRACTION_INDEX = 6
  PetscInt, parameter, public :: HYDRATE_ICE_MOLE_FRACTION_INDEX = 7
  PetscInt, parameter, public :: HYDRATE_TEMPERATURE_INDEX = 8
  PetscInt, parameter, public :: HYDRATE_GAS_SATURATION_INDEX = 9
  PetscInt, parameter, public :: HYDRATE_LIQ_SATURATION_INDEX = 10
  PetscInt, parameter, public :: HYDRATE_HYD_SATURATION_INDEX = 11
  PetscInt, parameter, public :: HYDRATE_ICE_SATURATION_INDEX = 12
  PetscInt, parameter, public :: HYDRATE_LIQUID_FLUX_INDEX = 13
  PetscInt, parameter, public :: HYDRATE_GAS_FLUX_INDEX = 14
  PetscInt, parameter, public :: HYDRATE_ENERGY_FLUX_INDEX = 15
  PetscInt, parameter, public :: HYDRATE_LIQUID_CONDUCTANCE_INDEX = 16
  PetscInt, parameter, public :: HYDRATE_GAS_CONDUCTANCE_INDEX = 17
  PetscInt, parameter, public :: HYDRATE_HYD2_SATURATION_INDEX = 18
  PetscInt, parameter, public :: HYDRATE_HYD2_MOLE_FRACTION_INDEX = 19
  PetscInt, parameter, public :: HYDRATE_LIQ_MOLE2_FRACTION_INDEX = 20
  PetscInt, parameter, public :: HYDRATE_CO2_PRESSURE_INDEX = 21
  PetscInt, parameter, public :: HYDRATE_GAS2_FLUX_INDEX = 22
  PetscInt, parameter, public :: HYDRATE_MAX_INDEX = 23
  
  PetscInt, parameter, public :: HYDRATE_UPDATE_FOR_DERIVATIVE = -1
  PetscInt, parameter, public :: HYDRATE_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: HYDRATE_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: HYDRATE_UPDATE_FOR_BOUNDARY = 2
  PetscInt, parameter, public :: HYDRATE_UPDATE_FOR_SS = 3

  PetscReal, parameter, public :: HYDRATE_IMMISCIBLE_VALUE = 1.d-10
  PetscReal, parameter, public :: HYDRATE_PRESSURE_SCALE = 1.d0
  
  ! these variables, which are global to hydrate, can be modified
  PetscInt, public, pointer :: dof_to_primary_variable(:,:)
  PetscInt, public :: hydrate_2ph_energy_dof = HYDRATE_TEMPERATURE_INDEX

  PetscInt, parameter, public :: HYD_MULTI_STATE = -2 
  PetscInt, parameter, public :: HYD_ANY_STATE = -1 
  PetscInt, parameter, public :: NULL_STATE = 0
  PetscInt, parameter, public :: L_STATE = 1
  PetscInt, parameter, public :: G_STATE = 2
  PetscInt, parameter, public :: H_STATE = 3 !5 (4 and 5 conflict with 
  PetscInt, parameter, public :: I_STATE = 4 !4 ANY_STATE and MULTI_STATE)
  PetscInt, parameter, public :: GA_STATE = 5
  PetscInt, parameter, public :: HG_STATE = 6
  PetscInt, parameter, public :: HA_STATE = 7
  PetscInt, parameter, public :: HI_STATE = 8
  PetscInt, parameter, public :: GI_STATE = 9
  PetscInt, parameter, public :: AI_STATE = 10
  PetscInt, parameter, public :: HGA_STATE = 11
  PetscInt, parameter, public :: HAI_STATE = 12
  PetscInt, parameter, public :: HGI_STATE = 13
  PetscInt, parameter, public :: GAI_STATE = 14
  PetscInt, parameter, public :: QUAD_STATE = 15
  PetscInt, parameter, public :: H2_STATE = 16
  PetscInt, parameter, public :: H2A_STATE = 17
  PetscInt, parameter, public :: H2G_STATE = 18
  PetscInt, parameter, public :: HHA_STATE = 19
  PetscInt, parameter, public :: HHGA_STATE = 20
  PetscInt, parameter, public :: HH_STATE = 21


  PetscInt, parameter :: lid = 1
  PetscInt, parameter :: gid = 2
  PetscInt, parameter :: hid = 3
  PetscInt, parameter :: iid = 4
  PetscInt, parameter :: h2id = 5

  !Structure 1 methane hydrate:
  PetscReal, public :: Nhyd = 6.d0
  PetscReal, public :: HYDRATE_DENSITY_KG = 920.d0 !kg/m^3
  PetscReal, public :: HYDRATE_DENSITY = 52.15551276d0 !mol/L
  PetscReal, public :: MW_CH4 = 16.04d0
  PetscReal, public :: MW_H2O = 18.01d0
  PetscReal, public :: MW_CO2 = 44.0098d0

  PetscReal, parameter, public :: MOL_RATIO_METH = 0.14285714285d0
  PetscReal, public :: MOL_RATIO_H2O = 1.d0 - MOL_RATIO_METH
  PetscReal, public :: MOL_RATIO_FORMER = 0.14285714285d0

  PetscReal, public :: TQD = 1.d-2 !0.d0 !1.0d-2 !Quad point temperature (C)

  !Ice: 
  PetscReal, parameter :: ICE_DENSITY_KG = 920.d0 !kg/m^3
  PetscReal, parameter :: ICE_DENSITY = 50.86d0 !mol/L
  PetscReal, parameter :: lambda_hyd = 0.49d0 !W/m-K

  PetscInt, parameter, public :: HYDRATE_FORMER_CH4 = ZERO_INTEGER
  PetscInt, parameter, public :: HYDRATE_FORMER_CO2 = ONE_INTEGER

  PetscInt, public :: hydrate_perm_scaling_function = 0
  PetscInt, public :: hydrate_phase_boundary = 1
  PetscInt, public :: hydrate_henrys_constant = 1
  PetscBool, public :: hydrate_perm_scaling = PETSC_FALSE
  PetscBool, public :: hydrate_eff_sat_scaling = PETSC_FALSE
  PetscBool, public :: hydrate_with_gibbs_thomson = PETSC_FALSE
  PetscBool, public :: hydrate_gt_3phase = PETSC_FALSE
  PetscBool, public :: hydrate_adjust_ghsz_solubility = PETSC_FALSE
  PetscBool, public :: hydrate_with_sedimentation = PETSC_FALSE
  PetscBool, public :: hydrate_no_pc = PETSC_FALSE
  PetscBool, public :: hydrate_with_methanogenesis = PETSC_FALSE
  PetscInt, public :: hydrate_max_states = UNINITIALIZED_INTEGER
  PetscInt, public :: hydrate_former = ZERO_INTEGER
  PetscInt, public :: hydrate_mass_balance_1 = ZERO_INTEGER
  PetscInt, public :: hydrate_mass_balance_2 = ONE_INTEGER
  PetscBool, public :: hydrate_salinity = PETSC_FALSE
  PetscBool, public :: hydrate_three_component = PETSC_FALSE 

  type, public :: hydrate_auxvar_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal, pointer :: den(:)    ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: den_kg(:) ! (iphase) kg/m^3 phase
    PetscReal :: temp
    PetscReal, pointer :: xmol(:,:) ! (icomp,iphase)
    PetscReal, pointer :: H(:) ! MJ/kmol
    PetscReal, pointer :: U(:) ! MJ/kmol
    PetscReal, pointer :: kr(:)
    PetscReal, pointer :: mobility(:) ! relative perm / kinematic viscosity
    PetscReal, pointer :: Pe(:) ! equilibrium phase boundary pressure
    PetscReal :: K_H_tilde
    PetscReal :: K_H_tilde2
    PetscReal :: K_H_tilde_hyd
    PetscReal :: Tf_ice ! freezing point of ice can be function of state vars
    PetscReal :: effective_porosity ! factors in compressibility
    PetscReal :: perm_base
    PetscReal :: v_sed
    PetscReal :: srl
    PetscReal :: srg
    PetscReal :: pert
    PetscReal :: h1_h2_ratio
    PetscBool :: istatechng
    type(hydrate_derivative_auxvar_type), pointer :: d
  end type hydrate_auxvar_type

  type, public :: hydrate_derivative_auxvar_type
    PetscReal :: pc_satg
    PetscReal :: por_p
    PetscReal :: denl_pl
    PetscReal :: denl_T
    PetscReal :: deng_pg
    PetscReal :: deng_pa
    PetscReal :: deng_T    
    PetscReal :: dengkg_pg
    PetscReal :: dengkg_T    
    PetscReal :: Ul_pl
    PetscReal :: Ul_T
    PetscReal :: Hl_pl
    PetscReal :: Hl_T
    PetscReal :: Ug_pg
    PetscReal :: Ug_pa
    PetscReal :: Ug_T
    PetscReal :: Hg_pg
    PetscReal :: Hg_pa
    PetscReal :: Hg_T 

    PetscReal :: Hv
    PetscReal :: Ha
    PetscReal :: Uv
    PetscReal :: Ua
    PetscReal :: Hv_pg
    PetscReal :: Hv_pa
    PetscReal :: Hv_T 
    PetscReal :: Ha_pg
    PetscReal :: Ha_pa
    PetscReal :: Ha_T 
    PetscReal :: Uv_pg
    PetscReal :: Uv_pa
    PetscReal :: Uv_T 
    PetscReal :: Ua_pg
    PetscReal :: Ua_pa
    PetscReal :: Ua_T 
    PetscReal :: denv
    PetscReal :: dena 
    PetscReal :: denv_T 
    PetscReal :: dena_T 
    PetscReal :: denv_pg
    PetscReal :: dena_pg
    PetscReal :: Hc 
    
    PetscReal :: psat_p
    PetscReal :: psat_T
    PetscReal :: pv_p
    PetscReal :: pv_pa
    PetscReal :: pv_T
    PetscReal :: Hc_p
    PetscReal :: Hc_T
    PetscReal :: mobilityl_pl
    PetscReal :: mobilityl_T
    PetscReal :: mobilityl_satg
    PetscReal :: mobilityg_pg
    PetscReal :: mobilityg_T
    PetscReal :: mobilityg_satg
    PetscReal :: mobilityg_pa
    PetscReal :: mug
    PetscReal :: mug_T
    PetscReal :: mug_pg
    PetscReal :: xmol_p(2,2)
    PetscReal :: xmol_T(2,2)    
  end type hydrate_derivative_auxvar_type
  
  type, public :: hydrate_parameter_type
    PetscReal, pointer :: diffusion_coefficient(:) ! (iphase)
    PetscReal :: newton_inf_scaled_res_tol
    PetscBool :: check_post_converged
    type(methanogenesis_type), pointer :: methanogenesis
  end type hydrate_parameter_type

  type, public :: methanogenesis_type
    character(len=MAXWORDLENGTH) :: source_name
    PetscReal :: alpha
    PetscReal :: k_alpha
    PetscReal :: lambda
    PetscReal :: omega
    PetscReal :: z_smt
    type(methanogenesis_type), pointer :: next
  end type methanogenesis_type
  
  type, public :: hydrate_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: inactive_cells_exist
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(hydrate_parameter_type), pointer :: hydrate_parameter
    type(hydrate_auxvar_type), pointer :: auxvars(:,:)
    type(hydrate_auxvar_type), pointer :: auxvars_bc(:)
    type(hydrate_auxvar_type), pointer :: auxvars_ss(:)
    type(matrix_zeroing_type), pointer :: matrix_zeroing
  end type hydrate_type

  interface HydrateAuxVarDestroy
    module procedure HydrateAuxVarSingleDestroy
    module procedure HydrateAuxVarArray1Destroy
    module procedure HydrateAuxVarArray2Destroy
  end interface HydrateAuxVarDestroy
  
  interface HydrateOutputAuxVars
    module procedure HydrateOutputAuxVars1
    module procedure HydrateOutputAuxVars2
  end interface HydrateOutputAuxVars
  
  public :: HydrateAuxCreate, &
            HydrateMethanogenesisCreate, &
            HydrateAuxDestroy, &
            HydrateAuxSetEnergyDOF, &
            HydrateAuxSetHydrateFormer, &
            HydrateAuxVarCompute3, &
            HydrateAuxVarCompute4, &
            HydrateComputeStateSpecificVars3, &
            HydrateComputeStateSpecificVars4, &
            HydrateAuxVarInit, &
            HydrateAuxVarCopy, &
            HydrateAuxVarDestroy, &
            HydrateAuxVarStrip, &
            HydrateAuxVarUpdateState, &
            HydrateAuxVarPerturb, &
            HydratePrintAuxVars, &
            HydrateOutputAuxVars, &
            HydrateCompositeThermalCond,&
            HydratePE, &
            HydrateMethanogenesis, &
            HenrysConstant, &
            HydrateGHSZSolubilityCorrection, &
            HydrateScalePermeability, &
            GibbsThomsonFreezing, &
            EOSIceEnergy, &
            EOSHydrateEnthalpy


contains

! ************************************************************************** !

function HydrateAuxCreate(option)
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Option_module

  implicit none

  type(option_type) :: option
    
  type(hydrate_type), pointer :: HydrateAuxCreate
  
  type(hydrate_type), pointer :: aux

  allocate(aux) 
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%inactive_cells_exist = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%matrix_zeroing)

  allocate(aux%hydrate_parameter)
  allocate(aux%hydrate_parameter%diffusion_coefficient(option%nphase))
  !geh: there is no point in setting default lquid diffusion coeffcient values 
  !     here as they will be overwritten by the fluid property defaults.
  aux%hydrate_parameter%diffusion_coefficient(LIQUID_PHASE) = &
                                                           UNINITIALIZED_DOUBLE
  aux%hydrate_parameter%diffusion_coefficient(GAS_PHASE) = 2.13d-5
  aux%hydrate_parameter%newton_inf_scaled_res_tol = 1.d-50
  aux%hydrate_parameter%check_post_converged = PETSC_FALSE

  HydrateAuxCreate => aux
  
end function HydrateAuxCreate

! ************************************************************************** !

function HydrateMethanogenesisCreate()

  ! 
  ! Allocate and initialize methanogenesis object
  ! 
  ! Author: Michael Nole
  ! Date: 11/21/19
  ! 

  type(methanogenesis_type), pointer :: HydrateMethanogenesisCreate
  type(methanogenesis_type), pointer :: methanogenesis

  allocate(methanogenesis)

  methanogenesis%source_name = ''
  methanogenesis%alpha = UNINITIALIZED_DOUBLE
  methanogenesis%k_alpha = UNINITIALIZED_DOUBLE
  methanogenesis%lambda = UNINITIALIZED_DOUBLE
  methanogenesis%omega = UNINITIALIZED_DOUBLE
  methanogenesis%z_smt = UNINITIALIZED_DOUBLE

  HydrateMethanogenesisCreate => methanogenesis

end function HydrateMethanogenesisCreate

! ************************************************************************** !
subroutine HydrateAuxVarInit(auxvar,allocate_derivative,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Option_module

  implicit none
  
  type(hydrate_auxvar_type) :: auxvar
  PetscBool :: allocate_derivative
  type(option_type) :: option

  auxvar%istate_store = NULL_STATE
  auxvar%temp = 0.d0
  auxvar%effective_porosity = 0.d0
  auxvar%perm_base = -999.9d0
  auxvar%v_sed = 0.d0
  auxvar%srl = 0.d0
  auxvar%srg = 0.d0
  auxvar%pert = 0.d0
  auxvar%h1_h2_ratio = 0.d0
  auxvar%Tf_ice = TQD
  auxvar%K_H_tilde = UNINITIALIZED_DOUBLE
  auxvar%K_H_tilde2 = UNINITIALIZED_DOUBLE
  auxvar%K_H_tilde_hyd = UNINITIALIZED_DOUBLE
  auxvar%istatechng = PETSC_FALSE
  
  allocate(auxvar%pres(option%nphase+FOUR_INTEGER))
  auxvar%pres = 0.d0
  allocate(auxvar%sat(option%nphase))
  auxvar%sat = 0.d0
  allocate(auxvar%den(option%nphase))
  auxvar%den = 0.d0
  allocate(auxvar%den_kg(option%nphase))
  auxvar%den_kg = 0.d0
  allocate(auxvar%xmol(option%nflowspec,option%nphase))
  auxvar%xmol = 0.d0
  allocate(auxvar%H(option%nphase))
  auxvar%H = 0.d0
  allocate(auxvar%U(option%nphase))
  auxvar%U = 0.d0
  allocate(auxvar%mobility(option%nphase))
  auxvar%mobility = 0.d0
  allocate(auxvar%kr(option%nphase))
  auxvar%kr = 0.d0
  allocate(auxvar%Pe(option%nflowspec-1))
  auxvar%Pe = 0.d0
  if (allocate_derivative) then
    allocate(auxvar%d)
    auxvar%d%pc_satg = 0.d0
    auxvar%d%por_p = 0.d0
    auxvar%d%denl_pl = 0.d0
    auxvar%d%denl_T = 0.d0
    auxvar%d%deng_pg = 0.d0
    auxvar%d%deng_pa = 0.d0
    auxvar%d%deng_T = 0.d0
    auxvar%d%dengkg_pg = 0.d0
    auxvar%d%dengkg_T = 0.d0
    auxvar%d%Ul_pl = 0.d0
    auxvar%d%Ul_T = 0.d0
    auxvar%d%Hl_pl = 0.d0
    auxvar%d%Hl_T = 0.d0
    auxvar%d%Ug_pg = 0.d0
    auxvar%d%Ug_pa = 0.d0
    auxvar%d%Ug_T = 0.d0
    auxvar%d%Hg_pg = 0.d0
    auxvar%d%Hg_pa = 0.d0
    auxvar%d%Hg_T = 0.d0

    auxvar%d%Hv = 0.d0
    auxvar%d%Ha = 0.d0
    auxvar%d%Uv = 0.d0
    auxvar%d%Ua = 0.d0
    auxvar%d%Hv_pg = 0.d0
    auxvar%d%Hv_pa = 0.d0
    auxvar%d%Hv_T = 0.d0
    auxvar%d%Ha_pg = 0.d0
    auxvar%d%Ha_pa = 0.d0
    auxvar%d%Ha_T = 0.d0
    auxvar%d%Uv_pg = 0.d0
    auxvar%d%Uv_pa = 0.d0
    auxvar%d%Uv_T = 0.d0
    auxvar%d%Ua_pg = 0.d0
    auxvar%d%Ua_pa = 0.d0
    auxvar%d%Ua_T = 0.d0
    auxvar%d%denv = 0.d0
    auxvar%d%dena = 0.d0
    auxvar%d%denv_T = 0.d0
    auxvar%d%dena_T = 0.d0
    auxvar%d%denv_pg = 0.d0
    auxvar%d%dena_pg = 0.d0
    auxvar%d%Hc = 0.d0
        
    auxvar%d%psat_p = 0.d0
    auxvar%d%psat_T = 0.d0
    auxvar%d%pv_p = 0.d0
    auxvar%d%pv_pa = 0.d0
    auxvar%d%pv_T = 0.d0
    auxvar%d%Hc_p = 0.d0
    auxvar%d%Hc_T = 0.d0
    auxvar%d%mug = 0.d0
    auxvar%d%mug_T = 0.d0
    auxvar%d%mug_pg = 0.d0
    auxvar%d%mobilityl_pl = 0.d0
    auxvar%d%mobilityl_T = 0.d0
    auxvar%d%mobilityl_satg = 0.d0
    auxvar%d%mobilityg_pg = 0.d0
    auxvar%d%mobilityg_T = 0.d0
    auxvar%d%mobilityg_satg = 0.d0
    auxvar%d%mobilityg_pa = 0.d0
    auxvar%d%xmol_p = 0.d0
    auxvar%d%xmol_T = 0.d0
  else
    nullify(auxvar%d)
  endif
  
end subroutine HydrateAuxVarInit

! ************************************************************************** !

subroutine HydrateAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Option_module

  implicit none
  
  type(hydrate_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%istate_store = auxvar%istate_store
  auxvar2%pres = auxvar%pres
  auxvar2%temp = auxvar%temp
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%xmol = auxvar%xmol
  auxvar2%H = auxvar%H
  auxvar2%U = auxvar%U
  auxvar2%mobility = auxvar%mobility
  auxvar2%kr = auxvar%kr
  auxvar2%Pe = auxvar%Pe
  auxvar2%perm_base = auxvar%perm_base
  auxvar2%v_sed = auxvar%v_sed
  auxvar2%srl = auxvar%srl
  auxvar2%srg = auxvar%srg
  auxvar2%h1_h2_ratio = auxvar%h1_h2_ratio
  auxvar2%Tf_ice = auxvar%Tf_ice
  auxvar2%K_H_tilde = auxvar%K_H_tilde
  auxvar2%K_H_tilde2 = auxvar%K_H_tilde2
  auxvar2%K_H_tilde_hyd = auxvar%K_H_tilde_hyd
  auxvar2%effective_porosity = auxvar%effective_porosity
  auxvar2%pert = auxvar%pert


end subroutine HydrateAuxVarCopy

! ************************************************************************** !

subroutine HydrateAuxSetEnergyDOF(energy_keyword,option)
  ! 
  ! Sets the two phase primary dependent variable for energy based on user 
  ! input.
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
  use Option_module
  use String_module

  implicit none
  
  character(len=MAXWORDLENGTH) :: energy_keyword
  type(option_type) :: option

  call StringToUpper(energy_keyword)
  select case(energy_keyword)
    case('TEMPERATURE')
      hydrate_2ph_energy_dof = HYDRATE_TEMPERATURE_INDEX
    case('AIR_PRESSURE')
      hydrate_2ph_energy_dof = HYDRATE_AIR_PRESSURE_INDEX
    case default
      option%io_buffer = 'Energy Keyword: ' // trim(energy_keyword) // &
                          ' not recognized in Hydrate Mode'    
      call PrintErrMsg(option)
  end select

end subroutine HydrateAuxSetEnergyDOF

! ************************************************************************** !

subroutine HydrateAuxSetHydrateFormer(hydrate_former_input,option)
  !
  ! Sets the molecule that will form hydrate (CO2 or CH4)
  !
  ! Author: David Fukuyama
  ! Date: 09/11/2019
  !
  use Option_module
  use String_module

  implicit none

  character(len=MAXWORDLENGTH) :: hydrate_former_input
  type(option_type) :: option

  call StringToUpper(hydrate_former_input)
  select case(hydrate_former_input)
    case('CO2')
      HYDRATE_DENSITY_KG = 1100.d0 !kg/m^3
      HYDRATE_DENSITY = 50.31596d0 !mol/L

      MOL_RATIO_FORMER = 0.148148148d0
      MOL_RATIO_H2O = 1.d0 - MOL_RATIO_FORMER

      TQD = 271.80 - 273.15  !Quad point temp, [C]
      !TQD2 = 283.10 - 273.15

      hydrate_former = HYDRATE_FORMER_CO2

    case('CH4')
      HYDRATE_DENSITY_KG = 920.d0 !kg/m^3
      HYDRATE_DENSITY = 52.15551276d0 !mol/L

      MOL_RATIO_FORMER = 0.14285714285d0
      MOL_RATIO_H2O = 1.d0 - MOL_RATIO_FORMER

      !TQD = 1.d-2 !Quad point temp, [C]

      hydrate_former = HYDRATE_FORMER_CH4

    case default
      HYDRATE_DENSITY_KG = 920.d0 !kg/m^3
      HYDRATE_DENSITY = 52.15551276d0 !mol/L

      MOL_RATIO_FORMER = 0.14285714285d0
      MOL_RATIO_H2O = 1.d0 - MOL_RATIO_FORMER

      TQD = 1.d-2 !Quad point temp, [C]

      hydrate_former = HYDRATE_FORMER_CH4
  end select

end subroutine HydrateAuxSetHydrateFormer

! ************************************************************************** !

subroutine HydrateAuxVarCompute3(x,hyd_auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell, with gas hydrate physics
  ! Author: Michael Nole
  ! Date: 04/04/19
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(hydrate_auxvar_type) :: hyd_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id

  PetscInt :: gid, lid, acid, wid, eid, hid, iid
  PetscReal :: cell_pressure, water_vapor_pressure
  PetscReal :: den_water_vapor, den_kg_water_vapor
  PetscReal :: u_water_vapor, h_water_vapor
  PetscReal :: den_air, h_air, u_air
  PetscReal :: xmol_air_in_gas, xmol_water_in_gas
  PetscReal :: krl, visl
  PetscReal :: dkrl_dsatl, dkrl_dsatg
  PetscReal :: dkrg_dsatl, dkrg_dsatg
  PetscReal :: krg, visg
  PetscReal :: K_H_tilde, K_H_tilde_hyd
  PetscInt :: apid, cpid, vpid, spid
  PetscReal :: xmass_air_in_gas
  PetscReal :: Ugas_J_kg, Hgas_J_kg
  PetscReal :: Uair_J_kg, Hair_J_kg
  PetscReal :: Uvapor_J_kg, Hvapor_J_kg
  PetscReal :: Hg_mixture_fractioned
  PetscReal :: H_hyd, U_ice, PE_hyd
  PetscReal :: aux(1)
  PetscReal :: hw, hw_dp, hw_dT
  PetscReal :: dpor_dp
  PetscReal :: dpc_dsatl
  PetscReal :: dden_ice_dT, dden_ice_dP
  character(len=8) :: state_char
  PetscErrorCode :: ierr, eos_henry_ierr
  PetscReal :: dTf, h_sat_eff, i_sat_eff, liq_sat_eff, gas_sat_eff
  PetscReal :: solid_sat_eff
  PetscReal :: sigma, dP

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id
  acid = option%air_id
  wid = option%water_id
  eid = option%energy_id


  !geh: do not initialize hyd_auxvar%temp as the previous value is used as the
  !     initial guess for two phase.
  hyd_auxvar%H = 0.d0
  hyd_auxvar%U = 0.d0
  hyd_auxvar%pres = 0.d0
  hyd_auxvar%sat = 0.d0
  hyd_auxvar%den = 0.d0
  hyd_auxvar%den_kg = 0.d0
  hyd_auxvar%xmol = 0.d0
  hyd_auxvar%effective_porosity = 0.d0
  hyd_auxvar%mobility = 0.d0
  hyd_auxvar%kr = 0.d0

#if 0
  if (option%iflag >= HYDRATE_UPDATE_FOR_ACCUM) then
    if (option%iflag == HYDRATE_UPDATE_FOR_ACCUM) then
      write(*,'(a,i3,3es17.8,a3)') 'before: ', &
        natural_id, x(1:3), trim(state_char)
    else
    endif
  endif
#endif

  hyd_auxvar%xmol(wid,hid) = MOL_RATIO_H2O
  hyd_auxvar%xmol(acid,hid) = MOL_RATIO_METH
  hyd_auxvar%den(hid) = HYDRATE_DENSITY
  hyd_auxvar%den_kg(hid) = HYDRATE_DENSITY_KG

  call HydrateComputeStateSpecificVars3(x,global_auxvar%istate, &
                                        global_auxvar%m_nacl(1), hyd_auxvar, &
                                        characteristic_curves, &
                                        material_auxvar,natural_id,option)
  !MAN Check EOS temperature
  !eos_henry_ierr = 0
  !call EOSGasHenry(hyd_auxvar%temp,hyd_auxvar%pres(spid),K_H_tilde, &
  !                       eos_henry_ierr)
  !if (eos_henry_ierr /= 0) then
  !   call HydrateEOSGasError(natural_id,eos_henry_ierr,hyd_auxvar,option)
  !endif

  cell_pressure = max(hyd_auxvar%pres(lid),hyd_auxvar%pres(gid), &
                      hyd_auxvar%pres(spid))

  ! calculate effective porosity as a function of pressure
  if (option%iflag /= HYDRATE_UPDATE_FOR_BOUNDARY) then
    dpor_dp = 0.d0
    hyd_auxvar%effective_porosity = material_auxvar%porosity_base
    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                hyd_auxvar%effective_porosity,dpor_dp)
    endif
    if (option%iflag /= HYDRATE_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = hyd_auxvar%effective_porosity
    endif
    !MAN: does this need to be in this conditional?
    solid_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid)

    if (hydrate_perm_scaling) then
      call HydrateScalePermeability(hyd_auxvar%perm_base, &
                                    material_auxvar%permeability, &
                                    solid_sat_eff,option)
    endif

  endif
  if (associated(hyd_auxvar%d)) then
    hyd_auxvar%d%por_p = dpor_dp
  endif

  ! ALWAYS UPDATE THERMODYNAMIC PROPERTIES

  ! Liquid phase thermodynamic properties
  ! must use cell_pressure as the pressure, not %pres(lid)
  if (.not.option%flow%density_depends_on_salinity) then
    if (associated(hyd_auxvar%d)) then
      call EOSWaterDensity(hyd_auxvar%temp,cell_pressure, &
                           hyd_auxvar%den_kg(lid),hyd_auxvar%den(lid), &
                           hyd_auxvar%d%denl_pl,hyd_auxvar%d%denl_T,ierr)
    else
      call EOSWaterDensity(hyd_auxvar%temp,cell_pressure, &
                           hyd_auxvar%den_kg(lid),hyd_auxvar%den(lid),ierr)
    endif
  else
    aux(1) = global_auxvar%m_nacl(1)
    call EOSWaterDensityExt(hyd_auxvar%temp,cell_pressure,aux, &
                              hyd_auxvar%den_kg(lid),hyd_auxvar%den(lid),ierr)
  endif

  call EOSWaterEnthalpy(hyd_auxvar%temp,cell_pressure,hw,ierr)

  hyd_auxvar%H(lid) = hw * 1.d-6 ! J/kmol -> MJ/kmol
  ! MJ/kmol comp
  hyd_auxvar%U(lid) = hyd_auxvar%H(lid) - &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        (cell_pressure / hyd_auxvar%den(lid) * &
                        1.d-6)
  if (global_auxvar%istate /= L_STATE) then
    if (global_auxvar%istate == G_STATE) then
      water_vapor_pressure = hyd_auxvar%pres(vpid)
    else
      water_vapor_pressure = hyd_auxvar%pres(spid)
    endif
    call EOSGasDensityEnergy(hyd_auxvar%temp,hyd_auxvar%pres(apid),den_air, &
                               h_air,u_air,ierr)
    h_air = h_air * 1.d-6 ! J/kmol -> MJ/kmol
    u_air = u_air * 1.d-6 ! J/kmol -> MJ/kmol

    call EOSWaterSteamDensityEnthalpy(hyd_auxvar%temp,water_vapor_pressure, &
                                        den_kg_water_vapor,den_water_vapor, &
                                        h_water_vapor,ierr)
    u_water_vapor = h_water_vapor - &
                    ! Pa / kmol/m^3 = J/kmol
                    water_vapor_pressure / den_water_vapor
    h_water_vapor = h_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    u_water_vapor = u_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    hyd_auxvar%den(gid) = den_water_vapor + den_air
    hyd_auxvar%den_kg(gid) = den_kg_water_vapor + den_air*fmw_comp(gid)
    xmol_air_in_gas = hyd_auxvar%xmol(acid,gid)
    xmol_water_in_gas = hyd_auxvar%xmol(wid,gid)

    ! MJ/kmol
    hyd_auxvar%U(gid) = xmol_water_in_gas * u_water_vapor + &
                        xmol_air_in_gas * u_air
    Hg_mixture_fractioned = xmol_water_in_gas*h_water_vapor + &
                            xmol_air_in_gas*h_air
    hyd_auxvar%H(gid) = hyd_auxvar%U(gid) + &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        hyd_auxvar%pres(gid)/hyd_auxvar%den(gid) * 1.d-6

  endif

  if (hydrate_eff_sat_scaling) then
    if (hyd_auxvar%sat(lid) > 0.d0 .or. hyd_auxvar%sat(gid) > 0.d0) then
      liq_sat_eff = hyd_auxvar%sat(lid) / (hyd_auxvar%sat(lid)+hyd_auxvar%sat(gid))
      gas_sat_eff = 1.d0 - liq_sat_eff
    else
      liq_sat_eff = 0.d0
      gas_sat_eff = 0.d0
    endif
  else
    liq_sat_eff = hyd_auxvar%sat(lid)
    gas_sat_eff = hyd_auxvar%sat(gid)
  endif

  if (liq_sat_eff > 0.d0) then
    if (liq_sat_eff >= 1.d0) then
      krl = 1.d0
    else
      call characteristic_curves%liq_rel_perm_function% &
           RelativePermeability(liq_sat_eff,krl,dkrl_dsatl,option)
      krl = max(0.d0,krl)
    endif
    call EOSWaterViscosity(hyd_auxvar%temp,cell_pressure, &
                               hyd_auxvar%pres(spid),visl,ierr)
    hyd_auxvar%mobility(lid) = krl/visl
    hyd_auxvar%kr(lid) = krl
  endif

  if (gas_sat_eff > 0.d0) then
    if (gas_sat_eff >=1.d0) then
      krg = 1.d0
    else
      call characteristic_curves%gas_rel_perm_function% &
           RelativePermeability(liq_sat_eff,krg,dkrg_dsatl,option)
      krg = max(0.d0,krg)
    endif
    call EOSGasViscosity(hyd_auxvar%temp,hyd_auxvar%pres(apid), &
                           hyd_auxvar%pres(gid),den_air,visg,ierr)
    hyd_auxvar%mobility(gid) = krg/visg
    hyd_auxvar%kr(gid) = krg
  endif

  call EOSHydrateEnthalpy(hyd_auxvar%temp, H_hyd, HYDRATE_FORMER_CH4)
  hyd_auxvar%U(hid) = H_hyd !- cell_pressure/hyd_auxvar%den(hid)*1.d-6
  hyd_auxvar%H(hid) = H_hyd
  hyd_auxvar%mobility(hid) = 0.d0

  call EOSIceEnergy(hyd_auxvar%temp, U_ice)
  hyd_auxvar%xmol(wid,iid) = 1.d0
  hyd_auxvar%xmol(gid,iid) = 0.d0
  !call EOSWaterDensityIcePainter(hyd_auxvar%temp,hyd_auxvar%pres(lid), &
  !                  PETSC_FALSE, hyd_auxvar%den(iid), &
  !                  dden_ice_dT, dden_ice_dP, ierr)
  hyd_auxvar%den(iid) = ICE_DENSITY
  hyd_auxvar%den_kg(iid) = ICE_DENSITY_KG
  hyd_auxvar%U(iid) = U_ice
  hyd_auxvar%H(iid) = U_ice
  hyd_auxvar%mobility(iid) = 0.d0

  hyd_auxvar%srl = characteristic_curves%gas_rel_perm_function%sr
  hyd_auxvar%srg = characteristic_curves%gas_rel_perm_function%srg

end subroutine HydrateAuxVarCompute3

! ************************************************************************** !

subroutine HydrateComputeStateSpecificVars3(x,state,m_nacl,hyd_auxvar, &
                  characteristic_curves,material_auxvar,natural_id,option)
  !
  ! Computes auxvars depending on phase state for 3 equations.
  !
  ! Author: Michael Nole
  ! Date: 08/26/20
  !

  use Option_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use co2eos_module
  use co2_span_wagner_module
  use co2_span_wagner_spline_module

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  PetscInt :: state
  PetscReal :: m_nacl
  type(hydrate_auxvar_type) :: hyd_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id

  PetscReal :: PE_hyd, dP, dP2, K_H_tilde, K_H_tilde2, K_H_tilde_hyd, &
               K_H_tilde_hyd2,sigma
  PetscInt :: gid, lid, acid, wid, eid, hid, iid, h2id, cid, &
              apid, cpid, vpid, spid
  PetscReal :: dTf, h_sat_eff, i_sat_eff, liq_sat_eff, gas_sat_eff
  PetscReal :: dpc_dsatl
  PetscInt :: ierr

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  h2id = option%hydrate2_phase ! co2 hydrate

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id
  wid = option%water_id
  cid = option%co2_id
  eid = option%energy_id

  PE_hyd = 0.d0

  select case(state)
    case(L_STATE)
!     ********* Aqueous State (A) ********************************
!     Primary variables: Pa, Xma, T
!
      hyd_auxvar%pres(lid) = x(HYDRATE_LIQUID_PRESSURE_DOF)
      hyd_auxvar%xmol(acid,lid) = x(HYDRATE_L_STATE_X_MOLE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%xmol(acid,lid) = max(0.d0,hyd_auxvar%xmol(acid,lid))

      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = 0.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0
      hyd_auxvar%sat(lid) = 1.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                        hyd_auxvar%pres(spid),ierr)
      call HydratePE(hyd_auxvar%temp, 0.d0, m_nacl, PE_hyd, &
                     dP, characteristic_curves, material_auxvar,option, &
                     hydrate_former)
      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      K_H_tilde_hyd = K_H_tilde

      if (hydrate_adjust_ghsz_solubility) then
        call HydrateGHSZSolubilityCorrection(hyd_auxvar%temp,hyd_auxvar% &
                                           pres(lid),dP,K_H_tilde_hyd)
      endif

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(gid) = max(hyd_auxvar%pres(lid),hyd_auxvar%pres(spid))
      hyd_auxvar%pres(apid) = K_H_tilde_hyd*hyd_auxvar%xmol(acid,lid)

      if (hyd_auxvar%pres(gid) <= 0.d0) then
        write(option%io_buffer,'(''Negative gas pressure at cell '', &
          & i8,'' in HydrateAuxVarCompute(L_STATE).  Attempting bailout.'')') &
          natural_id
        call PrintMsgByRank(option)
        hyd_auxvar%pres(vpid) = 0.5d0*hyd_auxvar%pres(spid)
        hyd_auxvar%pres(gid) = hyd_auxvar%pres(vpid) + hyd_auxvar%pres(apid)
      else
        hyd_auxvar%pres(vpid) = hyd_auxvar%pres(lid) - hyd_auxvar%pres(apid)
      endif
      hyd_auxvar%pres(cpid) = 0.d0
    case (G_STATE)
!     ********* Gas State (G) ********************************
!     Primary variables: Pg, Pa, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)

      hyd_auxvar%pres(apid) = x(HYDRATE_G_STATE_AIR_PRESSURE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 1.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0

      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / &
                                   hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                    hyd_auxvar%pres(spid),ierr)
      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 0.d0

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call characteristic_curves%saturation_function% &
             CapillaryPressure(hyd_auxvar%sat(lid), &
                               hyd_auxvar%pres(cpid),dpc_dsatl,option)
      endif

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - &
                             hyd_auxvar%pres(cpid)

    case (H_STATE)
!     ********* Hydrate State (H) ********************************
!     Primary variables: Pg, Xmh, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      x(HYDRATE_GAS_SATURATION_DOF) = MOL_RATIO_METH
      hyd_auxvar%xmol(acid,hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 1.d0
      hyd_auxvar%sat(iid) = 0.d0

      call HydratePE(hyd_auxvar%temp,hyd_auxvar%sat(hid), &
                     m_nacl, PE_hyd,dP,characteristic_curves, &
                     material_auxvar, option, hydrate_former)
      hyd_auxvar%pres(apid) = PE_hyd
      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

    case(I_STATE)
!     ********* Ice State (I) ********************************
!     Primary variables: Pg, Xmi, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      x(HYDRATE_GAS_SATURATION_DOF) = 0.d0
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 1.d0

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(apid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%xmol(acid,lid) = 0.d0
      hyd_auxvar%xmol(wid,lid) = 1.d0
      hyd_auxvar%xmol(wid,gid) = 1.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0

    case(GA_STATE)
!     ********* Gas & Aqueous State (GA) ********************************
!     Primary variables: Pg, Sg, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(gid) = max(0.d0,min(1.d0,hyd_auxvar%sat(gid)))
      hyd_auxvar%sat(lid) = 1.d0 - hyd_auxvar%sat(gid)
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)
      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(vpid)

      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call characteristic_curves%saturation_function% &
             CapillaryPressure(hyd_auxvar%sat(lid), hyd_auxvar%pres(cpid), &
                               dpc_dsatl,option)
      endif

      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(hyd_auxvar%temp,sigma)
      endif
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid)*sigma

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

    case(HG_STATE)
!     ********* Hydrate & Gas State (HG) ********************************
!     Primary variables: Pg, Sg, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(hid) = 1.d0 - hyd_auxvar%sat(gid)
      hyd_auxvar%sat(iid) = 0.d0

      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
        call HydratePE(hyd_auxvar%temp, hyd_auxvar%sat(hid)+ &
                hyd_auxvar%sat(gid), m_nacl, &
                PE_hyd, dP, characteristic_curves, &
                material_auxvar, option, hydrate_former)
      else
        call HydratePE(hyd_auxvar%temp, 2.d0 * hyd_auxvar%sat(hid), &
                m_nacl, PE_hyd, dP, characteristic_curves, &
                material_auxvar, option, hydrate_former)
      endif

      hyd_auxvar%pres(apid) = PE_hyd
      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

    case(HA_STATE)
!     ********* Hydrate & Aqueous State (HA) ********************************
!     Primary variables: Pg, Sh, T
!

      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(hid) = max(0.d0,min(1.d0,hyd_auxvar%sat(hid)))

      hyd_auxvar%sat(lid) = 1.d0 - hyd_auxvar%sat(hid)
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0

      call HydratePE(hyd_auxvar%temp, 0.d0, m_nacl, PE_hyd, &
                     dP, characteristic_curves, material_auxvar,option, &
                     hydrate_former)
      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      if (hydrate_adjust_ghsz_solubility) then
        call HydrateGHSZSolubilityCorrection(hyd_auxvar%temp,hyd_auxvar% &
                                           pres(gid),dP,K_H_tilde)
      endif

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)

      !hyd_auxvar%pres(spid) = 1.d-6
      hyd_auxvar%pres(cpid) = 0.d0
      ! Setting air pressure equal to gas pressure makes forming hydrate
      ! easier
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(spid)
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = 0.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0

    case(HI_STATE)
!     ********* Hydrate & Ice State (HI) ********************************
!     Primary variables: Pg, Sh, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(iid) = 1.d0 - hyd_auxvar%sat(hid)

      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
        call HydratePE(hyd_auxvar%temp, hyd_auxvar%sat(hid)+ &
                hyd_auxvar%sat(iid), m_nacl, PE_hyd, dP, &
                characteristic_curves, material_auxvar,option, hydrate_former)
      else
        call HydratePE(hyd_auxvar%temp, 2.d0 * hyd_auxvar%sat(hid), &
                       m_nacl, PE_hyd, dP, &
                       characteristic_curves, material_auxvar,option, &
                       hydrate_former)
      endif

      hyd_auxvar%pres(apid) = PE_hyd
      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)
    case(GI_STATE)
!     ********* Gas & Ice State (GI) ********************************
!     Primary variables: Pg, Si, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 1.d0 - hyd_auxvar%sat(iid)
      hyd_auxvar%sat(hid) = 0.d0

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(vpid)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = 0.d0
      hyd_auxvar%xmol(wid,lid) = 1.d0
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

    case(AI_STATE)
!     ********* Aqueous & Ice State (AI) ********************************
!     Primary variables: Pl, Xma, Sl
!
      hyd_auxvar%pres(lid) = x(HYDRATE_LIQUID_PRESSURE_DOF)
      hyd_auxvar%xmol(acid,lid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%sat(lid) = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%xmol(acid,lid) = max(0.d0,hyd_auxvar%xmol(acid,lid))

      if (hyd_auxvar%istatechng) then
        hyd_auxvar%sat(lid) = max(0.d0,min(1.d0,hyd_auxvar%sat(lid)))
      endif

      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 1.d0 - hyd_auxvar%sat(lid)

      if (hydrate_with_gibbs_thomson) then
        call GibbsThomsonFreezing(hyd_auxvar%sat(lid),6017.1d0,ICE_DENSITY,&
                                TQD, dTf,characteristic_curves, &
                                material_auxvar,option)
      else
        dTf = 0.d0
      endif

      hyd_auxvar%temp = TQD+dTf

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)
      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      !hyd_auxvar%pres(spid) = 1.d-6
      hyd_auxvar%pres(gid) = max(hyd_auxvar%pres(lid),hyd_auxvar%pres(spid))
      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(apid) = K_H_tilde*hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%pres(vpid) = 0.d0 !hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(acid,gid) = 0.d0
      hyd_auxvar%xmol(wid,gid) = 0.d0

    case(HGA_STATE)
!     ********* Hydrate, Gas, & Aqueous State (HGA) **************************
!     Primary variables: Sl, Sh, T
!
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      !hyd_auxvar%sat(lid) = max(0.d0,min(1.d0,hyd_auxvar%sat(lid)))
      !hyd_auxvar%sat(hid) = max(0.d0,min(1.d0,hyd_auxvar%sat(hid)))

      !if (hyd_auxvar%sat(lid) + hyd_auxvar%sat(hid) > 1.d0) then
      ! hyd_auxvar%sat(hid) = 1.d0 - hyd_auxvar%sat(lid)
      ! x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
      !endif

      hyd_auxvar%sat(gid) = max(0.d0,1.d0 - hyd_auxvar%sat(lid) - &
                                hyd_auxvar%sat(hid))
      hyd_auxvar%sat(iid) = 0.d0

      if (hydrate_gt_3phase) then
        if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
          h_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(gid)
          gas_sat_eff = 2.d0 * hyd_auxvar%sat(gid)
        else
          gas_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(gid)
          h_sat_eff = 2.d0 * hyd_auxvar%sat(hid)
        endif
      endif

      h_sat_eff = hyd_auxvar%sat(hid)

      if (hydrate_eff_sat_scaling) then
        liq_sat_eff = hyd_auxvar%sat(lid)/(hyd_auxvar%sat(lid)+ &
                    hyd_auxvar%sat(gid))
      else
        liq_sat_eff = hyd_auxvar%sat(lid)
      endif

      call HydratePE(hyd_auxvar%temp, h_sat_eff, m_nacl, &
                     PE_hyd, dP, characteristic_curves, material_auxvar, &
                     option, hydrate_former)
      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call characteristic_curves%saturation_function%CapillaryPressure( &
                liq_sat_eff, hyd_auxvar%pres(cpid), &
                dpc_dsatl,option)
      endif

      hyd_auxvar%pres(apid) = PE_hyd

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                    hyd_auxvar%pres(spid),ierr)
      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(gid) = hyd_auxvar%pres(apid) + hyd_auxvar%pres(vpid)

      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(hyd_auxvar%temp,sigma)
      endif
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid)*sigma

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

    case(HAI_STATE)
!     ********* Hydrate, Aqueous, & Ice State (HAI) **************************
!     Primary variables: Pg, Sl, Si
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 1.d0 - hyd_auxvar%sat(lid) - hyd_auxvar%sat(iid)

      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
        h_sat_eff = hyd_auxvar%sat(hid)+hyd_auxvar%sat(iid)
        i_sat_eff = 2.d0 * hyd_auxvar%sat(iid)
      else
        h_sat_eff = 2.d0 * hyd_auxvar%sat(hid)
        i_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid)
      endif

      if (hydrate_with_gibbs_thomson) then
        call GibbsThomsonFreezing(1.d0-i_sat_eff,6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,material_auxvar,option)
      else
        dTf = 0.d0
      endif

      hyd_auxvar%temp = TQD+dTf
      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)
      call HydratePE(hyd_auxvar%temp,h_sat_eff, m_nacl, &
                     PE_hyd, dP, characteristic_curves, material_auxvar, &
                     option, hydrate_former)
      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)

      hyd_auxvar%pres(cpid) = 0.d0

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(apid) = PE_hyd !hyd_auxvar%pres(gid)
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = 0.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0
    case(HGI_STATE)
!     ********* Hydrate, Gas, & Ice State (HGI) ******************************
!     Primary variables: Sh, Si, T
!
      hyd_auxvar%sat(iid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 1.d0 - hyd_auxvar%sat(hid) - hyd_auxvar%sat(iid)

      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
        if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
          call HydratePE(hyd_auxvar%temp, 1.d0, m_nacl, &
                     PE_hyd, dP, characteristic_curves, material_auxvar,  &
                     option, hydrate_former)
        else
          call HydratePE(hyd_auxvar%temp, 3.d0 * hyd_auxvar%sat(iid) + &
                  2.d0 * (hyd_auxvar%sat(hid)-hyd_auxvar%sat(iid)), &
                  m_nacl, PE_hyd, dP, characteristic_curves, &
                  material_auxvar, option, hydrate_former)
        endif
      elseif (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
        call HydratePE(hyd_auxvar%temp, 3.d0 * hyd_auxvar%sat(gid) + &
          2.d0 * (hyd_auxvar%sat(hid) - hyd_auxvar%sat(gid)), &
          m_nacl, PE_hyd, dP, characteristic_curves, &
          material_auxvar, option, hydrate_former)
      else
        call HydratePE(hyd_auxvar%temp, 3.d0 * hyd_auxvar%sat(hid), &
                       m_nacl, PE_hyd, dP, &
                       characteristic_curves,material_auxvar, option, &
                       hydrate_former)
      endif

      hyd_auxvar%pres(apid) = PE_hyd

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                    hyd_auxvar%pres(spid),ierr)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(gid) = hyd_auxvar%pres(apid) + hyd_auxvar%pres(vpid)

      hyd_auxvar%xmol(acid,lid) = 0.d0
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

    case(GAI_STATE)
!     ********* Gas, Aqueous, & Ice State (GAI) ******************************
!     Primary variables: Pg, Sl, Si
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(gid) = 1.d0 - hyd_auxvar%sat(lid) - hyd_auxvar%sat(iid)
      hyd_auxvar%sat(hid) = 0.d0

      call GibbsThomsonFreezing(1.d0-hyd_auxvar%sat(iid),6017.1d0, &
              ICE_DENSITY,TQD,dTf,characteristic_curves, material_auxvar,option)

      hyd_auxvar%temp = TQD+dTf

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                    hyd_auxvar%pres(spid),ierr)
      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(vpid)

      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call characteristic_curves%saturation_function% &
             CapillaryPressure(hyd_auxvar%sat(lid), &
                             hyd_auxvar%pres(cpid),dpc_dsatl,option)
      endif

      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(hyd_auxvar%temp,sigma)
      endif
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid)*sigma

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

    case(QUAD_STATE)
!     ********* Quadruple Point (HGAI) ********************************
!     Primary variables: Sg, Sl, Si
!
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_ENERGY_DOF)

      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
        h_sat_eff = hyd_auxvar%sat(hid)+hyd_auxvar%sat(iid)
        i_sat_eff = 2.d0 * hyd_auxvar%sat(iid)
      else
        h_sat_eff = 2.d0 * hyd_auxvar%sat(hid)
        i_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid)
      endif

      if (hydrate_with_gibbs_thomson) then
        call GibbsThomsonFreezing(1.d0-i_sat_eff,6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,material_auxvar,option)
      else
        dTf = 0.d0
      endif

      !if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
      !  if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
      !    call HydratePE(hyd_auxvar%temp, 1.d0 - hyd_auxvar%sat(lid), &
      !            PE_hyd, characteristic_curves, option)
      !  else
      !    call HydratePE(hyd_auxvar%temp, 3.d0 * hyd_auxvar%sat(iid) + &
      !            2.d0 * (hyd_auxvar%sat(hid) - hyd_auxvar%sat(iid)), PE_hyd, &
      !            characteristic_curves, option)
      !    call GibbsThomsonFreezing(3.d0 * hyd_auxvar%sat(iid), 6017.1d0, &
      !            ICE_DENSITY, TQD, dTf, characteristic_curves, option)
      !  endif
      !elseif (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
      !  call HydratePE(hyd_auxvar%temp, 3.d0 * hyd_auxvar%sat(gid) + &
      !          2.d0 * (hyd_auxvar%sat(hid) - hyd_auxvar%sat(gid)), &
      !          PE_hyd, characteristic_curves, option)
      !  call GibbsThomsonFreezing(1.d0-hyd_auxvar%sat(lid), 6017.1d0, &
      !          ICE_DENSITY, TQD, dTf, characteristic_curves, option)
      !else
      !  call HydratePE(hyd_auxvar%temp, 3.d0 *hyd_auxvar%sat(hid), &
      !          PE_hyd, characteristic_curves, option)
      !  if (hyd_auxvar%sat(iid) < hyd_auxvar%sat(gid)) then
      !    call GibbsThomsonFreezing(3.d0 * hyd_auxvar%sat(hid) + &
      !            2.d0 * (hyd_auxvar%sat(iid) - hyd_auxvar%sat(hid)), &
      !            6017.1d0, ICE_DENSITY, TQD, dTf, characteristic_curves, &
      !            option)
      !  endif
      !endif

      hyd_auxvar%temp = TQD + dTf
      call HydratePE(hyd_auxvar%temp,h_sat_eff, m_nacl, &
                     PE_hyd, dP, characteristic_curves, material_auxvar, &
                     option, hydrate_former)
      hyd_auxvar%pres(apid) = PE_hyd

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                    hyd_auxvar%pres(spid),ierr)

      call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(gid) = hyd_auxvar%pres(apid) + hyd_auxvar%pres(vpid)

      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call characteristic_curves%saturation_function% &
             CapillaryPressure(hyd_auxvar%sat(lid), &
                               hyd_auxvar%pres(cpid),dpc_dsatl,option)
      endif
      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(hyd_auxvar%temp,sigma)
      endif
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid)*sigma

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)


    case default
      write(option%io_buffer,*) state
      option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
        ') not recognized in HydrateAuxVarCompute.'
      call PrintErrMsgByRank(option)

  end select

  !Compute freezing point of water
  if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
    h_sat_eff = hyd_auxvar%sat(hid)+hyd_auxvar%sat(iid)
    i_sat_eff = 2.d0 * hyd_auxvar%sat(iid)
  else
    h_sat_eff = 2.d0 * hyd_auxvar%sat(hid)
    i_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid)
  endif

  if (hydrate_with_gibbs_thomson) then
    call GibbsThomsonFreezing(1.d0-i_sat_eff,6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,material_auxvar,option)
  else
    dTf = 0.d0
  endif

  hyd_auxvar%Tf_ice = TQD + dTf
  hyd_auxvar%Pe(ONE_INTEGER) = PE_hyd
  hyd_auxvar%K_H_tilde = K_H_tilde
  hyd_auxvar%K_H_tilde_hyd = K_H_tilde_hyd

end subroutine HydrateComputeStateSpecificVars3

! ************************************************************************** !


subroutine HydrateAuxVarCompute4(x,hyd_auxvar,global_auxvar,material_auxvar, &
                                characteristic_curves,natural_id,option)
  !
  ! Computes auxiliary variables for each grid cell, with gas hydrate physics
  ! Author: Michael Nole
  ! Date: 04/04/19
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use co2eos_module
  use co2_span_wagner_module
  use co2_span_wagner_spline_module

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(hydrate_auxvar_type) :: hyd_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id

  PetscInt :: gid, lid, acid, wid, eid, hid, iid, h2id, cid
  PetscReal :: cell_pressure, water_vapor_pressure
  PetscReal :: den_water_vapor, den_kg_water_vapor
  PetscReal :: u_water_vapor, h_water_vapor
  PetscReal :: den_air, h_air, u_air, fc, fug_air, xmole, x1m, poyn
  PetscReal :: den_co2, h_co2, u_co2
  PetscReal :: xmol_air_in_gas, xmol_water_in_gas, xmol_co2_in_gas
  PetscReal :: krl, visl
  PetscReal :: dkrl_dsatl, dkrl_dsatg
  PetscReal :: dkrg_dsatl, dkrg_dsatg
  PetscReal :: krg, visg
  PetscReal :: K_H_tilde, K_H_tilde2, K_H_tilde_hyd
  PetscInt :: apid, cpid, vpid, spid
  PetscReal :: xmass_air_in_gas
  PetscReal :: Ugas_J_kg, Hgas_J_kg
  PetscReal :: Uair_J_kg, Hair_J_kg
  PetscReal :: Uvapor_J_kg, Hvapor_J_kg
  PetscReal :: Hg_mixture_fractioned
  PetscReal :: H_hyd, U_ice, PE_hyd, PE_hyd2, H2_hyd
  PetscReal :: aux(1)
  PetscReal :: hw, hw_dp, hw_dT
  PetscReal :: dpor_dp
  PetscReal :: dpc_dsatl
  PetscReal :: dden_ice_dT, dden_ice_dP
  character(len=8) :: state_char
  PetscErrorCode :: ierr, ierr2, eos_henry_ierr, eos_henry_ierr2
  PetscReal :: dTf, h_sat_eff, i_sat_eff, liq_sat_eff, gas_sat_eff
  PetscReal :: solid_sat_eff
  PetscReal :: sigma, dP, dP2

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  h2id = option%hydrate2_phase ! co2 hydrate

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id
  wid = option%water_id
  cid = option%co2_id
  eid = option%energy_id

  !geh: do not initialize hyd_auxvar%temp as the previous value is used as the
  !     initial guess for two phase.

  hyd_auxvar%H = 0.d0
  hyd_auxvar%U = 0.d0
  hyd_auxvar%pres = 0.d0
  hyd_auxvar%sat = 0.d0
  hyd_auxvar%den = 0.d0
  hyd_auxvar%den_kg = 0.d0
  hyd_auxvar%xmol = 0.d0
  hyd_auxvar%effective_porosity = 0.d0
  hyd_auxvar%mobility = 0.d0
  hyd_auxvar%kr = 0.d0
  hyd_auxvar%Pe = 0.d0

  hyd_auxvar%xmol(wid,hid) = MOL_RATIO_H2O
  hyd_auxvar%xmol(acid,hid) = MOL_RATIO_FORMER
  if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
    hyd_auxvar%xmol(cid,h2id) = MOL_RATIO_FORMER
    hyd_auxvar%xmol(wid,h2id) = MOL_RATIO_H2O
    hyd_auxvar%den(h2id) = HYDRATE_DENSITY
    hyd_auxvar%den_kg(h2id) = HYDRATE_DENSITY_KG
    fmw_comp = [FMWH2O,FMWAIR,FMWAIR]
  elseif (hydrate_salinity) then
    fmw_comp = [FMWH2O,FMWAIR,FMWNACL]
    hyd_auxvar%xmol(cid,h2id) = MOL_RATIO_FORMER
    hyd_auxvar%xmol(wid,h2id) = MOL_RATIO_H2O
    hyd_auxvar%den(h2id) = HYDRATE_DENSITY
    hyd_auxvar%den_kg(h2id) = HYDRATE_DENSITY_KG
  else
    hyd_auxvar%xmol(cid,h2id) =MOL_RATIO_FORMER
    hyd_auxvar%xmol(wid,h2id) = MOL_RATIO_H2O
    hyd_auxvar%den(h2id) = HYDRATE_DENSITY
    hyd_auxvar%den_kg(h2id) = HYDRATE_DENSITY_KG
    fmw_comp = [FMWH2O,FMWAIR,FMWCO2]
  endif
  hyd_auxvar%den(hid) = HYDRATE_DENSITY
  hyd_auxvar%den_kg(hid) = HYDRATE_DENSITY_KG

  call HydrateComputeStateSpecificVars4(x,global_auxvar%istate,hyd_auxvar, &
                                        characteristic_curves, &
                                        material_auxvar,natural_id,option)
  !MAN Check EOS temperature
  eos_henry_ierr = 0
  eos_henry_ierr2 = 0
  call EOSGasHenry(hyd_auxvar%temp,hyd_auxvar%pres(spid),K_H_tilde, &
       eos_henry_ierr)
  call EOSGasHenry(hyd_auxvar%temp,hyd_auxvar%pres(spid),K_H_tilde2, &
       eos_henry_ierr2)
  if (eos_henry_ierr /= 0) then
     call HydrateEOSGasError(natural_id,eos_henry_ierr,hyd_auxvar,option)
  endif
  if (eos_henry_ierr2 /= 0) then
     call HydrateEOSGasError(natural_id,eos_henry_ierr2,hyd_auxvar,option)
  endif

  cell_pressure = max(hyd_auxvar%pres(lid),hyd_auxvar%pres(gid), &
                      hyd_auxvar%pres(spid))

  ! calculate effective porosity as a function of pressure
  if (option%iflag /= HYDRATE_UPDATE_FOR_BOUNDARY) then
    dpor_dp = 0.d0
    hyd_auxvar%effective_porosity = material_auxvar%porosity_base
    if (soil_compressibility_index > 0) then
      call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                hyd_auxvar%effective_porosity,dpor_dp)
    endif
    if (option%iflag /= HYDRATE_UPDATE_FOR_DERIVATIVE) then
      material_auxvar%porosity = hyd_auxvar%effective_porosity
    endif

    !MAN: does this need to be in this conditional?
    solid_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid)

    if (hydrate_perm_scaling) then
      call HydrateScalePermeability(hyd_auxvar%perm_base, &
                                    material_auxvar%permeability, &
                                    solid_sat_eff,option)
    endif

  endif
  if (associated(hyd_auxvar%d)) then
    hyd_auxvar%d%por_p = dpor_dp
  endif

  ! ALWAYS UPDATE THERMODYNAMIC PROPERTIES

  ! Liquid phase thermodynamic properties
  ! must use cell_pressure as the pressure, not %pres(lid)
  if (.not.option%flow%density_depends_on_salinity) then
    if (associated(hyd_auxvar%d)) then
      call EOSWaterDensity(hyd_auxvar%temp,cell_pressure, &
                           hyd_auxvar%den_kg(lid),hyd_auxvar%den(lid), &
                           hyd_auxvar%d%denl_pl,hyd_auxvar%d%denl_T,ierr)
    else
      if (hydrate_salinity) then
         ! call EOSWaterDensityExt(hyd_auxvar%temp,cell_pressure,&
         !                   hyd_auxvar%xmol(cid,lid),hyd_auxvar%den_kg(lid), &
         !                   hyd_auxvar%den(lid),ierr)
         !convert xmol to weight fraction in ppm
         aux(1) = hyd_auxvar%xmol(cid,lid)*FMWNACL / (hyd_auxvar%xmol(acid,lid) &
              * FMWAIR+ hyd_auxvar%xmol(cid,lid)*FMWNACL + &
              hyd_auxvar%xmol(wid,lid)*FMWH2O)
         call EOSWaterDensityExt(hyd_auxvar%temp,cell_pressure,aux, &
              hyd_auxvar%den_kg(lid),hyd_auxvar%den(lid),ierr)
         ! call EOSWaterCH4DuanMixture(hyd_auxvar%temp,cell_pressure, &
         !      hyd_auxvar%xmol(acid,lid),hyd_auxvar%xmol(cid,lid), &
         !      )
         ! call EOSWaterDensity(hyd_auxvar%temp,cell_pressure, &
         !      hyd_auxvar%den_kg(lid),hyd_auxvar%den(lid),ierr)
      else
      call EOSWaterDensity(hyd_auxvar%temp,cell_pressure, &
                           hyd_auxvar%den_kg(lid),hyd_auxvar%den(lid),ierr)
      endif
    endif
  else
    aux(1) = global_auxvar%m_nacl(1)
    call EOSWaterDensityExt(hyd_auxvar%temp,cell_pressure,aux, &
                              hyd_auxvar%den_kg(lid),hyd_auxvar%den(lid),ierr)
  endif

  call EOSWaterEnthalpy(hyd_auxvar%temp,cell_pressure,hw,ierr)

  hyd_auxvar%H(lid) = hw * 1.d-6 ! J/kmol -> MJ/kmol
  ! MJ/kmol comp
  hyd_auxvar%U(lid) = hyd_auxvar%H(lid) - &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        (cell_pressure / hyd_auxvar%den(lid) * &
                        1.d-6)
  if (global_auxvar%istate /= L_STATE) then
    if (global_auxvar%istate == G_STATE) then
      water_vapor_pressure = hyd_auxvar%pres(vpid)
    else
      water_vapor_pressure = hyd_auxvar%pres(spid)
    endif
    select case(hydrate_former)
      case(HYDRATE_FORMER_CO2)
        call EOSGasDensityEnergy(hyd_auxvar%temp,hyd_auxvar%pres(apid),den_air, &
             h_air,u_air,ierr)
      case(HYDRATE_FORMER_CH4)
        call EOSGasDensityEnergy(hyd_auxvar%temp,hyd_auxvar%pres(apid),den_air, &
             h_air,u_air,ierr)
        call EOSGasDensityEnergy(hyd_auxvar%temp,hyd_auxvar%pres(cid),den_co2, &
             h_co2, u_co2, ierr)
      case default
        call EOSGasDensityEnergy(hyd_auxvar%temp,hyd_auxvar%pres(apid),den_air, &
             h_air,u_air,ierr)
        call EOSGasDensityEnergy(hyd_auxvar%temp,hyd_auxvar%pres(cid),den_co2, &
             h_co2, u_co2, ierr)
    end select

    h_air = h_air * 1.d-6 ! J/kmol -> MJ/kmol
    u_air = u_air * 1.d-6 ! J/kmol -> MJ/kmol
    h_co2 = h_co2 * 1.d-6 ! J/kmol -> MJ/kmol
    u_co2 = u_co2 * 1.d-6 ! J/kmol -> MJ/kmol

    call EOSWaterSteamDensityEnthalpy(hyd_auxvar%temp,water_vapor_pressure, &
                                        den_kg_water_vapor,den_water_vapor, &
                                        h_water_vapor,ierr)
    u_water_vapor = h_water_vapor - &
                    ! Pa / kmol/m^3 = J/kmol
                    water_vapor_pressure / den_water_vapor
    u_water_vapor = (h_water_vapor - 8.314462618 *1.d3*&
                     (hyd_auxvar%temp+273.15)) *  1.d-6 ! J/kmol -> MJ/kmol
    h_water_vapor = u_water_vapor * 1.d-6 ! J/kmol -> MJ/kmol
    hyd_auxvar%den(gid) = den_water_vapor + den_air + den_co2
    hyd_auxvar%den_kg(gid) = den_kg_water_vapor + den_air*fmw_comp(gid)&
                             + den_co2 * fmw_comp(cid) !DF: Why are these summed?
    xmol_air_in_gas = hyd_auxvar%xmol(acid,gid)
    xmol_water_in_gas = hyd_auxvar%xmol(wid,gid)
    xmol_co2_in_gas = hyd_auxvar%xmol(cid,gid)
    ! MJ/kmol
    hyd_auxvar%U(gid) = xmol_water_in_gas * u_water_vapor + &
                            xmol_air_in_gas * u_air + &
                            xmol_co2_in_gas * u_co2
    Hg_mixture_fractioned = xmol_water_in_gas*h_water_vapor + &
                            xmol_air_in_gas*h_air + &
                            xmol_co2_in_gas*h_co2
    hyd_auxvar%H(gid) = hyd_auxvar%U(gid) + &
                        ! Pa / kmol/m^3 * 1.e-6 = MJ/kmol
                        hyd_auxvar%pres(gid)/hyd_auxvar%den(gid) * 1.d-6

  endif

  if (hydrate_eff_sat_scaling) then
    if (hyd_auxvar%sat(lid) > 0.d0 .or. hyd_auxvar%sat(gid) > 0.d0) then
      liq_sat_eff = hyd_auxvar%sat(lid) / (hyd_auxvar%sat(lid)+hyd_auxvar%sat(gid))
      gas_sat_eff = 1.d0 - liq_sat_eff
    else
      liq_sat_eff = 0.d0
      gas_sat_eff = 0.d0
    endif
  else
    liq_sat_eff = hyd_auxvar%sat(lid)
    gas_sat_eff = hyd_auxvar%sat(gid)
  endif

  if (liq_sat_eff > 0.d0) then
    if (liq_sat_eff >= 1.d0) then
      krl = 1.d0
    else
      call characteristic_curves%liq_rel_perm_function% &
           RelativePermeability(liq_sat_eff,krl,dkrl_dsatl,option)
      krl = max(0.d0,krl)
    endif
    call EOSWaterViscosity(hyd_auxvar%temp,cell_pressure, &
                               hyd_auxvar%pres(spid),visl,ierr)
    hyd_auxvar%mobility(lid) = krl/visl
    hyd_auxvar%kr(lid) = krl
  endif

  if (gas_sat_eff > 0.d0) then
    if (gas_sat_eff >=1.d0) then
      krg = 1.d0
    else
      call characteristic_curves%gas_rel_perm_function% &
           RelativePermeability(liq_sat_eff,krg,dkrg_dsatl,option)
      krg = max(0.d0,krg)
    endif
    select case(hydrate_former)
      CASE(HYDRATE_FORMER_CO2)!DF: for now, keep as methane
        !call VISCO2(hyd_auxvar%temp,hyd_auxvar%den_kg(gid), visg)
        call EOSGasViscosity(hyd_auxvar%temp,hyd_auxvar%pres(apid), &
                              hyd_auxvar%pres(gid),den_air,visg,ierr)
      CASE(HYDRATE_FORMER_CH4)
        call EOSGasViscosity(hyd_auxvar%temp,hyd_auxvar%pres(apid), &
                              hyd_auxvar%pres(gid),den_air,visg,ierr)
      case default
        call EOSGasViscosity(hyd_auxvar%temp,hyd_auxvar%pres(apid), &
                              hyd_auxvar%pres(gid),den_air,visg,ierr)
    END SELECT
    hyd_auxvar%mobility(gid) = krg/visg
    hyd_auxvar%kr(gid) = krg
  endif
  if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
    call EOSHydrateEnthalpy(hyd_auxvar%temp, H_hyd, HYDRATE_FORMER_CH4)
    hyd_auxvar%U(hid) = H_hyd !- cell_pressure/hyd_auxvar%den(hid)*1.d-6
    hyd_auxvar%H(hid) = H_hyd
    hyd_auxvar%mobility(hid) = 0.d0

    call EOSHydrateEnthalpy(hyd_auxvar%temp, H2_hyd, HYDRATE_FORMER_CH4)
    hyd_auxvar%U(h2id) = H2_hyd
    hyd_auxvar%H(h2id) = H2_hyd
    hyd_auxvar%mobility(h2id) = 0.d0
  else
    call EOSHydrateEnthalpy(hyd_auxvar%temp, H_hyd, HYDRATE_FORMER_CH4)
    hyd_auxvar%U(hid) = H_hyd !- cell_pressure/hyd_auxvar%den(hid)*1.d-6
    hyd_auxvar%H(hid) = H_hyd
    hyd_auxvar%mobility(hid) = 0.d0

    call EOSHydrateEnthalpy(hyd_auxvar%temp, H2_hyd, HYDRATE_FORMER_CH4)
    hyd_auxvar%U(h2id) = H2_hyd
    hyd_auxvar%H(h2id) = H2_hyd
    hyd_auxvar%mobility(h2id) = 0.d0
  endif

  call EOSIceEnergy(hyd_auxvar%temp, U_ice)
  hyd_auxvar%xmol(wid,iid) = 1.d0
  hyd_auxvar%xmol(gid,iid) = 0.d0
  !call EOSWaterDensityIcePainter(hyd_auxvar%temp,hyd_auxvar%pres(lid), &
  !                  PETSC_FALSE, hyd_auxvar%den(iid), &
  !                  dden_ice_dT, dden_ice_dP, ierr)
  hyd_auxvar%den(iid) = ICE_DENSITY
  hyd_auxvar%den_kg(iid) = ICE_DENSITY_KG
  hyd_auxvar%U(iid) = U_ice
  hyd_auxvar%H(iid) = U_ice
  hyd_auxvar%mobility(iid) = 0.d0

  hyd_auxvar%srl = characteristic_curves%gas_rel_perm_function%sr
  hyd_auxvar%srg = characteristic_curves%gas_rel_perm_function%srg

  !hyd_auxvar%h1_h2_ratio = hyd_auxvar%xmol(acid,lid) / &
  !                     (hyd_auxvar%xmol(acid,lid) + hyd_auxvar%xmol(cid,lid))
   ! call HydratePrintAuxVars(hyd_auxvar,global_auxvar,material_auxvar, &
  !      natural_id,'after_auxvars: ',option)

end subroutine HydrateAuxVarCompute4

! ************************************************************************** !

subroutine HydrateComputeStateSpecificVars4(x,state,hyd_auxvar, &
                  characteristic_curves,material_auxvar,natural_id,option)
  !
  ! Computes auxvars depending on phase state for 4 equations.
  !
  ! Author: Michael Nole
  ! Date: 08/26/20
  !

  use Option_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use EOS_Water_module
  use EOS_Gas_module
  use Characteristic_Curves_module
  use Material_Aux_class
  use co2eos_module
  use co2_span_wagner_module
  use co2_span_wagner_spline_module

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  PetscInt :: state
  type(hydrate_auxvar_type) :: hyd_auxvar
  class(characteristic_curves_type) :: characteristic_curves
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id

  PetscReal :: PE_hyd, PE_hyd2, dP, dP2, K_H_tilde, K_H_tilde2, K_H_tilde_hyd, &
               K_H_tilde_hyd2, sigma
  PetscInt :: gid, lid, acid, wid, eid, hid, iid, h2id, cid, &
              apid, cpid, vpid, spid
  PetscReal :: dTf, h_sat_eff, i_sat_eff, liq_sat_eff, gas_sat_eff
  PetscReal :: dpc_dsatl
  PetscInt :: ierr

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  h2id = option%hydrate2_phase ! co2 hydrate

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id
  wid = option%water_id
  cid = option%co2_id
  eid = option%energy_id

  PE_hyd = 0.d0
  PE_hyd2 = 0.d0

  select case(state)
    case(L_STATE)
!     ********* Aqueous State (A) ********************************
!     Primary variables: Pa, Xma, Xma2, T
!

      hyd_auxvar%pres(lid) = x(HYDRATE_LIQUID_PRESSURE_DOF)
      hyd_auxvar%xmol(acid,lid) = x(HYDRATE_L_STATE_X_MOLE_DOF)
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)


      hyd_auxvar%xmol(acid,lid) = max(0.d0,hyd_auxvar%xmol(acid,lid))
      hyd_auxvar%xmol(cid,lid) = max(0.d0,hyd_auxvar%xmol(cid,lid))

      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid) &
                                 - hyd_auxvar%xmol(cid,lid)
      hyd_auxvar%xmol(wid,gid) = 0.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0
      hyd_auxvar%xmol(cid,gid) = 0.d0
      hyd_auxvar%sat(lid) = 1.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0
      hyd_auxvar%sat(h2id) = 0.d0

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                        hyd_auxvar%pres(spid),ierr)
      call HydratePE(hyd_auxvar%temp, 0.d0,hyd_auxvar%xmol(cid,lid),&
                     PE_hyd, dP, characteristic_curves, &
                     material_auxvar,option, HYDRATE_FORMER_CH4)
      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)
      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CO2,K_H_tilde2)

      K_H_tilde_hyd = K_H_tilde

      if (hydrate_adjust_ghsz_solubility) then
        call HydrateGHSZSolubilityCorrection(hyd_auxvar%temp,hyd_auxvar% &
                                           pres(lid),dP,K_H_tilde_hyd)
      endif

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(gid) = max(hyd_auxvar%pres(lid),hyd_auxvar%pres(spid))
      hyd_auxvar%pres(cid) = K_H_tilde2*hyd_auxvar%xmol(cid,lid)
      hyd_auxvar%pres(apid) = K_H_tilde_hyd*hyd_auxvar%xmol(acid,lid)

      if (hyd_auxvar%pres(gid) <= 0.d0) then
        write(option%io_buffer,'(''Negative gas pressure at cell '', &
          & i8,'' in HydrateAuxVarCompute(L_STATE).  Attempting bailout.'')') &
          natural_id
        call PrintMsgByRank(option)
        hyd_auxvar%pres(vpid) = 0.5d0*hyd_auxvar%pres(spid)
        hyd_auxvar%pres(gid) = hyd_auxvar%pres(vpid) + hyd_auxvar%pres(apid)
      else
        hyd_auxvar%pres(vpid) = hyd_auxvar%pres(lid) - hyd_auxvar%pres(apid) - &
                                hyd_auxvar%pres(cid)
      endif
      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%h1_h2_ratio = hyd_auxvar%xmol(acid,lid) / &
                       (hyd_auxvar%xmol(acid,lid) + hyd_auxvar%xmol(cid,lid))
    case (G_STATE)
!     ********* Gas State (G) ********************************
!     Primary variables: Pg, Pa, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%pres(apid) = x(HYDRATE_G_STATE_AIR_PRESSURE_DOF)
      hyd_auxvar%pres(cid) = x(HYDRATE_G_STATE_CO2_PRESSURE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 1.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0
      hyd_auxvar%sat(h2id) = 0.d0

      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / &
                                    hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(cid,gid) = hyd_auxvar%pres(cid) / &
                                   hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid) - &
                                   hyd_auxvar%xmol(cid,gid)
      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                    hyd_auxvar%pres(spid),ierr)

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 0.d0

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call characteristic_curves%saturation_function% &
             CapillaryPressure(hyd_auxvar%sat(lid), &
                               hyd_auxvar%pres(cpid),dpc_dsatl,option)
      endif

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - &
                             hyd_auxvar%pres(cpid)

    case (H_STATE)
!     ********* Hydrate State (H) ********************************
!     Primary variables: Pg, Xmh, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      x(HYDRATE_GAS_SATURATION_DOF) = MOL_RATIO_FORMER
      hyd_auxvar%xmol(acid,hid) = x(HYDRATE_GAS_SATURATION_DOF)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = 0.d0
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 1.d0
      hyd_auxvar%sat(iid) = 0.d0
      hyd_auxvar%sat(h2id) = 0.d0

      call HydratePE(hyd_auxvar%temp,hyd_auxvar%sat(hid),&
                     hyd_auxvar%xmol(cid,lid),PE_hyd,dP, &
                     characteristic_curves, material_auxvar, option, HYDRATE_FORMER_CH4)
      hyd_auxvar%pres(apid) = PE_hyd


      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

   case (H2_STATE)
!     ********* CO2 Hydrate State (H2) *************************
!     Primary variables: Pg, Xmh, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      x(HYDRATE_GAS_SATURATION_DOF) = MOL_RATIO_FORMER
      hyd_auxvar%xmol(cid,hid) = x(HYDRATE_GAS_SATURATION_DOF)
      X(HYDRATE_L_STATE_X2_MOLE_DOF) = 0.d0
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0
      hyd_auxvar%sat(h2id) = 1.d0

      call HydratePE(hyd_auxvar%temp,hyd_auxvar%sat(h2id),&
                     hyd_auxvar%xmol(cid,lid),PE_hyd,dP, &
                     characteristic_curves, material_auxvar,&
                     option,HYDRATE_FORMER_CO2)
      hyd_auxvar%pres(apid) = PE_hyd

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

   case (HH_STATE)
!     ********* CH4 and CO2 Hydrate State (HH) *************************
!     Primary variables: Pg, Xmh, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      x(HYDRATE_GAS_SATURATION_DOF) = MOL_RATIO_FORMER
      hyd_auxvar%xmol(acid,hid) = x(HYDRATE_GAS_SATURATION_DOF)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = 0.d0
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = max(0.d0,min(1.d0,hyd_auxvar%sat(hid)))
      hyd_auxvar%sat(iid) = 0.d0
      hyd_auxvar%sat(h2id) = 1.d0 - hyd_auxvar%sat(hid)

      call HydratePE(hyd_auxvar%temp,hyd_auxvar%sat(h2id),0.d0,PE_hyd,dP, &
           characteristic_curves, material_auxvar, option,HYDRATE_FORMER_CH4)
      hyd_auxvar%pres(apid) = PE_hyd

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid) - &
                                   hyd_auxvar%xmol(cid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

    case(I_STATE)
!     ********* Ice State (I) ********************************
!     Primary variables: Pg, Xmi, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      x(HYDRATE_GAS_SATURATION_DOF) = 0.d0
      x(HYDRATE_HYD2_SATURATION_DOF) =  0.d0
      hyd_auxvar%sat(h2id) = x(HYDRATE_HYD2_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 1.d0
      hyd_auxvar%sat(h2id) = 0.d0

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(apid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%xmol(acid,lid) = 0.d0
      hyd_auxvar%xmol(wid,lid) = 1.d0
      hyd_auxvar%xmol(wid,gid) = 1.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0

    case(GA_STATE)
!     ********* Gas & Aqueous State (GA) ********************************
!     Primary variables: Pg, Sg, T, Xma2
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(gid) = max(0.d0,min(1.d0,hyd_auxvar%sat(gid)))
      hyd_auxvar%sat(lid) = 1.d0 - hyd_auxvar%sat(gid)
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0
      hyd_auxvar%sat(h2id) = 0.d0

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)
      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(vpid)

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)

      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call characteristic_curves%saturation_function% &
             CapillaryPressure(hyd_auxvar%sat(lid), hyd_auxvar%pres(cpid), &
                               dpc_dsatl,option)
      endif

      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(hyd_auxvar%temp,sigma)
      endif
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid)*sigma

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
        hyd_auxvar%pres(cid) = K_H_tilde * hyd_auxvar%xmol(cid,lid)
        hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(vpid) -&
                                hyd_auxvar%pres(cid)
        hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      endif

      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid) - &
                                 hyd_auxvar%xmol(cid,lid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(cid,gid) = hyd_auxvar%pres(cid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid) - &
                                 hyd_auxvar%xmol(cid,gid)

    case(HG_STATE)
!     ********* Hydrate & Gas State (HG) ********************************
!     Primary variables: Pg, Sg, T, Pco2
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%pres(cid) = x(HYDRATE_G_STATE_CO2_PRESSURE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(hid) = 1.d0 - hyd_auxvar%sat(gid)
      hyd_auxvar%sat(iid) = 0.d0
      hyd_auxvar%sat(h2id) = 0.d0

      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
        call HydratePE(hyd_auxvar%temp, hyd_auxvar%sat(hid)+ &
              hyd_auxvar%sat(gid),0.d0,PE_hyd, dP, characteristic_curves, &
              material_auxvar, option,HYDRATE_FORMER_CH4)
      else
        call HydratePE(hyd_auxvar%temp, 2.d0 * hyd_auxvar%sat(hid),0.d0, &
              PE_hyd, dP, characteristic_curves, material_auxvar, option, &
              HYDRATE_FORMER_CH4)
      endif

      hyd_auxvar%pres(apid) = PE_hyd

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)
      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CO2,K_H_tilde2)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid) - &
                                hyd_auxvar%pres(cid)
      !DF check vpid
      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(cid,lid) = hyd_auxvar%pres(cid) / K_H_tilde2
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid) - &
                                   hyd_auxvar%xmol(cid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(acid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(cid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid) - &
                                    hyd_auxvar%xmol(acid,gid)

    case(H2G_STATE)
!     ********* CO2 Hydrate & Gas State (H2G) *******************************
!     Primary variables: Pg, Sg, Pco2, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%pres(cid) = x(HYDRATE_G_STATE_CO2_PRESSURE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(h2id) = 1.d0 - hyd_auxvar%sat(gid)
      hyd_auxvar%sat(iid) = 0.d0
      hyd_auxvar%sat(hid) = 0.d0

      if (hyd_auxvar%sat(h2id) > hyd_auxvar%sat(gid)) then
        call HydratePE(hyd_auxvar%temp, hyd_auxvar%sat(h2id)+ &
                hyd_auxvar%sat(gid),0.d0,PE_hyd, dP, characteristic_curves, &
                material_auxvar, option,HYDRATE_FORMER_CO2)
      else
        call HydratePE(hyd_auxvar%temp, 2.d0 * hyd_auxvar%sat(h2id),0.d0, &
                PE_hyd, dP, characteristic_curves, material_auxvar, option, &
                HYDRATE_FORMER_CO2)
      endif

      hyd_auxvar%pres(apid) = PE_hyd

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)
      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CO2,K_H_tilde2)

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid) - &
                                hyd_auxvar%pres(cid)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(cid,lid) = hyd_auxvar%pres(cid) / K_H_tilde2
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(cid,gid) = hyd_auxvar%pres(cid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid) - &
                                    hyd_auxvar%xmol(cid,gid)

    case(HA_STATE)
!     ********* Hydrate & Aqueous State (HA) ********************************
!     Primary variables: Pg, Sh, T, XMA2

      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)

      hyd_auxvar%sat(hid) = max(0.d0,min(1.d0,hyd_auxvar%sat(hid)))
      hyd_auxvar%sat(lid) = 1.d0 - hyd_auxvar%sat(hid)
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0
      hyd_auxvar%sat(h2id) = 0.d0

      call HydratePE(hyd_auxvar%temp,hyd_auxvar%sat(hid),&
                     hyd_auxvar%xmol(cid,lid), PE_hyd, dP,&
              characteristic_curves, material_auxvar, option,HYDRATE_FORMER_CH4)

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)
      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CO2,K_H_tilde2)
      if (hydrate_adjust_ghsz_solubility) then
        call HydrateGHSZSolubilityCorrection(hyd_auxvar%temp,hyd_auxvar% &
                                           pres(gid),dP,K_H_tilde)
      endif

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)

      hyd_auxvar%pres(cpid) = 0.d0

      hyd_auxvar%pres(apid) = PE_hyd!hyd_auxvar%pres(gid) !- hyd_auxvar%pres(spid)
      hyd_auxvar%pres(cid) = hyd_auxvar%xmol(cid,lid)*K_H_tilde2
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid)


      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
        hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      else
        hyd_auxvar%xmol(acid,lid) = PE_hyd / K_H_tilde!hyd_auxvar%pres(apid) / K_H_tilde
      endif

      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid) - &
                                 hyd_auxvar%xmol(cid,lid)
      hyd_auxvar%xmol(wid,gid) = 0.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid) - hyd_auxvar%pres(cid)
      hyd_auxvar%h1_h2_ratio = hyd_auxvar%xmol(acid,lid)/ &
                          (hyd_auxvar%xmol(acid,lid) + hyd_auxvar%xmol(cid,lid))

    case(H2A_STATE)
!     ********* CO2 Hydrate & Aqueous State (H2A) ********************************
!     Primary variables: Pg, XMA, T, SH2
!break point

      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%xmol(acid,lid) = x(HYDRATE_L_STATE_X_MOLE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)
      hyd_auxvar%sat(h2id) = x(HYDRATE_HYD2_SATURATION_DOF)

      hyd_auxvar%sat(h2id) = max(0.d0,min(1.d0,hyd_auxvar%sat(h2id)))

      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(lid) = 1.d0 - hyd_auxvar%sat(h2id)
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0

      call HydratePE(hyd_auxvar%temp,hyd_auxvar%sat(h2id),&
                     hyd_auxvar%xmol(cid,lid),PE_hyd,&
                     dP,characteristic_curves, material_auxvar, option, &
                     HYDRATE_FORMER_CO2)
      hyd_auxvar%pres(apid) = PE_hyd

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)
      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CO2,K_H_tilde2)

      if (HYDRATE_ADJUST_GHSZ_SOLUBILITY) then
        call HydrateGHSZSolubilityCorrection(hyd_auxvar%temp,hyd_auxvar% &
                                           pres(gid),dP,K_H_tilde2)
      endif

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)


      hyd_auxvar%pres(cpid) = 0.d0

      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid)
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
        hyd_auxvar%xmol(cid,lid) = hyd_auxvar%pres(apid) / K_H_tilde2
      else
        hyd_auxvar%xmol(cid,lid) = hyd_auxvar%pres(cid) / K_H_tilde2
      endif

      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(cid,lid) - &
                                 hyd_auxvar%xmol(acid,lid)

      hyd_auxvar%xmol(wid,gid) = 0.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0

      hyd_auxvar%h1_h2_ratio = hyd_auxvar%xmol(acid,lid)/ &
                          (hyd_auxvar%xmol(acid,lid) + hyd_auxvar%xmol(cid,lid))


    case(HI_STATE)
!     ********* Hydrate & Ice State (HI) ********************************
!     Primary variables: Pg, Sh, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = 0.d0
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(iid) = 1.d0 - hyd_auxvar%sat(hid)
      hyd_auxvar%sat(h2id) = 0.d0

      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
        call HydratePE(hyd_auxvar%temp, hyd_auxvar%sat(hid)+ &
              hyd_auxvar%sat(iid),hyd_auxvar%xmol(cid,lid),PE_hyd,dP,&
              characteristic_curves, &
              material_auxvar,option,HYDRATE_FORMER_CH4)
      else
        call HydratePE(hyd_auxvar%temp, 2.d0 * hyd_auxvar%sat(hid),&
              hyd_auxvar%xmol(cid,lid),PE_hyd, &
              dP, characteristic_curves, material_auxvar,&
              option,HYDRATE_FORMER_CH4)
      endif

      hyd_auxvar%pres(apid) = PE_hyd

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)
      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CO2,K_H_tilde2)

      hyd_auxvar%pres(cpid) = 0.d0
!      hyd_auxvar%pres(cid) = 
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

    case(GI_STATE)
!     ********* Gas & Ice State (GI) ********************************
!     Primary variables: Pg, Si, T
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%pres(cid) = x(HYDRATE_G_STATE_CO2_PRESSURE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(gid) = 1.d0 - hyd_auxvar%sat(iid)
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(h2id) = 0.d0

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(vpid) - &
                                hyd_auxvar%pres(cid)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)

      hyd_auxvar%xmol(acid,lid) = 0.d0
      hyd_auxvar%xmol(wid,lid) = 1.d0
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(cid,gid) = hyd_auxvar%pres(cid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

    case(AI_STATE)
!     ********* Aqueous & Ice State (AI) ********************************
!     Primary variables: Pl, Xma, Sl
!
      hyd_auxvar%pres(lid) = x(HYDRATE_LIQUID_PRESSURE_DOF)
      hyd_auxvar%xmol(acid,lid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%sat(lid) = x(HYDRATE_ENERGY_DOF)
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)

      hyd_auxvar%xmol(acid,lid) = max(0.d0,hyd_auxvar%xmol(acid,lid))

      if (hyd_auxvar%istatechng) then
        hyd_auxvar%sat(lid) = max(0.d0,min(1.d0,hyd_auxvar%sat(lid)))
      endif

      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(iid) = 1.d0 - hyd_auxvar%sat(lid)
      hyd_auxvar%sat(h2id) = 0.d0

      if (hydrate_with_gibbs_thomson) then
        call GibbsThomsonFreezing(hyd_auxvar%sat(lid),6017.1d0,ICE_DENSITY,&
                                TQD, dTf,characteristic_curves, &
                                material_auxvar,option)
      else
        dTf = 0.d0
      endif

      hyd_auxvar%temp = TQD+dTf

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)

      !hyd_auxvar%pres(spid) = 1.d-6
      hyd_auxvar%pres(gid) = max(hyd_auxvar%pres(lid),hyd_auxvar%pres(spid))
      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(apid) = K_H_tilde*hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%pres(vpid) = 0.d0 !hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)

!      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
         hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde - &
              hyd_auxvar%xmol(cid,lid)
         hyd_auxvar%xmol(acid,lid) = max(0.d0,hyd_auxvar%xmol(acid,lid))
      else
         hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      endif
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid) - &
                                   hyd_auxvar%xmol(cid,lid)
      hyd_auxvar%xmol(acid,gid) = 0.d0
      hyd_auxvar%xmol(wid,gid) = 0.d0

    case(HGA_STATE)
!     ********* Hydrate, Gas, & Aqueous State (HGA) **************************
!     Primary variables: Sl, Sh, T, XMA2
!
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)

      !hyd_auxvar%sat(lid) = max(0.d0,min(1.d0,hyd_auxvar%sat(lid)))
      !hyd_auxvar%sat(hid) = max(0.d0,min(1.d0,hyd_auxvar%sat(hid)))

      !if (hyd_auxvar%sat(lid) + hyd_auxvar%sat(hid) > 1.d0) then
      ! hyd_auxvar%sat(hid) = 1.d0 - hyd_auxvar%sat(lid)
      ! x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
      !endif

      hyd_auxvar%sat(gid) = max(0.d0,1.d0 - hyd_auxvar%sat(lid) - &
                                hyd_auxvar%sat(hid))
      hyd_auxvar%sat(iid) = 0.d0
      hyd_auxvar%sat(h2id) = 0.d0

      if (hydrate_gt_3phase) then
        if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
          h_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(gid)
          gas_sat_eff = 2.d0 * hyd_auxvar%sat(gid)
        else
          gas_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(gid)
          h_sat_eff = 2.d0 * hyd_auxvar%sat(hid)
        endif
      endif

      h_sat_eff = hyd_auxvar%sat(hid)

      if (hydrate_eff_sat_scaling) then
        liq_sat_eff = hyd_auxvar%sat(lid)/(hyd_auxvar%sat(lid)+ &
                    hyd_auxvar%sat(gid))
      else
        liq_sat_eff = hyd_auxvar%sat(lid)
      endif
      call HydratePE(hyd_auxvar%temp, h_sat_eff,hyd_auxvar%xmol(cid,lid),&
           PE_hyd, dP,&
           characteristic_curves, material_auxvar,option, HYDRATE_FORMER_CH4)
      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call characteristic_curves%saturation_function%CapillaryPressure( &
                liq_sat_eff, hyd_auxvar%pres(cpid), &
                dpc_dsatl,option)
      endif

      hyd_auxvar%pres(apid) = PE_hyd

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                    hyd_auxvar%pres(spid),ierr)

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(gid) = hyd_auxvar%pres(apid) + hyd_auxvar%pres(vpid)

      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(hyd_auxvar%temp,sigma)
      endif
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid)*sigma

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
        hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde - &
        hyd_auxvar%xmol(cid,lid)
        hyd_auxvar%xmol(acid,lid) = max(0.d0,hyd_auxvar%xmol(acid,lid))
      else
        hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      endif
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid) - &
                                 hyd_auxvar%xmol(cid,lid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

    case(HAI_STATE)
!     ********* Hydrate, Aqueous, & Ice State (HAI) **************************
!     Primary variables: Pg, Sl, Si
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(hid) = 1.d0 - hyd_auxvar%sat(lid) - hyd_auxvar%sat(iid)
      hyd_auxvar%sat(h2id) = 0.d0

      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
        h_sat_eff = hyd_auxvar%sat(hid)+hyd_auxvar%sat(iid)
        i_sat_eff = 2.d0 * hyd_auxvar%sat(iid)
      else
        h_sat_eff = 2.d0 * hyd_auxvar%sat(hid)
        i_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid)
      endif

      if (hydrate_with_gibbs_thomson) then
        call GibbsThomsonFreezing(1.d0-i_sat_eff,6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,material_auxvar,option)
      else
        dTf = 0.d0
      endif

      hyd_auxvar%temp = TQD+dTf
      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)
      call HydratePE(hyd_auxvar%temp,h_sat_eff,hyd_auxvar%xmol(cid,lid),PE_hyd,&
           dP, characteristic_curves, material_auxvar,option,HYDRATE_FORMER_CH4)
      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)

      hyd_auxvar%pres(cpid) = 0.d0

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(apid) = PE_hyd !hyd_auxvar%pres(gid)
      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)
      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
         hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde - &
              hyd_auxvar%xmol(cid,lid)
         hyd_auxvar%xmol(acid,lid) = max(0.d0,hyd_auxvar%xmol(acid,lid))
      else
         hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      endif

!      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = 0.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0

    case(HGI_STATE)
!     ********* Hydrate, Gas, & Ice State (HGI) ******************************
!     Primary variables: Sh, Si, T
!
      hyd_auxvar%sat(iid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%pres(cid) = x(HYDRATE_G_STATE_CO2_PRESSURE_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(lid) = 0.d0
      hyd_auxvar%sat(h2id) = 0.d0
      hyd_auxvar%sat(gid) = 1.d0 - hyd_auxvar%sat(hid) - hyd_auxvar%sat(iid)

      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
        if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
          call HydratePE(hyd_auxvar%temp, 1.d0,0.d0, PE_hyd, dP, &
                characteristic_curves, material_auxvar, option,HYDRATE_FORMER_CH4)
        else
          call HydratePE(hyd_auxvar%temp, 3.d0 * hyd_auxvar%sat(iid) + &
                2.d0 * (hyd_auxvar%sat(hid)-hyd_auxvar%sat(iid)),0.d0,PE_hyd, &
                dP, characteristic_curves, material_auxvar, option, &
                HYDRATE_FORMER_CH4)
        endif
      elseif (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
        call HydratePE(hyd_auxvar%temp, 3.d0 * hyd_auxvar%sat(gid) + &
               2.d0 * (hyd_auxvar%sat(hid) - hyd_auxvar%sat(gid)), 0.d0,PE_hyd, &
               dP, characteristic_curves, material_auxvar, option,&
               HYDRATE_FORMER_CH4)
      else
        call HydratePE(hyd_auxvar%temp, 3.d0 * hyd_auxvar%sat(hid),0.d0,PE_hyd,&
              dP,characteristic_curves,material_auxvar, option,HYDRATE_FORMER_CH4)
      endif

      hyd_auxvar%pres(apid) = PE_hyd

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                    hyd_auxvar%pres(spid),ierr)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(gid) = hyd_auxvar%pres(apid) + hyd_auxvar%pres(vpid) + &
                               hyd_auxvar%pres(cid)

      hyd_auxvar%xmol(acid,lid) = 0.d0
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(wid,gid) = hyd_auxvar%pres(vpid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(cid,gid) = hyd_auxvar%pres(cid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(acid,gid) = 1.d0 - hyd_auxvar%xmol(wid,gid)

    case(GAI_STATE)
!     ********* Gas, Aqueous, & Ice State (GAI) ******************************
!     Primary variables: Pg, Sl, Si
!
      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_ENERGY_DOF)


      hyd_auxvar%sat(gid) = 1.d0 - hyd_auxvar%sat(lid) - hyd_auxvar%sat(iid)
      hyd_auxvar%sat(hid) = 0.d0
      hyd_auxvar%sat(h2id) = 0.d0

      call GibbsThomsonFreezing(1.d0-hyd_auxvar%sat(iid),6017.1d0, &
              ICE_DENSITY,TQD,dTf,characteristic_curves, material_auxvar,option)

      hyd_auxvar%temp = TQD+dTf

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                    hyd_auxvar%pres(spid),ierr)
      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(vpid)

      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call characteristic_curves%saturation_function% &
             CapillaryPressure(hyd_auxvar%sat(lid), &
                             hyd_auxvar%pres(cpid),dpc_dsatl,option)
      endif

      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(hyd_auxvar%temp,sigma)
      endif
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid)*sigma

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
         hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde - &
              hyd_auxvar%xmol(cid,lid)
         hyd_auxvar%xmol(acid,lid) = max(0.d0,hyd_auxvar%xmol(acid,lid))
      else
         hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      endif

!      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

    case(QUAD_STATE)
!     ********* Quadruple Point (HGAI) ********************************
!     Primary variables: Sg, Sl, Si
!
      hyd_auxvar%sat(lid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(gid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%xmol(cid,lid) = x(HYDRATE_L_STATE_X2_MOLE_DOF)
      hyd_auxvar%sat(iid) = x(HYDRATE_ENERGY_DOF)

      if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
        h_sat_eff = hyd_auxvar%sat(hid)+hyd_auxvar%sat(iid)
        i_sat_eff = 2.d0 * hyd_auxvar%sat(iid)
      else
        h_sat_eff = 2.d0 * hyd_auxvar%sat(hid)
        i_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid)
      endif

      if (hydrate_with_gibbs_thomson) then
        call GibbsThomsonFreezing(1.d0-i_sat_eff,6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,material_auxvar,option)
      else
        dTf = 0.d0
      endif

      !if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
      !  if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
      !    call HydratePE(hyd_auxvar%temp, 1.d0 - hyd_auxvar%sat(lid), &
      !            PE_hyd, characteristic_curves, option)
      !  else
      !    call HydratePE(hyd_auxvar%temp, 3.d0 * hyd_auxvar%sat(iid) + &
      !            2.d0 * (hyd_auxvar%sat(hid) - hyd_auxvar%sat(iid)), PE_hyd, &
      !            characteristic_curves, option)
      !    call GibbsThomsonFreezing(3.d0 * hyd_auxvar%sat(iid), 6017.1d0, &
      !            ICE_DENSITY, TQD, dTf, characteristic_curves, option)
      !  endif
      !elseif (hyd_auxvar%sat(hid) > hyd_auxvar%sat(gid)) then
      !  call HydratePE(hyd_auxvar%temp, 3.d0 * hyd_auxvar%sat(gid) + &
      !          2.d0 * (hyd_auxvar%sat(hid) - hyd_auxvar%sat(gid)), &
      !          PE_hyd, characteristic_curves, option)
      !  call GibbsThomsonFreezing(1.d0-hyd_auxvar%sat(lid), 6017.1d0, &
      !          ICE_DENSITY, TQD, dTf, characteristic_curves, option)
      !else
      !  call HydratePE(hyd_auxvar%temp, 3.d0 *hyd_auxvar%sat(hid), &
      !          PE_hyd, characteristic_curves, option)
      !  if (hyd_auxvar%sat(iid) < hyd_auxvar%sat(gid)) then
      !    call GibbsThomsonFreezing(3.d0 * hyd_auxvar%sat(hid) + &
      !            2.d0 * (hyd_auxvar%sat(iid) - hyd_auxvar%sat(hid)), &
      !            6017.1d0, ICE_DENSITY, TQD, dTf, characteristic_curves, &
      !            option)
      !  endif      
      !endif

      !hyd_auxvar%temp = TQD + dTf
      call HydratePE(hyd_auxvar%temp,h_sat_eff,hyd_auxvar%xmol(cid,lid),PE_hyd, dP, &
           characteristic_curves, material_auxvar, option,HYDRATE_FORMER_CH4)
      hyd_auxvar%pres(apid) = PE_hyd

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                    hyd_auxvar%pres(spid),ierr)

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)

      !hyd_auxvar%pres(spid) = 1.d-6

      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(spid)
      hyd_auxvar%pres(gid) = hyd_auxvar%pres(apid) + hyd_auxvar%pres(vpid) + &
                             hyd_auxvar%pres(cid)

      if (hydrate_no_pc) then
        hyd_auxvar%pres(cpid) = 0.d0
      else
        call characteristic_curves%saturation_function% &
             CapillaryPressure(hyd_auxvar%sat(lid), &
                               hyd_auxvar%pres(cpid),dpc_dsatl,option)
      endif
      !IFT calculation
      sigma=1.d0
      if (characteristic_curves%saturation_function%calc_int_tension) then
       call characteristic_curves%saturation_function% &
           CalcInterfacialTension(hyd_auxvar%temp,sigma)
      endif
      hyd_auxvar%pres(cpid) = hyd_auxvar%pres(cpid)*sigma

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(cpid)
      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
         hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde - &
              hyd_auxvar%xmol(cid,lid)
         hyd_auxvar%xmol(acid,lid) = max(0.d0,hyd_auxvar%xmol(acid,lid))
      else
         hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      endif

!      hyd_auxvar%xmol(acid,lid) = hyd_auxvar%pres(apid) / K_H_tilde
      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid) - &
                                   hyd_auxvar%xmol(cid,lid)
      hyd_auxvar%xmol(acid,gid) = hyd_auxvar%pres(apid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(cid,gid) = hyd_auxvar%pres(cid) / hyd_auxvar%pres(gid)
      hyd_auxvar%xmol(wid,gid) = 1.d0 - hyd_auxvar%xmol(acid,gid)

! States specific to hydrate_num_eqns = 4:

    case (HHA_STATE)
!     ********* H1, H2, & Aqueous State (HHA) ********************************
!     Primary variables: Pg, Sh, T, SH2

      hyd_auxvar%pres(gid) = x(HYDRATE_GAS_PRESSURE_DOF)
      hyd_auxvar%sat(hid) = x(HYDRATE_GAS_SATURATION_DOF)
      hyd_auxvar%sat(h2id) = x(HYDRATE_HYD2_SATURATION_DOF)
      hyd_auxvar%temp = x(HYDRATE_ENERGY_DOF)

      hyd_auxvar%sat(hid) = max(0.d0,min(1.d0,hyd_auxvar%sat(hid)))
      hyd_auxvar%sat(h2id) = max(0.d0,min(1.d0,hyd_auxvar%sat(h2id)))

      if (hyd_auxvar%sat(hid) < 1.d-10) hyd_auxvar%sat(hid) = 0.d0
      if (hyd_auxvar%sat(h2id) < 1.d-10) hyd_auxvar%sat(h2id) = 0.d0

      hyd_auxvar%sat(lid) = 1.d0 - hyd_auxvar%sat(hid) - hyd_auxvar%sat(h2id)
      hyd_auxvar%sat(gid) = 0.d0
      hyd_auxvar%sat(iid) = 0.d0

      call HydratePE(hyd_auxvar%temp,0.d0,hyd_auxvar%xmol(cid,lid), PE_hyd,&
                     dP,characteristic_curves, &
                     material_auxvar, option, HYDRATE_FORMER_CH4)
      call HydratePE(hyd_auxvar%temp,(hyd_auxvar%sat(hid)+ &
                     hyd_auxvar%sat(h2id)),hyd_auxvar%xmol(cid,lid),&
                     PE_hyd2,dP2,characteristic_curves, &
                     material_auxvar, option, HYDRATE_FORMER_CO2)

      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CH4,K_H_tilde)
      call HenrysConstant(hyd_auxvar%temp,HYDRATE_FORMER_CO2,K_H_tilde2)

      if (hydrate_adjust_ghsz_solubility) then
        call HydrateGHSZSolubilityCorrection(hyd_auxvar%temp,hyd_auxvar% &
             pres(gid),dP,K_H_tilde)
        call HydrateGHSZSolubilityCorrection(hyd_auxvar%temp,hyd_auxvar% &
             pres(gid),dP2,K_H_tilde2)
      endif

      call EOSWaterSaturationPressure(hyd_auxvar%temp, &
                                          hyd_auxvar%pres(spid),ierr)

      hyd_auxvar%pres(cpid) = 0.d0
      hyd_auxvar%pres(apid) = hyd_auxvar%pres(gid) !- hyd_auxvar%pres(spid)

      hyd_auxvar%pres(lid) = hyd_auxvar%pres(gid)
      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid) - hyd_auxvar%pres(apid)
!      hyd_auxvar%pres(vpid) = hyd_auxvar%pres(gid)
      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
        if (dabs(hyd_auxvar%sat(hid) + hyd_auxvar%sat(h2id)) < 1.d-40) then
          hyd_auxvar%xmol(acid,lid) = hyd_auxvar%h1_h2_ratio * &
                                      (hyd_auxvar%pres(apid) / K_H_tilde)
          hyd_auxvar%xmol(cid,lid) = (hyd_auxvar%pres(apid) / K_H_tilde) - &
                                     hyd_auxvar%xmol(acid,lid)
        else
          hyd_auxvar%h1_h2_ratio = hyd_auxvar%sat(hid)/(hyd_auxvar%sat(hid) + &
                                    hyd_auxvar%sat(h2id))
          hyd_auxvar%xmol(acid,lid) = (hyd_auxvar%pres(apid) / K_H_tilde) * &
                                    hyd_auxvar%h1_h2_ratio
          hyd_auxvar%xmol(cid,lid) = (hyd_auxvar%pres(apid) / K_H_tilde) - &
                                   hyd_auxvar%xmol(acid,lid)
        endif
      else
         hyd_auxvar%xmol(acid,lid) = PE_hyd / K_H_tilde
         hyd_auxvar%xmol(cid,lid) = PE_hyd2  / K_H_tilde2
      endif

      hyd_auxvar%xmol(wid,lid) = 1.d0 - hyd_auxvar%xmol(acid,lid) - &
                                 hyd_auxvar%xmol(cid,lid)
      hyd_auxvar%xmol(wid,gid) = 0.d0
      hyd_auxvar%xmol(acid,gid) = 0.d0


    case default
      write(option%io_buffer,*) state
      option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
        ') not recognized in HydrateAuxVarCompute.'
      call PrintErrMsgByRank(option)

  end select

  !Compute freezing point of water
  if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
    h_sat_eff = hyd_auxvar%sat(hid)+hyd_auxvar%sat(iid)
    i_sat_eff = 2.d0 * hyd_auxvar%sat(iid)
  else
    h_sat_eff = 2.d0 * hyd_auxvar%sat(hid)
    i_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid)
  endif

  if (hydrate_with_gibbs_thomson) then
    call GibbsThomsonFreezing(1.d0-i_sat_eff,6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,material_auxvar,option)
  else
    dTf = 0.d0
  endif

  hyd_auxvar%Tf_ice = TQD + dTf
  hyd_auxvar%Pe(ONE_INTEGER) = PE_hyd
  hyd_auxvar%Pe(TWO_INTEGER) = PE_hyd2
  hyd_auxvar%K_H_tilde = K_H_tilde
  hyd_auxvar%K_H_tilde2 = K_H_tilde2
  hyd_auxvar%K_H_tilde_hyd = K_H_tilde_hyd

end subroutine HydrateComputeStateSpecificVars4

! ************************************************************************** !

subroutine HydrateEOSGasError(natural_id,ierr,hyd_auxvar,option)

  !
  !  HydrateEOSGasError: Elaborates when variables exceeds the bounds of
  !                      the equation of state.
  !
  !  Author: Heeho Park
  !  Date: 5/29/19
  !

  use Option_module
  use String_module
  use EOS_Gas_module

  implicit none

  PetscInt :: natural_id
  PetscInt :: ierr
  type(hydrate_auxvar_type) :: hyd_auxvar
  type(option_type) :: option
  
  
  call PrintMsgByCell(option,natural_id, &
                      'Error in HydrateAuxVarCompute->EOSGasHenry')
  if (ierr == EOS_GAS_TEMP_BOUND_EXCEEDED) then
    option%io_buffer = 'Temperature at cell ID ' // trim(StringWrite(natural_id)) // &
                               ' exceeds the equation of state temperature bound with ' // &
                               trim(StringWrite(hyd_auxvar%temp)) // ' [C].'
    call PrintErrMsgByRank(option)
    hydrate_high_temp_ts_cut = PETSC_TRUE
  endif

  
end subroutine HydrateEOSGasError

! ************************************************************************** !

subroutine HydrateAuxVarUpdateState(x,hyd_auxvar,global_auxvar, &
                                    material_auxvar, &
                                    characteristic_curves,natural_id, &
                                    option)

  !
  ! Decides on state changes and adds epsilons to new primary variables
  ! accordingly.
  !
  ! Author: Michael Nole
  ! Date: 01/28/18
  !

  use Option_module
  use Global_Aux_module
  use EOS_Water_module
  use Characteristic_Curves_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  class(characteristic_curves_type) :: characteristic_curves
  type(hydrate_auxvar_type) :: hyd_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscReal, parameter :: epsilon = 0.d0
  PetscReal :: x(option%nflowdof)
  PetscReal :: PE_hyd,PE_hyd2, dP, K_H, Tf_ice, dTf, h_sat_eff, i_sat_eff
  PetscReal :: K_H_tilde, K_H_tilde2, K_H_tilde_hyd
  PetscReal :: xmol_compare
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, hid, iid, acid, wid, h2id, cid
  PetscInt :: old_state,new_state
  PetscBool :: istatechng
  PetscErrorCode :: ierr
  PetscReal :: epsilons(4)
  character(len=MAXSTRINGLENGTH) :: state_change_string, append

  if (hydrate_immiscible .or. hyd_auxvar%istatechng) return

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = 3
  iid = 4
  h2id = 5

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id
  cid = option%co2_id

  hyd_auxvar%istate_store(PREV_IT) = global_auxvar%istate
  istatechng = PETSC_FALSE

  !man: right now comparing hydrate equilib pressure to gas
  !pressure (assuming low water solubility in methane).
  !Ideally would compare to partial pressure of methane.

  old_state = global_auxvar%istate

  if (option%nflowdof == 3) then
    if (global_auxvar%istate == ZERO_INTEGER) global_auxvar% &
                                istate = global_auxvar%istate
    hyd_auxvar%istate_store(PREV_IT) = global_auxvar%istate

    if (hyd_auxvar%sat(hid) > hyd_auxvar%sat(iid)) then
      h_sat_eff = hyd_auxvar%sat(hid)+hyd_auxvar%sat(iid)
      i_sat_eff = 2.d0 * hyd_auxvar%sat(iid)
    else
      h_sat_eff = 2.d0 * hyd_auxvar%sat(hid)
      i_sat_eff = hyd_auxvar%sat(hid) + hyd_auxvar%sat(iid)
    endif

    call HydratePE(hyd_auxvar%temp,h_sat_eff, global_auxvar%m_nacl(1), &
                   PE_hyd, dP, characteristic_curves, material_auxvar, &
                   option, hydrate_former)
    call HenrysConstant(hyd_auxvar%temp,hydrate_former,K_H_tilde)
    K_H_tilde_hyd = K_H_tilde
    if (hydrate_adjust_ghsz_solubility) then
          call HydrateGHSZSolubilityCorrection(hyd_auxvar%temp,hyd_auxvar% &
                                             pres(lid),dP,K_H_tilde_hyd)
    endif
    if (hydrate_with_gibbs_thomson) then
    call GibbsThomsonFreezing(1.d0-i_sat_eff,6017.1d0,ICE_DENSITY,TQD,dTf, &
                            characteristic_curves,material_auxvar,option)
    else
      dTf = 0.d0
    endif

    Tf_ice = TQD - dTf

    hyd_auxvar%Pe(ONE_INTEGER) = PE_hyd 
    hyd_auxvar%Tf_ice = Tf_ice 
    hyd_auxvar%K_H_tilde = K_H_tilde 
    hyd_auxvar%K_H_tilde_hyd = K_H_tilde_hyd


    call HydrateGetNewState3(hyd_auxvar,istatechng,global_auxvar%istate, &
                             epsilons,option)
  elseif (option%nflowdof == 4) then
    call HydrateGetNewState4(hyd_auxvar,istatechng,global_auxvar%istate, &
                             epsilons,option)
  endif

  hydrate_state_changed = istatechng

  new_state = global_auxvar%istate

  !Update primary variables

  if (istatechng) then

    select case (old_state)
      case(L_STATE)
        state_change_string = 'Liquid --> '
      case(G_STATE)
        state_change_string = 'Gas --> '
      case(H_STATE)
        state_change_string = 'Hydrate --> '
      case(I_STATE)
        state_change_string = 'Ice --> '
      case(GA_STATE)
        state_change_string = 'Gas-Aqueous --> '
      case(HG_STATE)
        state_change_string = 'Hydrate-Gas --> '
      case(HA_STATE)
        state_change_string = 'Hydrate-Aqueous --> '
      case(HI_STATE)
        state_change_string = 'Hydrate-Ice --> '
      case(GI_STATE)
        state_change_string = 'Gas-Ice --> '
      case(AI_STATE)
        state_change_string = 'Aqueous-Ice --> '
      case(HGA_STATE)
        state_change_string = 'Hydrate-Gas-Aqueous --> '
      case(HAI_STATE)
        state_change_string = 'Hydrate-Aqueous-Ice --> '
      case(HGI_STATE)
        state_change_string = 'Hydrate-Gas-Ice --> '
      case(GAI_STATE)
        state_change_string = 'Gas-Aqueous-Ice --> '
      case(QUAD_STATE)
        state_change_string = '4-Phase --> '
      case(H2A_STATE)
        state_change_string = 'H2-Aqueous --> '
      case(HHA_STATE)
        state_change_string = 'H1-H2-Aqueous --> '
      case(HHGA_STATE)
        state_change_string = 'H1-H2-Gas-Aqueous -->'

    end select
    select case (new_state)
      case(L_STATE)
        state_change_string = trim(state_change_string) // ' Liquid'
      case(G_STATE)
        state_change_string = trim(state_change_string) // ' Gas'
      case(H_STATE)
        state_change_string = trim(state_change_string) // ' Hydrate'
      case(I_STATE)
        state_change_string = trim(state_change_string) // ' Ice'
      case(GA_STATE)
        state_change_string = trim(state_change_string) // ' Gas-Aqueous'
      case(HG_STATE)
        state_change_string = trim(state_change_string) // ' Hydrate-Gas'
      case(HA_STATE)
        state_change_string = trim(state_change_string) // ' Hydrate-Aqueous'
      case(HI_STATE)
        state_change_string = trim(state_change_string) // ' Hydrate-Ice'
      case(GI_STATE)
        state_change_string = trim(state_change_string) // ' Gas-Ice'
      case(AI_STATE)
        state_change_string = trim(state_change_string) // ' Aqueous-Ice'
      case(HGA_STATE)
        state_change_string = trim(state_change_string) // &
                              ' Hydrate-Gas-Aqueous'
      case(HAI_STATE)
        state_change_string = trim(state_change_string) // &
                              ' Hydrate-Aqueous-Ice'
      case(HGI_STATE)
        state_change_string = trim(state_change_string) // ' Hydrate-Gas-Ice'
      case(GAI_STATE)
        state_change_string = trim(state_change_string) // ' Gas-Aqueous-Ice'
      case(QUAD_STATE)
        state_change_string = trim(state_change_string) // ' 4-Phase'
      case(H2A_STATE)
        state_change_string = trim(state_change_string) // ' H2-Aqueous'
      case(HHA_STATE)
        state_change_string = trim(state_change_string) // ' H1-H2-Aqueous'
      case(HHGA_STATE)
        state_change_string = trim(state_change_string) // ' H1-H2-Gas-Aqueous'
    end select

    if (option%iflag == HYDRATE_UPDATE_FOR_ACCUM) then
      write(append,'('' at Cell '',i8)') natural_id
    else if (option%iflag == HYDRATE_UPDATE_FOR_DERIVATIVE) then
      write(append, &
            '(''(due to perturbation) '',i8)') natural_id
    else
      write(append, &
             '('' at Boundary Face '', i8)') natural_id
    endif

    state_change_string = trim(state_change_string) // trim(append)

    if (hydrate_restrict_state_chng) hyd_auxvar%istatechng = PETSC_TRUE

    if (option%nflowdof == 3) then
      call HydrateUpdatePrimaryVars3(x,global_auxvar%istate,hyd_auxvar, &
                                     epsilons,option)
      call HydrateAuxVarCompute3(x,hyd_auxvar,global_auxvar,material_auxvar, &
          characteristic_curves,natural_id,option)
    elseif (option%nflowdof == 4) then
      call HydrateUpdatePrimaryVars4(x,global_auxvar%istate,hyd_auxvar, &
                                     epsilons,option)
      call HydrateAuxVarCompute4(x,hyd_auxvar,global_auxvar,material_auxvar, &
          characteristic_curves,natural_id,option)
    endif

    state_change_string = 'State Transition: ' // trim(state_change_string)
    if (hydrate_print_state_transition) then
      call PrintMsgByRank(option,state_change_string)
    endif

  endif

end subroutine HydrateAuxVarUpdateState

! ************************************************************************** !

subroutine HydrateGetNewState3(hyd_auxvar,istatechng,state,epsilons,option)
  !
  ! Finds the new phase state for a system of 3 equations
  !
  ! Author: Michael Nole
  ! Date: 08/26/20
  !

  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: hyd_auxvar
  PetscBool :: istatechng
  PetscInt :: state
  type(option_type) :: option

  PetscReal :: PE_hyd, Tf_ice
  PetscInt :: gid, lid, acid, wid, eid, hid, iid, h2id, cid, &
              apid, cpid, vpid, spid
  PetscReal :: liq_epsilon, gas_epsilon, two_phase_epsilon, &
               ha_epsilon
  PetscReal :: K_H_tilde, K_H_tilde_hyd
  PetscReal :: epsilons(4)

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  h2id = option%hydrate2_phase ! co2 hydrate

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id
  wid = option%water_id
  cid = option%co2_id
  eid = option%energy_id

  liq_epsilon = 0.d0
  gas_epsilon = 0.d0
  two_phase_epsilon = 0.d0
  ha_epsilon = 0.d0

  PE_hyd = hyd_auxvar%Pe(ONE_INTEGER)
  Tf_ice = hyd_auxvar%Tf_ice
  K_H_tilde = hyd_auxvar%K_H_tilde
  K_H_tilde_hyd = hyd_auxvar%K_H_tilde_hyd

  select case(state)
    case(L_STATE)
      if (hyd_auxvar%temp > Tf_ice) then
        if (hyd_auxvar%pres(lid) >= PE_hyd .and. &
            K_H_tilde_hyd*hyd_auxvar%xmol(acid,lid) >= hyd_auxvar% &
             pres(lid)*(1.d0-window_epsilon)) then
            istatechng = PETSC_TRUE
            state = HA_STATE
        elseif (hyd_auxvar%pres(lid) <= PE_hyd .and. &
            K_H_tilde*hyd_auxvar%xmol(acid,lid) >= hyd_auxvar% &
             pres(lid)*(1.d0-window_epsilon)) then
            istatechng = PETSC_TRUE
            state = GA_STATE
            liq_epsilon = hydrate_phase_chng_epsilon
        else
          istatechng = PETSC_FALSE
        endif
      elseif (hyd_auxvar%pres(apid) >= hyd_auxvar%pres(lid)* &
              (1.d0-window_epsilon)) then
        if (hyd_auxvar%pres(gid) < PE_hyd) then
          istatechng = PETSC_TRUE
          state = GAI_STATE
          liq_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%pres(gid) > PE_hyd) then
          istatechng = PETSC_TRUE
          state = HAI_STATE
        else
          istatechng = PETSC_TRUE
          state = QUAD_STATE
          liq_epsilon = hydrate_phase_chng_epsilon
        endif
      else
        istatechng = PETSC_TRUE
        state = AI_STATE
      endif

    case(G_STATE)
      if (hyd_auxvar%pres(vpid) >= hyd_auxvar%pres(spid)* &
         (1.d0-window_epsilon)) then
        if (hyd_auxvar%pres(apid) < PE_hyd) then
          if (hyd_auxvar%temp > Tf_ice) then
            istatechng = PETSC_TRUE
            state = GA_STATE
            gas_epsilon = hydrate_phase_chng_epsilon

          elseif (hyd_auxvar%temp == Tf_ice) then
            istatechng = PETSC_TRUE
            state = GAI_STATE
            gas_epsilon = hydrate_phase_chng_epsilon

          else
            istatechng = PETSC_TRUE
            state = GI_STATE
            gas_epsilon = hydrate_phase_chng_epsilon
          endif
        else
          istatechng = PETSC_TRUE
          state = HG_STATE
          gas_epsilon = hydrate_phase_chng_epsilon

        endif
      else
        istatechng = PETSC_FALSE
      endif

    case(H_STATE)
      if (hyd_auxvar%pres(apid) < PE_hyd) then
      !if (hyd_auxvar%pres(gid) < PE_hyd) then
        if (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          state = HGA_STATE
        elseif (hyd_auxvar%temp == Tf_ice) then
          istatechng = PETSC_TRUE
          state = QUAD_STATE
        else
          istatechng = PETSC_TRUE
          state = HGI_STATE
        endif
      else
          istatechng = PETSC_FALSE
      endif

    case(I_STATE)
      if (hyd_auxvar%temp > Tf_ice) then
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      else
        istatechng = PETSC_FALSE
      endif

   case(GA_STATE)
      if (hyd_auxvar%pres(apid) < PE_hyd) then
      !if (hyd_auxvar%pres(gid) < PE_hyd) then
        if (hyd_auxvar%temp > Tf_ice) then
          if (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) > 0.d0) then
            istatechng = PETSC_FALSE
          elseif (hyd_auxvar%sat(gid) <= 0.d0) then
            istatechng = PETSC_TRUE
            state = L_STATE
            two_phase_epsilon = hydrate_phase_chng_epsilon
          elseif (hyd_auxvar%sat(gid) >= 1.d0) then
            istatechng = PETSC_TRUE
            state = G_STATE
            two_phase_epsilon = hydrate_phase_chng_epsilon
          endif
        else
          istatechng = PETSC_TRUE
          state = GAI_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        endif
      else
        if (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          state = HGA_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        else
          istatechng = PETSC_TRUE
          state = QUAD_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        endif
      endif

    case(HG_STATE)
      if (hyd_auxvar%pres(apid) > PE_hyd) then
      !if (hyd_auxvar%pres(gid) > PE_hyd) then
        if (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(hid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = H_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = G_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        endif
      else
        if (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          state = HGA_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%temp == Tf_ice) then
          istatechng = PETSC_TRUE
          state = QUAD_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        else
          istatechng = PETSC_TRUE
          state = HGI_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        endif
      endif

    case(HA_STATE)
      !if (hyd_auxvar%pres(gid) > PE_hyd .and. hyd_auxvar%temp > Tf_ice) then
      if (hyd_auxvar%pres(apid) > PE_hyd .and. hyd_auxvar%temp > Tf_ice) then
        if (hyd_auxvar%sat(hid) >0.d0 .and. hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(hid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = H_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = L_STATE
        endif
      elseif (hyd_auxvar%temp > Tf_ice) then
        istatechng = PETSC_TRUE
        state = HGA_STATE
        ha_epsilon = hydrate_phase_chng_epsilon
      elseif (hyd_auxvar%pres(apid) > PE_hyd) then
       istatechng = PETSC_TRUE
       state = HAI_STATE
      else
       istatechng = PETSC_TRUE
       state = QUAD_STATE
      endif

    case(HI_STATE)
      if (hyd_auxvar%pres(apid) > PE_hyd) then
      !if (hyd_auxvar%pres(gid) > PE_hyd) then
        if (hyd_auxvar%temp < Tf_ice) then
          istatechng = PETSC_FALSE
        else
          istatechng = PETSC_TRUE
          state = HAI_STATE
        endif
      else
        if (hyd_auxvar%temp < Tf_ice) then
          istatechng = PETSC_TRUE
          state = HGI_STATE
        else
          istatechng = PETSC_TRUE
          state = QUAD_STATE
        endif
      endif

    case(GI_STATE)
      if (hyd_auxvar%temp < Tf_ice .and. hyd_auxvar%pres(apid) < PE_hyd) then
      !if (hyd_auxvar%temp < Tf_ice .and. hyd_auxvar%pres(gid) < PE_hyd) then
        if (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = G_STATE
        else
          istatechng = PETSC_TRUE
          state = I_STATE
        endif
      elseif (hyd_auxvar%temp < Tf_ice) then
        istatechng = PETSC_TRUE
        state = HGI_STATE
      elseif (hyd_auxvar%pres(apid) < PE_hyd) then
        istatechng = PETSC_TRUE
        state = GAI_STATE
      else
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      endif

    case(AI_STATE)
      !if (hyd_auxvar%pres(apid) >= hyd_auxvar% &
      !       pres(lid)*(1.d0-window_epsilon)) then 
      !  if (hyd_auxvar%pres(apid) < PE_hyd) then
      if (hyd_auxvar%pres(apid) > PE_hyd*(1.d0-window_epsilon)) then
        if (hyd_auxvar%pres(lid) < PE_hyd*(1.d0+window_epsilon)) then
      !  if (hyd_auxvar%pres(lid) >  PE_hyd*(1.d0-window_epsilon)) then
          istatechng = PETSC_TRUE
          state = GAI_STATE
        else
          istatechng = PETSC_TRUE
          state = HAI_STATE
        endif
      else
        if (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = L_STATE
        else
          istatechng = PETSC_TRUE
          state = I_STATE
        endif
      endif

    case(HGA_STATE)
      if (hyd_auxvar%temp > Tf_ice) then
        if (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) > 0.d0 &
            .and. hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HA_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HG_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = GA_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = L_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = G_STATE
        else
          istatechng = PETSC_TRUE
          state = H_STATE
        endif
      else
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      endif

    case(HAI_STATE)
      if (hyd_auxvar%pres(apid) >= PE_hyd*(1.d0-window_epsilon)) then
      !if (hyd_auxvar%pres(gid) > PE_hyd*(1.d0-window_epsilon)) then
        if (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(hid) > 0.d0 &
            .and. hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = AI_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(hid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HA_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HI_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = L_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = H_STATE
        elseif (hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = I_STATE
        endif
      else
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      endif

    case(HGI_STATE)
      if (hyd_auxvar%temp < Tf_ice) then
        if (hyd_auxvar%sat(iid) > 0.d0 .and. hyd_auxvar%sat(hid) > 0.d0 &
            .and. hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(hid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HG_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HI_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = GI_STATE
        elseif (hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = I_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = G_STATE
        else
          istatechng = PETSC_TRUE
          state = H_STATE
        endif
      else
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      endif

    case(GAI_STATE)
      if (hyd_auxvar%pres(apid) < PE_hyd) then
      !if (hyd_auxvar%pres(gid) < PE_hyd) then
        if (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) > 0.d0 &
            .and. hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) &
                > 0.d0) then
          istatechng = PETSC_TRUE
          state = GA_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) &
                > 0.d0) then
          istatechng = PETSC_TRUE
          state = GI_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) &
                > 0.d0) then
          istatechng = PETSC_TRUE
          state = AI_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = G_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = L_STATE
        else
          istatechng = PETSC_TRUE
          state = I_STATE
        endif
      else
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      endif

    case(QUAD_STATE)
      if (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(gid) > 0.d0 &
          .and. hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) &
          > 0.d0) then
        istatechng = PETSC_FALSE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) &
               > 0.d0 .and. hyd_auxvar%sat(lid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = HGA_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = HAI_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = GAI_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) &
              > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = HGI_STATE
      elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = AI_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = GA_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = HA_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = HI_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = GI_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(hid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = HG_STATE
      elseif (hyd_auxvar%sat(lid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = L_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = H_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = G_STATE
      else
        istatechng = PETSC_TRUE
        state = I_STATE
      endif
  end select
  epsilons(1) = liq_epsilon
  epsilons(2) = gas_epsilon
  epsilons(3) = two_phase_epsilon
  epsilons(4) = ha_epsilon

end subroutine HydrateGetNewState3
! ************************************************************************** !

subroutine HydrateGetNewState4(hyd_auxvar,istatechng,state,epsilons,option)
  !
  ! Finds the new phase state for a system of 4 equations
  !
  ! Author: Michael Nole
  ! Date: 08/26/20
  !

  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: hyd_auxvar
  PetscBool :: istatechng
  PetscInt :: state
  type(option_type) :: option

  PetscReal :: PE_hyd, PE_hyd2, Tf_ice, xmol_compare
  PetscInt :: gid, lid, acid, wid, eid, hid, iid, h2id, cid, &
              apid, cpid, vpid, spid
  PetscReal :: liq_epsilon, gas_epsilon, two_phase_epsilon, &
               ha_epsilon
  PetscReal :: K_H_tilde, K_H_tilde2, K_H_tilde_hyd
  PetscReal :: epsilons(4)

  liq_epsilon = 0.d0
  gas_epsilon = 0.d0
  two_phase_epsilon = 0.d0
  ha_epsilon = 0.d0

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  h2id = option%hydrate2_phase ! co2 hydrate

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id
  wid = option%water_id
  cid = option%co2_id
  eid = option%energy_id

  PE_hyd = hyd_auxvar%Pe(ONE_INTEGER)
  PE_hyd2 = hyd_auxvar%Pe(TWO_INTEGER)
  Tf_ice = hyd_auxvar%Tf_ice
  K_H_tilde = hyd_auxvar%K_H_tilde
  K_H_tilde2 = hyd_auxvar%K_H_tilde2
  K_H_tilde_hyd = hyd_auxvar%K_H_tilde_hyd

  if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
    xmol_compare = hyd_auxvar%xmol(acid,lid) + hyd_auxvar%xmol(cid,lid)
  else
    xmol_compare = hyd_auxvar%xmol(acid,lid)
  endif

  select case(state)
    case(L_STATE)
      if (hyd_auxvar%temp > Tf_ice) then
        if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
          if (hyd_auxvar%pres(lid) >= PE_hyd .and. &
              K_H_tilde_hyd*xmol_compare >= hyd_auxvar% &
              pres(lid)*(1.d0-window_epsilon)) then

            if (hyd_auxvar%xmol(acid,lid) > 0.d0 .and. hyd_auxvar% &
                xmol(cid,lid) > 0.d0 .and. hyd_auxvar%h1_h2_ratio < 9.99d-1 &
                .and. hyd_auxvar%h1_h2_ratio > 1.d-3) then
              istatechng = PETSC_TRUE
              state = HHA_STATE
            elseif (hyd_auxvar%xmol(acid,lid) > 0.d0 .and. &
                    hyd_auxvar%h1_h2_ratio >= 1.d3) then
              istatechng = PETSC_TRUE
              state = HA_STATE
            elseif (hyd_auxvar%xmol(cid,lid) > 0.d0 .and. &
                    hyd_auxvar%h1_h2_ratio < 1.d-3) then
              istatechng = PETSC_TRUE
              state = H2A_STATE
            endif

          elseif (hyd_auxvar%pres(lid) < PE_hyd .and. &
                  K_H_tilde*xmol_compare >= hyd_auxvar% &
                  pres(lid)*(1.d0-window_epsilon)) then

            istatechng = PETSC_TRUE
            state = GA_STATE
            liq_epsilon = hydrate_phase_chng_epsilon

          else
            istatechng = PETSC_FALSE
          endif
        ! elseif (hydrate_salinity) then
        !   if (hyd_auxvar%pres(lid) >= PE_hyd .and. &
        !       K_H_tilde_hyd * hyd_auxvar%xmol(acid,lid) > &
        !       hyd_auxvar%pres(lid)) then
        !     istatechng = PETSC_TRUE
        !     state = HA_STATE
        !   endif
        else !different hydrate species
          if (hyd_auxvar%pres(lid) >= PE_hyd .and. &
               K_H_tilde_hyd*xmol_compare >= PE_hyd*(1.d0-window_epsilon)) then
            istatechng = PETSC_TRUE
            if (hyd_auxvar%xmol(acid,lid) > 0.d0 .and. hyd_auxvar% &
                 xmol(cid,lid) > 0.d0) then
              state = HA_STATE
            elseif (hyd_auxvar%xmol(acid,lid) > 0.d0) then
              state = HA_STATE
            else
              state = H2A_STATE
            endif
          elseif (hyd_auxvar%pres(lid) <= PE_hyd .and. &
               K_H_tilde*xmol_compare >= hyd_auxvar% &
               pres(lid)*(1.d0-window_epsilon)) then

            istatechng = PETSC_TRUE
            state = GA_STATE
            liq_epsilon = hydrate_phase_chng_epsilon
          else
             istatechng = PETSC_FALSE
          endif
        endif
      elseif (hyd_auxvar%pres(apid) >= hyd_auxvar%pres(lid)* &
              (1.d0-window_epsilon)) then
        if (hyd_auxvar%pres(gid) < PE_hyd) then
          istatechng = PETSC_TRUE
          state = GAI_STATE
          liq_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%pres(gid) > PE_hyd) then
          istatechng = PETSC_TRUE
          state = HAI_STATE
        else
          istatechng = PETSC_TRUE
          state = QUAD_STATE
          liq_epsilon = hydrate_phase_chng_epsilon
        endif
      else
        istatechng = PETSC_TRUE
        state = AI_STATE
      endif

    case(G_STATE)
      if (hyd_auxvar%pres(vpid) >= hyd_auxvar%pres(spid)* &
         (1.d0-window_epsilon)) then
        if (hyd_auxvar%pres(apid) < PE_hyd) then
          if (hyd_auxvar%temp > Tf_ice) then
            istatechng = PETSC_TRUE
            state = GA_STATE
            gas_epsilon = hydrate_phase_chng_epsilon

          elseif (hyd_auxvar%temp == Tf_ice) then
            istatechng = PETSC_TRUE
            state = GAI_STATE
            gas_epsilon = hydrate_phase_chng_epsilon

          else
            istatechng = PETSC_TRUE
            state = GI_STATE
            gas_epsilon = hydrate_phase_chng_epsilon
          endif
        else
          istatechng = PETSC_TRUE
          state = HG_STATE
          gas_epsilon = hydrate_phase_chng_epsilon

        endif
      else
        istatechng = PETSC_FALSE
      endif

    case(H_STATE)
      if (hyd_auxvar%pres(apid) < PE_hyd) then
        if (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          state = HGA_STATE

       elseif (hyd_auxvar%temp == Tf_ice) then
          istatechng = PETSC_TRUE
          state = QUAD_STATE
        else
          istatechng = PETSC_TRUE
          state = HGI_STATE
        endif
      else
          istatechng = PETSC_FALSE
      endif

    case(I_STATE)
      if (hyd_auxvar%temp > Tf_ice) then
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      else
        istatechng = PETSC_FALSE
      endif

   case(GA_STATE)
      if (hyd_auxvar%pres(apid) < PE_hyd) then
        if (hyd_auxvar%temp > Tf_ice) then
          if (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) > 0.d0) then
            istatechng = PETSC_FALSE
          elseif (hyd_auxvar%sat(gid) <= 0.d0) then
            istatechng = PETSC_TRUE
            state = L_STATE
            two_phase_epsilon = hydrate_phase_chng_epsilon
          elseif (hyd_auxvar%sat(gid) >= 1.d0) then
            istatechng = PETSC_TRUE
            state = G_STATE
            two_phase_epsilon = hydrate_phase_chng_epsilon
          endif
        else
          istatechng = PETSC_TRUE
          state = GAI_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        endif
      else
        if (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          state = HGA_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        else
          istatechng = PETSC_TRUE
          state = QUAD_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        endif
      endif

    case(HG_STATE)
      if (hyd_auxvar%pres(apid) > PE_hyd) then
        if (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(hid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = H_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = G_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        endif
      else
        if (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          state = HGA_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%temp == Tf_ice) then
          istatechng = PETSC_TRUE
          state = QUAD_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        else
          istatechng = PETSC_TRUE
          state = HGI_STATE
          two_phase_epsilon = hydrate_phase_chng_epsilon
        endif
      endif
    case(HA_STATE)
      if (hyd_auxvar%pres(gid) > PE_hyd .and. hyd_auxvar%temp > Tf_ice) then
        if (hyd_auxvar%sat(hid) >0.d0 .and. hyd_auxvar%sat(lid) > 0.d0) then
          if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then! .and. &
            if (hyd_auxvar%xmol(cid,lid) > 0.d0 .and. &
                hyd_auxvar%h1_h2_ratio < 9.99d-1) then
              istatechng = PETSC_TRUE
              state = HHA_STATE
            else
              istatechng = PETSC_FALSE
            endif
          elseif (hydrate_salinity) then
            istatechng = PETSC_FALSE
          else !different species hydrate
            if (hyd_auxvar%xmol(cid,lid)*K_H_tilde2 >PE_hyd2) then
              istatechng = PETSC_TRUE
              state = HHA_STATE
            else
              istatechng = PETSC_FALSE
            endif
          endif
        elseif (hyd_auxvar%sat(hid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = H_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = L_STATE
        endif
      ! elseif (hyd_auxvar%temp > Tf_ice) then
      !   istatechng = PETSC_TRUE
      !   state = HGA_STATE
      !   ha_epsilon = hydrate_phase_chng_epsilon
      ! elseif (hyd_auxvar%pres(apid) > PE_hyd) then
      !  istatechng = PETSC_TRUE
      !  state = HAI_STATE
      ! else
      !  istatechng = PETSC_TRUE
      !  state = QUAD_STATE
      endif

    case(H2A_STATE)
      !if (hyd_auxvar%pres(gid) > PE_hyd .and. hyd_auxvar%temp > Tf_ice) then
      !MAN: this needs to check the partial pressure of CO2 against CO2 hydrate
      !     equilibrium pressure when mass balance 2 is CO2
      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
        if (hyd_auxvar%pres(apid) > PE_hyd .and. hyd_auxvar%temp > Tf_ice) then
           if (hyd_auxvar%sat(h2id) >0.d0 .and. hyd_auxvar%sat(lid) > 0.d0) then
             if (hyd_auxvar%xmol(acid,lid) > 0.d0 .and. &
                 hyd_auxvar%h1_h2_ratio > 1.d-3) then
               istatechng = PETSC_TRUE
               state = HHA_STATE
             else
               istatechng = PETSC_FALSE
             endif
           elseif (hyd_auxvar%sat(h2id) > 0.d0) then
             istatechng = PETSC_TRUE
             state = H2_STATE
           elseif (hyd_auxvar%sat(lid) > 0.d0) then
             istatechng = PETSC_TRUE
             state = L_STATE
           endif
        elseif (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          state = HGA_STATE
          ha_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%pres(apid) > PE_hyd) then
          istatechng = PETSC_TRUE
          state = HAI_STATE
        else
          istatechng = PETSC_TRUE
          state = QUAD_STATE
        endif
      endif

    case(HI_STATE)
      if (hyd_auxvar%pres(apid) > PE_hyd) then
      !if (hyd_auxvar%pres(gid) > PE_hyd) then
        if (hyd_auxvar%temp < Tf_ice) then
          istatechng = PETSC_FALSE
        else
          istatechng = PETSC_TRUE
          state = HAI_STATE
        endif
      else
        if (hyd_auxvar%temp < Tf_ice) then
          istatechng = PETSC_TRUE
          state = HGI_STATE
        else
          istatechng = PETSC_TRUE
          state = QUAD_STATE
        endif
      endif
    case(GI_STATE)
      if (hyd_auxvar%temp < Tf_ice .and. hyd_auxvar%pres(apid) < PE_hyd) then
      !if (hyd_auxvar%temp < Tf_ice .and. hyd_auxvar%pres(gid) < PE_hyd) then
        if (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = G_STATE
        else
          istatechng = PETSC_TRUE
          state = I_STATE
        endif
      elseif (hyd_auxvar%temp < Tf_ice) then
        istatechng = PETSC_TRUE
        state = HGI_STATE
      elseif (hyd_auxvar%pres(apid) < PE_hyd) then
        istatechng = PETSC_TRUE
        state = GAI_STATE
      else
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      endif
    case(AI_STATE)
      !if (hyd_auxvar%pres(apid) >= hyd_auxvar% &
      !       pres(lid)*(1.d0-window_epsilon)) then
      !  if (hyd_auxvar%pres(apid) < PE_hyd) then
      if (hyd_auxvar%pres(apid) > PE_hyd*(1.d0-window_epsilon)) then
        if (hyd_auxvar%pres(lid) < PE_hyd*(1.d0+window_epsilon)) then
      !  if (hyd_auxvar%pres(lid) >  PE_hyd*(1.d0-window_epsilon)) then
          istatechng = PETSC_TRUE
          state = GAI_STATE
        else
          istatechng = PETSC_TRUE
          state = HAI_STATE
        endif
      else
        if (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = L_STATE
        else
          istatechng = PETSC_TRUE
          state = I_STATE
        endif
      endif
    case(HGA_STATE)
      if (hyd_auxvar%temp > Tf_ice) then
        if (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) > 0.d0 &
            .and. hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HA_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HG_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = GA_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = L_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = G_STATE
        else
          istatechng = PETSC_TRUE
          state = H_STATE
        endif
      else
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      endif
    case(HAI_STATE)
      if (hyd_auxvar%pres(apid) >= PE_hyd*(1.d0-window_epsilon)) then
      !if (hyd_auxvar%pres(gid) > PE_hyd*(1.d0-window_epsilon)) then
        if (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(hid) > 0.d0 &
            .and. hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = AI_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(hid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HA_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HI_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = L_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = H_STATE
        elseif (hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = I_STATE
        endif
      else
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      endif
    case(HGI_STATE)
      if (hyd_auxvar%temp < Tf_ice) then
        if (hyd_auxvar%sat(iid) > 0.d0 .and. hyd_auxvar%sat(hid) > 0.d0 &
            .and. hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(hid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HG_STATE
        elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = HI_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) > &
                0.d0) then
          istatechng = PETSC_TRUE
          state = GI_STATE
        elseif (hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = I_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = G_STATE
        else
          istatechng = PETSC_TRUE
          state = H_STATE
        endif
      else
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      endif

    case(GAI_STATE)
      if (hyd_auxvar%pres(apid) < PE_hyd) then
      !if (hyd_auxvar%pres(gid) < PE_hyd) then
        if (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) > 0.d0 &
            .and. hyd_auxvar%sat(iid) > 0.d0) then
          istatechng = PETSC_FALSE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) &
                > 0.d0) then
          istatechng = PETSC_TRUE
          state = GA_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) &
                > 0.d0) then
          istatechng = PETSC_TRUE
          state = GI_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) &
                > 0.d0) then
          istatechng = PETSC_TRUE
          state = AI_STATE
        elseif (hyd_auxvar%sat(gid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = G_STATE
        elseif (hyd_auxvar%sat(lid) > 0.d0) then
          istatechng = PETSC_TRUE
          state = L_STATE
        else
          istatechng = PETSC_TRUE
          state = I_STATE
        endif
      else
        istatechng = PETSC_TRUE
        state = QUAD_STATE
      endif
    case(QUAD_STATE)
      if (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(gid) > 0.d0 &
          .and. hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) &
          > 0.d0) then
        istatechng = PETSC_FALSE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) &
               > 0.d0 .and. hyd_auxvar%sat(lid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = HGA_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = HAI_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = GAI_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(gid) &
              > 0.d0 .and. hyd_auxvar%sat(iid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = HGI_STATE
      elseif (hyd_auxvar%sat(lid) > 0.d0 .and. hyd_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = AI_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = GA_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = HA_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = HI_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0 .and. hyd_auxvar%sat(iid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = GI_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(hid) &
              > 0.d0) then
        istatechng = PETSC_TRUE
        state = HG_STATE
      elseif (hyd_auxvar%sat(lid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = L_STATE
      elseif (hyd_auxvar%sat(hid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = H_STATE
      elseif (hyd_auxvar%sat(gid) > 0.d0) then
        istatechng = PETSC_TRUE
        state = G_STATE
      else
        istatechng = PETSC_TRUE
        state = I_STATE
      endif
    case(HHA_STATE)
      if (hydrate_mass_balance_1 == hydrate_mass_balance_2) then
        if (hyd_auxvar%pres(apid) > PE_hyd .and. hyd_auxvar%temp > Tf_ice) then
          if ((hyd_auxvar%sat(hid) >0.d0 .or. hyd_auxvar%sat(h2id) > 0.d0) &
              .and. hyd_auxvar%sat(lid) > 0.d0) then
            istatechng = PETSC_FALSE
          !elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(h2id) &
          !         > 0.d0) then
          !  istatechng = PETSC_TRUE
          !  !MAN: need to implement a saturated hydrate1-hydrate2 state
          !  !state = HH_STATE
          ! MAN: I don't think you can go from HHA to HA or H2A when 
          !      both mass balances are the same species.
          !elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) > &
          !        0.d0) then
          !  istatechng = PETSC_TRUE
          !  state = HA_STATE
          !elseif (hyd_auxvar%sat(h2id) > 0.d0 .and. hyd_auxvar%sat(lid) > &
          !        0.d0) then
          !  istatechng = PETSC_TRUE
          !  state = H2A_STATE
          elseif (hyd_auxvar%sat(lid) > 0.d0) then
            istatechng = PETSC_TRUE
            state = L_STATE
          endif
        elseif (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          state = HHGA_STATE
          ha_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%pres(apid) > PE_hyd) then
         istatechng = PETSC_TRUE
         state = HAI_STATE
        else
         istatechng = PETSC_TRUE
         state = QUAD_STATE
        endif
      else
        if (hyd_auxvar%pres(apid) > PE_hyd .and. hyd_auxvar%temp > Tf_ice) then
          if ((hyd_auxvar%sat(hid) >0.d0 .or. hyd_auxvar%sat(h2id) > 0.d0) &
              .and. hyd_auxvar%sat(lid) > 0.d0) then
            istatechng = PETSC_FALSE
          !elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(h2id) &
          !         > 0.d0) then
          !  istatechng = PETSC_TRUE
          !  !MAN: need to implement a saturated hydrate1-hydrate2 state
          !  !state = HH_STATE
          ! MAN: I don't think you can go from HHA to HA or H2A when 
          !      both mass balances are the same species.
          !elseif (hyd_auxvar%sat(hid) > 0.d0 .and. hyd_auxvar%sat(lid) > &
          !        0.d0) then
          !  istatechng = PETSC_TRUE
          !  state = HA_STATE
          !elseif (hyd_auxvar%sat(h2id) > 0.d0 .and. hyd_auxvar%sat(lid) > &
          !        0.d0) then
          !  istatechng = PETSC_TRUE
          !  state = H2A_STATE
          elseif (hyd_auxvar%sat(lid) > 0.d0) then
            istatechng = PETSC_TRUE
            state = L_STATE
          endif
        elseif (hyd_auxvar%temp > Tf_ice) then
          istatechng = PETSC_TRUE
          state = HGA_STATE
          ha_epsilon = hydrate_phase_chng_epsilon
        elseif (hyd_auxvar%pres(apid) > PE_hyd) then
         istatechng = PETSC_TRUE
         state = HAI_STATE
        else
         istatechng = PETSC_TRUE
         state = QUAD_STATE
        endif
      endif
    case(H2G_STATE)
      !MAN: needs filling in
    case(H2_STATE)
      !MAN: needs filling in
  end select

  epsilons(1) = liq_epsilon
  epsilons(2) = gas_epsilon
  epsilons(3) = two_phase_epsilon
  epsilons(4) = ha_epsilon

end subroutine HydrateGetNewState4

! ************************************************************************** !

subroutine HydrateUpdatePrimaryVars3(x,state,hyd_auxvar,epsilons,option)
  !
  ! Links the solution vector with its primary variables, for 3 equations.
  !
  ! Author: Michael Nole
  ! Date: 08/26/20
  !

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  PetscInt :: state
  type(hydrate_auxvar_type) :: hyd_auxvar
  PetscReal :: epsilons(4)

  PetscInt :: gid, lid, acid, wid, eid, hid, iid, h2id, cid, &
              apid, cpid, vpid, spid
  PetscReal :: liq_epsilon, gas_epsilon, two_phase_epsilon, ha_epsilon

  liq_epsilon = epsilons(1)
  gas_epsilon = epsilons(2)
  two_phase_epsilon = epsilons(3)
  ha_epsilon = epsilons(4)

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  h2id = option%hydrate2_phase ! co2 hydrate

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id
  wid = option%water_id
  cid = option%co2_id
  eid = option%energy_id

    select case(state)

      case(L_STATE)
!     ********* Aqueous State (A) ********************************
!     Primary variables: Pg, Xma, T
!
        x(HYDRATE_LIQUID_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_L_STATE_X_MOLE_DOF) = max(0.d0,hyd_auxvar%xmol(acid,lid))
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(G_STATE)
!     ********* Gas State (G) ********************************
!     Primary variables: Pg, Pa, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_G_STATE_AIR_PRESSURE_DOF) = hyd_auxvar%pres(apid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(H_STATE)
!     ********* Hydrate State (H) ********************************
!     Primary variables: Pg, Xmh, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = MOL_RATIO_METH
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(I_STATE)
!     ********* Ice State (I) ********************************
!     Primary variables: Pg, Xmi, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = 0.d0
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(GA_STATE)
!     ********* Gas & Aqueous State (GA) ********************************
!     Primary variables: Pg, Sg, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(HG_STATE)
!     ********* Hydrate & Gas State (HG) ********************************
!     Primary variables: Pg, Sg, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(HA_STATE)
!     ********* Hydrate & Aqueous State (HA) ********************************
!     Primary variables: Pg, Sh, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(HI_STATE)
!     ********* Hydrate & Ice State (HI) ********************************
!     Primary variables: Pg, Sh, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp
      case(GI_STATE)
!     ********* Gas & Ice State (GI) ********************************
!     Primary variables: Pg, Si, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(iid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(AI_STATE)
!     ********* Aqueous & Ice State (AI) ********************************
!     Primary variables: Pg, Xma, Sl
!
        x(HYDRATE_LIQUID_PRESSURE_DOF) = hyd_auxvar%pres(lid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%xmol(acid,lid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(lid)

      case(HGA_STATE)
!     ********* Hydrate, Gas, & Aqueous State (HGA) **************************
!     Primary variables: Sl, Sh, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%sat(lid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(HAI_STATE)
!     ********* Hydrate, Aqueous, & Ice State (HAI) **************************
!     Primary variables: Pg, Sl, Si
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(lid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(iid)

      case(HGI_STATE)
!     ********* Hydrate, Gas, & Ice State (HGI) ******************************
!     Primary variables: Sh, Si, T
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%sat(iid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

      case(GAI_STATE)
!     ********* Gas, Aqueous, & Ice State (GAI) ******************************
!     Primary variables: Pg, Sl, Si
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(lid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(iid)

      case(QUAD_STATE)
!     ********* Quadruple Point (HGAI) ********************************
!     Primary variables: Sg, Sl, Si
!
        x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%sat(lid)
        x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
        x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(iid)

      case default
        write(option%io_buffer,*) state
        option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
          ') not recognized in HydrateAuxVarUpdateState.'
        call PrintErrMsgByRank(option)

    end select


end subroutine HydrateUpdatePrimaryVars3

! ************************************************************************** !

subroutine HydrateUpdatePrimaryVars4(x,state,hyd_auxvar,epsilons,option)
  !
  ! Links the solution vector with its primary variables, for 4 equations.
  !
  ! Author: Michael Nole
  ! Date: 08/26/20
  !

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal :: x(option%nflowdof)
  PetscInt :: state
  PetscReal :: epsilons(4)
  type(hydrate_auxvar_type) :: hyd_auxvar

  PetscInt :: gid, lid, acid, wid, eid, hid, iid, h2id, cid, &
              apid, cpid, vpid, spid
  PetscReal :: liq_epsilon, gas_epsilon, two_phase_epsilon, ha_epsilon

  liq_epsilon = epsilons(1)
  gas_epsilon = epsilons(2)
  two_phase_epsilon = epsilons(3)
  ha_epsilon = epsilons(4)

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  h2id = option%hydrate2_phase ! co2 hydrate

  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id
  wid = option%water_id
  cid = option%co2_id
  eid = option%energy_id

  select case(state)

    case(L_STATE)
!   ********* Aqueous State (A) ********************************
!   Primary variables: Pg, Xma, T

      x(HYDRATE_LIQUID_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_L_STATE_X_MOLE_DOF) = max(0.d0,hyd_auxvar%xmol(acid,lid))
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = max(0.d0,hyd_auxvar%xmol(cid,lid))

    case(G_STATE)
!   ********* Gas State (G) ********************************
!   Primary variables: Pg, Pa, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_G_STATE_AIR_PRESSURE_DOF) = hyd_auxvar%pres(apid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp
      x(HYDRATE_G_STATE_CO2_PRESSURE_DOF) = hyd_auxvar%pres(cid)

    case(H_STATE)
!   ********* Hydrate State (H) ********************************
!   Primary variables: Pg, Xmh, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = MOL_RATIO_FORMER
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = hyd_auxvar%xmol(cid,lid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(H2_STATE)
!   ********* CO2 Hydrate State (H2) ********************************
!   Primary variables: Pg, Xmh, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = MOL_RATIO_FORMER
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = hyd_auxvar%xmol(cid,lid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp


    case(I_STATE)
!   ********* Ice State (I) ********************************
!   Primary variables: Pg, Xmi, Sh2, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = 0.d0
      x(HYDRATE_HYD2_SATURATION_DOF) = hyd_auxvar%sat(h2id)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(HH_STATE)
!   ********* CH4 and CO2 Hydrate State (HH) **************************
!   Primary variables: Pg, Sg, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = hyd_auxvar%xmol(cid,lid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(GA_STATE)
!   ********* Gas & Aqueous State (GA) ********************************
!   Primary variables: Pg, Sg, T, Xma2

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = hyd_auxvar%xmol(cid,lid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(HG_STATE)
!   ********* Hydrate & Gas State (HG) ********************************
!   Primary variables: Pg, Sg, Pco2, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
      x(HYDRATE_G_STATE_CO2_PRESSURE_DOF) = hyd_auxvar%pres(cid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(H2G_STATE)
!   ********* CO2 Hydrate & Gas State (H2G) ********************************
!   Primary variables: Pg, Sg, Pco2, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
      x(HYDRATE_G_STATE_CO2_PRESSURE_DOF) = hyd_auxvar%pres(cid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(HA_STATE)
!   ********* Hydrate & Aqueous State (HA) ********************************
!   Primary variables: Pg, Sh, X2a, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = max(0.d0,hyd_auxvar%xmol(cid,lid))
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(HI_STATE)
!   ********* Hydrate & Ice State (HI) ********************************
!   Primary variables: Pg, Sh, X2a, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = 0.d0
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(GI_STATE)
!   ********* Gas & Ice State (GI) ********************************
!   Primary variables: Pg, Si, Pco2, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(iid)
      x(HYDRATE_G_STATE_CO2_PRESSURE_DOF) = hyd_auxvar%pres(cid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(AI_STATE)
!   ********* Aqueous & Ice State (AI) ********************************
!   Primary variables: Pg, Xma, Sl, X2a

      x(HYDRATE_LIQUID_PRESSURE_DOF) = hyd_auxvar%pres(lid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%xmol(acid,lid)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = hyd_auxvar%sat(lid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(lid)

    case(HGA_STATE)
!   ********* Hydrate, Gas, & Aqueous State (HGA) **************************
!   Primary variables: Sl, Sh, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%sat(lid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = max(0.d0,hyd_auxvar%xmol(cid,lid))
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(HAI_STATE)
!   ********* Hydrate, Aqueous, & Ice State (HAI) **************************
!   Primary variables: Pg, Sl, X2a, Si

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(lid)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = hyd_auxvar%xmol(cid,lid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(iid)

    case(HGI_STATE)
!   ********* Hydrate, Gas, & Ice State (HGI) ******************************
!   Primary variables: Sh, Si, Pco2, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%sat(iid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(hid)
      x(HYDRATE_G_STATE_CO2_PRESSURE_DOF) = hyd_auxvar%pres(cid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(GAI_STATE)
!   ********* Gas, Aqueous, & Ice State (GAI) ******************************
!   Primary variables: Pg, Sl, Si

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(lid)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = hyd_auxvar%xmol(cid,lid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(iid)
    case(QUAD_STATE)
!   ********* Quadruple Point (HGAI) ********************************
!   Primary variables: Sg, Sl, Si

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%sat(lid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%sat(gid)
      x(HYDRATE_L_STATE_X2_MOLE_DOF) = hyd_auxvar%xmol(cid,lid)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%sat(iid)

    case(H2A_STATE)
!   ********* Hydrate2 & Aqueous State (H2A) ********************************
!   Primary variables: Pg, Sh, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_L_STATE_X_MOLE_DOF) = max(0.d0,hyd_auxvar%xmol(acid,lid))
      x(HYDRATE_HYD2_SATURATION_DOF) = hyd_auxvar%sat(h2id)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(HHA_STATE)
!   ********* Hydrate1, Hydrate2, & Aqueous State (HHA) *********************
!   Primary variables: Pg, Sh1, Sh2, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%h1_h2_ratio * &
                            hydrate_phase_chng_epsilon !hyd_auxvar%sat(hid)
      x(HYDRATE_HYD2_SATURATION_DOF) = (1.d0 - hyd_auxvar%h1_h2_ratio)* &
                            hydrate_phase_chng_epsilon !hyd_auxvar%sat(h2id)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case(HHGA_STATE)
!   ********* Hydrate1, Hydrate2, Gas, & Aqueous State (HHA) *********************
!   Primary variables: Sl, Sh1, Sh2, T

      x(HYDRATE_GAS_PRESSURE_DOF) = hyd_auxvar%sat(lid)
      x(HYDRATE_GAS_SATURATION_DOF) = hyd_auxvar%h1_h2_ratio * &
                            hydrate_phase_chng_epsilon !hyd_auxvar%sat(hid)
      x(HYDRATE_HYD2_SATURATION_DOF) = (1.d0 - hyd_auxvar%h1_h2_ratio)* &
                            hydrate_phase_chng_epsilon !hyd_auxvar%sat(h2id)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar%temp

    case default
      write(option%io_buffer,*) state
      option%io_buffer = 'State (' // trim(adjustl(option%io_buffer)) // &
        ') not recognized in HydrateAuxVarUpdateState.'
      call PrintErrMsgByRank(option)

  end select

end subroutine HydrateUpdatePrimaryVars4

! ************************************************************************** !

subroutine HydrateAuxVarPerturb(hyd_auxvar,global_auxvar, &
                                material_auxvar, &
                                characteristic_curves,natural_id, &
                                option)
  !
  ! Calculates auxiliary variables for perturbed system
  !
  ! Author: Michael Nole
  ! Date: 01/30/19
  !

  use Option_module
  use Characteristic_Curves_module
  use Global_Aux_module
  use Material_Aux_class

  implicit none

  type(option_type) :: option
  PetscInt :: natural_id
  type(hydrate_auxvar_type) :: hyd_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(characteristic_curves_type) :: characteristic_curves

  PetscReal :: x(option%nflowdof), x_pert_plus(option%nflowdof), &
               pert(option%nflowdof), x_pert_minus(option%nflowdof)

  PetscReal :: tempreal
  PetscInt :: lid, gid, hid
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
!  PetscReal, parameter :: perturbation_tolerance = 1.d-11
  PetscReal, parameter :: min_mole_fraction_pert = 1.d-12
  PetscReal, parameter :: min_perturbation = 1.d-10
  
  PetscReal, parameter :: min_pres_pert = 1.d-3
  PetscReal, parameter :: min_temp_pert = 8.66d-9
  PetscReal, parameter :: min_xmol_pert = 1.d-14
  PetscReal, parameter :: min_sat_pert = 3.16d-11
  
  PetscInt :: idof

  lid = 1
  gid = 2
  hid = 3

  select case(global_auxvar%istate)
    case(L_STATE)
      x(HYDRATE_LIQUID_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%liquid_phase)
      x(HYDRATE_L_STATE_X_MOLE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%xmol(option%air_id,option%liquid_phase)
      x(HYDRATE_ENERGY_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_LIQUID_PRESSURE_DOF) = &
        perturbation_tolerance*x(HYDRATE_LIQUID_PRESSURE_DOF) + &
        min_perturbation
      if (x(HYDRATE_L_STATE_X_MOLE_DOF) > &
          1.d3 * perturbation_tolerance) then
        pert(HYDRATE_L_STATE_X_MOLE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_L_STATE_X_MOLE_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = -1.d0 * &
        (perturbation_tolerance*x(HYDRATE_ENERGY_DOF) + min_perturbation)
    case(G_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_G_STATE_AIR_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%air_pressure_id)
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar(ZERO_INTEGER)%temp
      ! gas pressure [p(g)] must always be perturbed down as p(v) = p(g) - p(a)
      ! and p(v) >= Psat (i.e. an increase in p(v)) results in two phase.

      pert(HYDRATE_GAS_PRESSURE_DOF) = -1.d0 * &
        (perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF) + min_perturbation)
      ! perturb air pressure towards gas pressure unless the perturbed
      ! air pressure exceeds the gas pressure
      tempreal = perturbation_tolerance* &
                 x(HYDRATE_G_STATE_AIR_PRESSURE_DOF) + min_perturbation
      if (x(HYDRATE_GAS_PRESSURE_DOF) - &
          x(HYDRATE_G_STATE_AIR_PRESSURE_DOF) > tempreal) then
        pert(HYDRATE_G_STATE_AIR_PRESSURE_DOF) = tempreal
      else
        pert(HYDRATE_G_STATE_AIR_PRESSURE_DOF) = -1.d0 * tempreal
      endif
      pert(HYDRATE_ENERGY_DOF) = &
        perturbation_tolerance*x(HYDRATE_ENERGY_DOF) + min_perturbation

    case(H_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = MOL_RATIO_METH
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = -1.d0 * &
         (perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF) + min_perturbation)
      pert(HYDRATE_GAS_SATURATION_DOF) = x(HYDRATE_GAS_SATURATION_DOF)
      pert(HYDRATE_ENERGY_DOF) = &
         perturbation_tolerance*x(HYDRATE_ENERGY_DOF) + min_perturbation

    case(I_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = 0.d0
      x(HYDRATE_ENERGY_DOF) = hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
        perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)
      pert(HYDRATE_GAS_SATURATION_DOF) = 0.d0
      pert(HYDRATE_ENERGY_DOF) = &
          perturbation_tolerance*x(HYDRATE_ENERGY_DOF)

    case(GA_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
       hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%sat(option%gas_phase)
      x(HYDRATE_ENERGY_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
        perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)+min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
        perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(HG_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(option%gas_phase)
      x(HYDRATE_ENERGY_DOF) = &
           hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)+min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
           perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(HA_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(hid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
        perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)+min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
         perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(HI_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(hid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)+min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
           perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(GI_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(option%gas_phase)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(iid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF)+min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
           perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(AI_STATE)
      x(HYDRATE_LIQUID_PRESSURE_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%pres(option%liquid_phase)
      x(HYDRATE_L_STATE_X_MOLE_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%xmol(option%air_id,option%liquid_phase)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(lid)
      pert(HYDRATE_LIQUID_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_LIQUID_PRESSURE_DOF) + &
         min_perturbation

      if (x(HYDRATE_L_STATE_X_MOLE_DOF) > &
           1.d3 * perturbation_tolerance) then
        pert(HYDRATE_L_STATE_X_MOLE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_L_STATE_X_MOLE_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_ENERGY_DOF) > 0.5d0) then
        pert(HYDRATE_ENERGY_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_ENERGY_DOF) = perturbation_tolerance
      endif

    case(HGA_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%sat(lid)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(hid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      pert(HYDRATE_ENERGY_DOF) = &
           perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_GAS_PRESSURE_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_PRESSURE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_PRESSURE_DOF) = perturbation_tolerance
      endif

    case(HAI_STATE)

      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(lid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(iid)

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF) + &
         min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_ENERGY_DOF) > 0.5d0) then
        pert(HYDRATE_ENERGY_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_ENERGY_DOF) = perturbation_tolerance
      endif

    case(HGI_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%sat(hid)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(iid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%temp

      if (x(HYDRATE_GAS_PRESSURE_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_PRESSURE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_PRESSURE_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      pert(HYDRATE_ENERGY_DOF) = &
           perturbation_tolerance*x(HYDRATE_ENERGY_DOF)+min_perturbation

    case(GAI_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%pres(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(lid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(iid)

      pert(HYDRATE_GAS_PRESSURE_DOF) = &
         perturbation_tolerance*x(HYDRATE_GAS_PRESSURE_DOF) + &
         min_perturbation
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_ENERGY_DOF) > 0.5d0) then
        pert(HYDRATE_ENERGY_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_ENERGY_DOF) = perturbation_tolerance
      endif

    case(QUAD_STATE)
      x(HYDRATE_GAS_PRESSURE_DOF) = &
        hyd_auxvar(ZERO_INTEGER)%sat(gid)
      x(HYDRATE_GAS_SATURATION_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(lid)
      x(HYDRATE_ENERGY_DOF) = &
         hyd_auxvar(ZERO_INTEGER)%sat(iid)

      if (x(HYDRATE_GAS_PRESSURE_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_PRESSURE_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_PRESSURE_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_GAS_SATURATION_DOF) > 0.5d0) then
        pert(HYDRATE_GAS_SATURATION_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_GAS_SATURATION_DOF) = perturbation_tolerance
      endif
      if (x(HYDRATE_ENERGY_DOF) > 0.5d0) then
        pert(HYDRATE_ENERGY_DOF) = -1.d0 * perturbation_tolerance
      else
        pert(HYDRATE_ENERGY_DOF) = perturbation_tolerance
      endif

  end select
  ! HYDRATE_UPDATE_FOR_DERIVATIVE indicates call from perturbation  

  option%iflag = HYDRATE_UPDATE_FOR_DERIVATIVE
 
  do idof = 1, option%nflowdof
   
    if (hydrate_central_diff_jacobian) then 
      pert(idof) = max(1.d-7 * x(idof),1.d-7)

      x_pert_minus = x
      x_pert_minus(idof) = x(idof) - pert(idof)
      if (option%nflowdof == 3) then
        call HydrateAuxVarCompute3(x_pert_minus, &
               hyd_auxvar(idof+option%nflowdof),global_auxvar,material_auxvar, &
               characteristic_curves,natural_id,option)
      elseif (option%nflowdof == 4) then
        call HydrateAuxVarCompute4(x_pert_minus, &
               hyd_auxvar(idof+option%nflowdof),global_auxvar,material_auxvar, &
               characteristic_curves,natural_id,option)
      endif

    endif

    hyd_auxvar(idof)%pert = pert(idof)
    hyd_auxvar(idof+option%nflowdof)%pert = pert(idof)
    x_pert_plus = x
    x_pert_plus(idof) = x(idof) + pert(idof)
    if (option%nflowdof == 3) then
      call HydrateAuxVarCompute3(x_pert_plus,hyd_auxvar(idof),global_auxvar, &
                                material_auxvar,characteristic_curves, &
                                natural_id,option)
    elseif (option%nflowdof == 4) then
      call HydrateAuxVarCompute4(x_pert_plus,hyd_auxvar(idof),global_auxvar, &
                                material_auxvar,characteristic_curves, &
                                natural_id,option)
    endif
  enddo

  select case(global_auxvar%istate)
    case(L_STATE)
      hyd_auxvar(HYDRATE_LIQUID_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_LIQUID_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
    case(G_STATE)
      hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
      hyd_auxvar(HYDRATE_G_STATE_AIR_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_G_STATE_AIR_PRESSURE_DOF)%pert / &
        HYDRATE_PRESSURE_SCALE
    case(H_STATE)
      hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
    case(I_STATE)
      hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
    case(GA_STATE)
      hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
      if (hydrate_2ph_energy_dof == HYDRATE_AIR_PRESSURE_INDEX) then
        hyd_auxvar(HYDRATE_2PH_STATE_AIR_PRESSURE_DOF)%pert = &
          hyd_auxvar(HYDRATE_2PH_STATE_AIR_PRESSURE_DOF)%pert / &
          HYDRATE_PRESSURE_SCALE
      endif
    case(HG_STATE)
      hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
    case(HA_STATE)
      hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
    case(HI_STATE)
      hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
    case(GI_STATE)
      hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
    case(AI_STATE)
      hyd_auxvar(HYDRATE_LIQUID_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_LIQUID_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
    case(HGA_STATE)
    case(HAI_STATE)
      hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
    case(HGI_STATE)

    case(GAI_STATE)
      hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert = &
        hyd_auxvar(HYDRATE_GAS_PRESSURE_DOF)%pert / HYDRATE_PRESSURE_SCALE
    case(QUAD_STATE)
  end select

end subroutine HydrateAuxVarPerturb


! ************************************************************************** !

subroutine HydratePrintAuxVars(hydrate_auxvar,global_auxvar,material_auxvar, &
                               natural_id,string,option)
  ! 
  ! Prints out the contents of an auxvar
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: hydrate_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  type(option_type) :: option

  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_energy, gas_energy
  PetscReal :: liquid_saturation, gas_saturation

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id

  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_energy = 0.d0
  gas_energy = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0

  print *, '--------------------------------------------------------'
  print *, trim(string)
  print *, '                 cell id: ', natural_id
  select case(global_auxvar%istate)
    case(L_STATE)
      print *, '     Thermodynamic state: Liquid phase'
      liquid_density = hydrate_auxvar%den(lid)
      liquid_energy = hydrate_auxvar%U(lid)
      liquid_saturation = hydrate_auxvar%sat(lid)
    case(G_STATE)
      print *, '     Thermodynamic state: Gas phase'
      gas_density = hydrate_auxvar%den(gid)
      gas_energy = hydrate_auxvar%U(gid)
      gas_saturation = hydrate_auxvar%sat(gid)
    case(GA_STATE)
      print *, '     Thermodynamic state: Two phase'
      liquid_density = hydrate_auxvar%den(lid)
      gas_density = hydrate_auxvar%den(gid)
      liquid_energy = hydrate_auxvar%U(lid)
      gas_energy = hydrate_auxvar%U(gid)
      liquid_saturation = hydrate_auxvar%sat(lid)
      gas_saturation = hydrate_auxvar%sat(gid)
  end select
  liquid_mass = (liquid_density*hydrate_auxvar%xmol(lid,lid)* & 
                 liquid_saturation+ &
                 gas_density*hydrate_auxvar%xmol(lid,gid)* & 
                 gas_saturation)* & 
                 hydrate_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*hydrate_auxvar%xmol(gid,lid)* & 
              liquid_saturation+ &
              gas_density*hydrate_auxvar%xmol(gid,gid)* & 
              gas_saturation)* & 
              hydrate_auxvar%effective_porosity*material_auxvar%volume
  print *, 'tot liq comp mass [kmol]: ', liquid_mass
  print *, 'tot gas comp mass [kmol]: ', gas_mass
  print *, '             energy [MJ]: ', liquid_mass*liquid_energy + &
                                         gas_mass*gas_energy
  print *, '         liquid pressure: ', hydrate_auxvar%pres(lid)
  print *, '            gas pressure: ', hydrate_auxvar%pres(gid)
  print *, '            air pressure: ', hydrate_auxvar%pres(apid)
  print *, '      capillary pressure: ', hydrate_auxvar%pres(cpid)
  print *, '          vapor pressure: ', hydrate_auxvar%pres(vpid)
  print *, '     saturation pressure: ', hydrate_auxvar%pres(spid)
  print *, '       liquid saturation: ', hydrate_auxvar%sat(lid)
  print *, '          gas saturation: ', hydrate_auxvar%sat(gid)
  print *, '   liquid density [kmol]: ', hydrate_auxvar%den(lid)
  print *, '      gas density [kmol]: ', hydrate_auxvar%den(gid)
  print *, '     liquid density [kg]: ', hydrate_auxvar%den_kg(lid)
  print *, '        gas density [kg]: ', hydrate_auxvar%den_kg(gid)
  print *, '         temperature [C]: ', hydrate_auxvar%temp
  print *, '      liquid H [MJ/kmol]: ', hydrate_auxvar%H(lid)
  print *, '         gas H [MJ/kmol]: ', hydrate_auxvar%H(gid)
  print *, '      liquid U [MJ/kmol]: ', hydrate_auxvar%U(lid)
  print *, '         gas U [MJ/kmol]: ', hydrate_auxvar%U(gid)
  print *, '     X (water in liquid): ', hydrate_auxvar%xmol(lid,lid)
  print *, '       X (air in liquid): ', hydrate_auxvar%xmol(gid,lid)
  print *, '        X (water in gas): ', hydrate_auxvar%xmol(lid,gid)
  print *, '          X (air in gas): ', hydrate_auxvar%xmol(gid,gid)
  print *, '         liquid mobility: ', hydrate_auxvar%mobility(lid)
  print *, '            gas mobility: ', hydrate_auxvar%mobility(gid)
  print *, '      effective porosity: ', hydrate_auxvar%effective_porosity
  print *, '--------------------------------------------------------'

end subroutine HydratePrintAuxVars

! ************************************************************************** !

subroutine HydrateOutputAuxVars1(hydrate_auxvar,global_auxvar,material_auxvar, &
                                 natural_id,string,append,option)
  ! 
  ! Prints out the contents of an auxvar to a file
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Global_Aux_module
  use Material_Aux_class
  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: hydrate_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  character(len=*) :: string
  PetscBool :: append
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string2
  PetscInt :: apid, cpid, vpid, spid
  PetscInt :: gid, lid, acid, wid, eid
  PetscReal :: liquid_mass, gas_mass
  PetscReal :: liquid_density, gas_density
  PetscReal :: liquid_energy, gas_energy
  PetscReal :: liquid_saturation, gas_saturation

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id
  spid = option%saturation_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id
  
  liquid_density = 0.d0
  gas_density = 0.d0
  liquid_energy = 0.d0
  gas_energy = 0.d0
  liquid_saturation = 0.d0
  gas_saturation = 0.d0

  write(string2,*) natural_id
  string2 = trim(adjustl(string)) // '_' // trim(adjustl(string2)) // '.txt'
  if (append) then
    open(unit=86,file=string2,position='append')
  else
    open(unit=86,file=string2)
  endif

  write(86,*) '--------------------------------------------------------'
  write(86,*) trim(string)
  write(86,*) '             cell id: ', natural_id
  select case(global_auxvar%istate)
    case(L_STATE)
      write(86,*) ' Thermodynamic state: Liquid phase'
      liquid_density = hydrate_auxvar%den(lid)
      liquid_energy = hydrate_auxvar%U(lid)
      liquid_saturation = hydrate_auxvar%sat(lid)
    case(G_STATE)
      write(86,*) ' Thermodynamic state: Gas phase'
      gas_density = hydrate_auxvar%den(gid)
      gas_energy = hydrate_auxvar%U(gid)
      gas_saturation = hydrate_auxvar%sat(gid)
    case(GA_STATE)
      write(86,*) ' Thermodynamic state: Two phase'
      liquid_density = hydrate_auxvar%den(lid)
      gas_density = hydrate_auxvar%den(gid)
      liquid_energy = hydrate_auxvar%U(lid)
      gas_energy = hydrate_auxvar%U(gid)
      liquid_saturation = hydrate_auxvar%sat(lid)
      gas_saturation = hydrate_auxvar%sat(gid)
  end select
  liquid_mass = (liquid_density*hydrate_auxvar%xmol(lid,lid)* & 
                 liquid_saturation+ &
                 gas_density*hydrate_auxvar%xmol(lid,gid)* & 
                 gas_saturation)* & 
                 hydrate_auxvar%effective_porosity*material_auxvar%volume
  gas_mass = (liquid_density*hydrate_auxvar%xmol(gid,lid)* & 
              liquid_saturation+ &
              gas_density*hydrate_auxvar%xmol(gid,gid)* & 
              gas_saturation)* & 
              hydrate_auxvar%effective_porosity*material_auxvar%volume
  write(86,*) 'tot liq comp mass [kmol]: ', liquid_mass
  write(86,*) 'tot gas comp mass [kmol]: ', gas_mass
  write(86,*) '             energy [MJ]: ', liquid_mass*liquid_energy + &
                                            gas_mass*gas_energy
  write(86,*) '         liquid pressure: ', hydrate_auxvar%pres(lid)
  write(86,*) '            gas pressure: ', hydrate_auxvar%pres(gid)
  write(86,*) '            air pressure: ', hydrate_auxvar%pres(apid)
  write(86,*) '      capillary pressure: ', hydrate_auxvar%pres(cpid)
  write(86,*) '          vapor pressure: ', hydrate_auxvar%pres(vpid)
  write(86,*) '     saturation pressure: ', hydrate_auxvar%pres(spid)
  write(86,*) '         temperature [C]: ', hydrate_auxvar%temp
  write(86,*) '       liquid saturation: ', hydrate_auxvar%sat(lid)
  write(86,*) '          gas saturation: ', hydrate_auxvar%sat(gid)
  write(86,*) '   liquid density [kmol]: ', hydrate_auxvar%den(lid)
  write(86,*) '     liquid density [kg]: ', hydrate_auxvar%den_kg(lid)
  write(86,*) '      gas density [kmol]: ', hydrate_auxvar%den(gid)
  write(86,*) '        gas density [kg]: ', hydrate_auxvar%den_kg(gid)
  write(86,*) '     X (water in liquid): ', hydrate_auxvar%xmol(lid,lid)
  write(86,*) '       X (air in liquid): ', hydrate_auxvar%xmol(gid,lid)
  write(86,*) '        X (water in gas): ', hydrate_auxvar%xmol(lid,gid)
  write(86,*) '          X (air in gas): ', hydrate_auxvar%xmol(gid,gid)
  write(86,*) '      liquid H [MJ/kmol]: ', hydrate_auxvar%H(lid)
  write(86,*) '         gas H [MJ/kmol]: ', hydrate_auxvar%H(gid)
  write(86,*) '      liquid U [MJ/kmol]: ', hydrate_auxvar%U(lid)
  write(86,*) '         gas U [MJ/kmol]: ', hydrate_auxvar%U(gid)
  write(86,*) '         liquid mobility: ', hydrate_auxvar%mobility(lid)
  write(86,*) '            gas mobility: ', hydrate_auxvar%mobility(gid)
  write(86,*) '      effective porosity: ', hydrate_auxvar%effective_porosity
  write(86,*) '...'
  write(86,*) liquid_mass
  write(86,*) gas_mass
  write(86,*) liquid_mass*hydrate_auxvar%U(lid) + &
              gas_mass*hydrate_auxvar%U(gid)
  write(86,*) hydrate_auxvar%pres(lid)
  write(86,*) hydrate_auxvar%pres(gid)
  write(86,*) hydrate_auxvar%pres(apid)
  write(86,*) hydrate_auxvar%pres(cpid)
  write(86,*) hydrate_auxvar%pres(vpid)
  write(86,*) hydrate_auxvar%pres(spid)
  write(86,*) hydrate_auxvar%temp
  write(86,*) hydrate_auxvar%sat(lid)
  write(86,*) hydrate_auxvar%sat(gid)
  write(86,*) hydrate_auxvar%den(lid)
  write(86,*) hydrate_auxvar%den_kg(lid)
  write(86,*) hydrate_auxvar%den(gid)
  write(86,*) hydrate_auxvar%den_kg(gid)
  write(86,*) hydrate_auxvar%xmol(lid,lid)
  write(86,*) hydrate_auxvar%xmol(gid,lid)
  write(86,*) hydrate_auxvar%xmol(lid,gid)
  write(86,*) hydrate_auxvar%xmol(gid,gid)
  write(86,*) hydrate_auxvar%H(lid)
  write(86,*) hydrate_auxvar%H(gid)
  write(86,*) hydrate_auxvar%U(lid)
  write(86,*) hydrate_auxvar%U(gid)
  write(86,*) ''
  write(86,*) hydrate_auxvar%mobility(lid)
  write(86,*) hydrate_auxvar%mobility(gid)
  write(86,*) hydrate_auxvar%effective_porosity
  write(86,*) '--------------------------------------------------------'
  
  close(86)

end subroutine HydrateOutputAuxVars1

! ************************************************************************** !

subroutine HydrateOutputAuxVars2(hydrate_auxvars,global_auxvars,option)
  ! 
  ! Prints out the contents of an auxvar to a file
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Global_Aux_module
  use Option_module

  implicit none

  type(hydrate_auxvar_type) :: hydrate_auxvars(0:,:)
  type(global_auxvar_type) :: global_auxvars(:)
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: apid, cpid, vpid
  PetscInt :: gid, lid, acid, wid, eid
  PetscInt :: i, n, idof

  lid = option%liquid_phase
  gid = option%gas_phase
  apid = option%air_pressure_id
  cpid = option%capillary_pressure_id
  vpid = option%vapor_pressure_id

  acid = option%air_id ! air component id
  wid = option%water_id
  eid = option%energy_id
  
  string = 'hydrate_auxvar.txt'
  open(unit=86,file=string)
  
  n = size(global_auxvars)

100 format(a,100('','',i9))
  
  write(86,'(a,100('','',i9))') '             cell id: ', &
    ((i,i=1,n),idof=0,3)
  write(86,'(a,100('','',i2))') '                idof: ', &
    ((idof,i=1,n),idof=0,3)
  write(86,'(a,100('','',i2))') '               state: ', &
    (global_auxvars(i)%istate,i=1,n)
  write(86,100) '      liquid pressure: ', &
    ((hydrate_auxvars(idof,i)%pres(lid),i=1,n),idof=0,3)
  write(86,100) '         gas pressure: ', &
    ((hydrate_auxvars(idof,i)%pres(gid),i=1,n),idof=0,3)
  write(86,100) '         air pressure: ', &
    ((hydrate_auxvars(idof,i)%pres(apid),i=1,n),idof=0,3)
  write(86,100) '   capillary pressure: ', &
    ((hydrate_auxvars(idof,i)%pres(cpid),i=1,n),idof=0,3)
  write(86,100) '       vapor pressure: ', &
    ((hydrate_auxvars(idof,i)%pres(vpid),i=1,n),idof=0,3)
  write(86,100) '      temperature [C]: ', &
    ((hydrate_auxvars(idof,i)%temp,i=1,n),idof=0,3)
  write(86,100) '    liquid saturation: ', &
    ((hydrate_auxvars(idof,i)%sat(lid),i=1,n),idof=0,3)
  write(86,100) '       gas saturation: ', &
    ((hydrate_auxvars(idof,i)%sat(gid),i=1,n),idof=0,3)
  write(86,100) 'liquid density [kmol]: ', &
    ((hydrate_auxvars(idof,i)%den(lid),i=1,n),idof=0,3)
  write(86,100) '  liquid density [kg]: ', &
    ((hydrate_auxvars(idof,i)%den_kg(lid),i=1,n),idof=0,3)
  write(86,100) '   gas density [kmol]: ', &
    ((hydrate_auxvars(idof,i)%den(gid),i=1,n),idof=0,3)
  write(86,100) '     gas density [kg]: ', &
    ((hydrate_auxvars(idof,i)%den_kg(gid),i=1,n),idof=0,3)
  write(86,100) '  X (water in liquid): ', &
    ((hydrate_auxvars(idof,i)%xmol(lid,lid),i=1,n),idof=0,3)
  write(86,100) '    X (air in liquid): ', &
    ((hydrate_auxvars(idof,i)%xmol(gid,lid),i=1,n),idof=0,3)
  write(86,100) '     X (water in gas): ', &
    ((hydrate_auxvars(idof,i)%xmol(lid,gid),i=1,n),idof=0,3)
  write(86,100) '       X (air in gas): ', &
    ((hydrate_auxvars(idof,i)%xmol(gid,gid),i=1,n),idof=0,3)
  write(86,100) '   liquid H [MJ/kmol]: ', &
    ((hydrate_auxvars(idof,i)%H(lid),i=1,n),idof=0,3)
  write(86,100) '      gas H [MJ/kmol]: ', &
    ((hydrate_auxvars(idof,i)%H(gid),i=1,n),idof=0,3)
  write(86,100) '   liquid U [MJ/kmol]: ', &
    ((hydrate_auxvars(idof,i)%U(lid),i=1,n),idof=0,3)
  write(86,100) '      gas U [MJ/kmol]: ', &
    ((hydrate_auxvars(idof,i)%U(gid),i=1,n),idof=0,3)
  write(86,*)
  write(86,100) '      liquid mobility: ', &
    ((hydrate_auxvars(idof,i)%mobility(lid),i=1,n),idof=0,3)
  write(86,100) '         gas mobility: ', &
    ((hydrate_auxvars(idof,i)%mobility(gid),i=1,n),idof=0,3)
  write(86,100) '   effective porosity: ', &
    ((hydrate_auxvars(idof,i)%effective_porosity,i=1,n),idof=0,3)
  
  close(86)

end subroutine HydrateOutputAuxVars2

! ************************************************************************** !

subroutine HydrateSetBogusDeriv(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  implicit none

  type(hydrate_auxvar_type) :: auxvar

  auxvar%d%pc_satg = 1234.56789d0
  auxvar%d%por_p = 1234.56789d0
  auxvar%d%denl_pl = 1234.56789d0
  auxvar%d%denl_T = 1234.56789d0
  auxvar%d%deng_pg = 1234.56789d0
  auxvar%d%deng_pa = 1234.56789d0
  auxvar%d%deng_T = 1234.56789d0
  auxvar%d%dengkg_pg = 1234.56789d0
  auxvar%d%dengkg_T = 1234.56789d0
  auxvar%d%Ul_pl = 1234.56789d0
  auxvar%d%Ul_T = 1234.56789d0
  auxvar%d%Hl_pl = 1234.56789d0
  auxvar%d%Hl_T = 1234.56789d0
  auxvar%d%Ug_pg = 1234.56789d0
  auxvar%d%Ug_pa = 1234.56789d0
  auxvar%d%Ug_T = 1234.56789d0
  auxvar%d%Hg_pg = 1234.56789d0
  auxvar%d%Hg_pa = 1234.56789d0
  auxvar%d%Hg_T = 1234.56789d0

  auxvar%d%Hv = 1234.56789d0
  auxvar%d%Ha = 1234.56789d0
  auxvar%d%Uv = 1234.56789d0
  auxvar%d%Ua = 1234.56789d0
  auxvar%d%Hv_pg = 1234.56789d0
  auxvar%d%Hv_pa = 1234.56789d0
  auxvar%d%Hv_T = 1234.56789d0
  auxvar%d%Ha_pg = 1234.56789d0
  auxvar%d%Ha_pa = 1234.56789d0
  auxvar%d%Ha_T = 1234.56789d0
  auxvar%d%Uv_pg = 1234.56789d0
  auxvar%d%Uv_pa = 1234.56789d0
  auxvar%d%Uv_T = 1234.56789d0
  auxvar%d%Ua_pg = 1234.56789d0
  auxvar%d%Ua_pa = 1234.56789d0
  auxvar%d%Ua_T = 1234.56789d0
  auxvar%d%denv = 1234.56789d0
  auxvar%d%dena = 1234.56789d0
  auxvar%d%denv_T = 1234.56789d0
  auxvar%d%dena_T = 1234.56789d0
  auxvar%d%denv_pg = 1234.56789d0
  auxvar%d%dena_pg = 1234.56789d0
  auxvar%d%Hc = 1234.56789d0
    
  auxvar%d%psat_p = 1234.56789d0
  auxvar%d%psat_T = 1234.56789d0
  auxvar%d%pv_p = 1234.56789d0
  auxvar%d%pv_pa = 1234.56789d0
  auxvar%d%pv_T = 1234.56789d0
  auxvar%d%Hc_p = 1234.56789d0
  auxvar%d%Hc_T = 1234.56789d0
  auxvar%d%mobilityl_pl = 1234.56789d0
  auxvar%d%mobilityl_T = 1234.56789d0
  auxvar%d%mobilityl_satg = 1234.56789d0
  auxvar%d%mobilityg_pg = 1234.56789d0
  auxvar%d%mobilityg_T = 1234.56789d0
  auxvar%d%mobilityg_satg = 1234.56789d0
  auxvar%d%mobilityg_pa = 1234.56789d0
  auxvar%d%mug = 1234.56789d0
  auxvar%d%mug_T = 1234.56789d0
  auxvar%d%mug_pg = 1234.56789d0
  auxvar%d%xmol_p(:,:) = 1234.56789d0
  auxvar%d%xmol_T(:,:) = 1234.56789d0

end subroutine HydrateSetBogusDeriv

! ************************************************************************** !

subroutine HydrateZeroAnalyticalDeriv(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  implicit none

  type(hydrate_auxvar_type) :: auxvar

  auxvar%d%pc_satg = 0.d0
  auxvar%d%por_p = 0.d0
  auxvar%d%denl_pl = 0.d0
  auxvar%d%denl_T = 0.d0
  auxvar%d%deng_pg = 0.d0
  auxvar%d%deng_pa = 0.d0
  auxvar%d%deng_T = 0.d0
  auxvar%d%dengkg_pg = 0.d0
  auxvar%d%dengkg_T = 0.d0
  auxvar%d%Ul_pl = 0.d0
  auxvar%d%Ul_T = 0.d0
  auxvar%d%Hl_pl = 0.d0
  auxvar%d%Hl_T = 0.d0
  auxvar%d%Ug_pg = 0.d0
  auxvar%d%Ug_pa = 0.d0
  auxvar%d%Ug_T = 0.d0
  auxvar%d%Hg_pg = 0.d0
  auxvar%d%Hg_pa = 0.d0
  auxvar%d%Hg_T = 0.d0

  auxvar%d%Hv = 0.d0
  auxvar%d%Ha = 0.d0
  auxvar%d%Uv = 0.d0
  auxvar%d%Ua = 0.d0
  auxvar%d%Hv_pg = 0.d0
  auxvar%d%Hv_pa = 0.d0
  auxvar%d%Hv_T = 0.d0
  auxvar%d%Ha_pg = 0.d0
  auxvar%d%Ha_pa = 0.d0
  auxvar%d%Ha_T = 0.d0
  auxvar%d%Uv_pg = 0.d0
  auxvar%d%Uv_pa = 0.d0
  auxvar%d%Uv_T = 0.d0
  auxvar%d%Ua_pg = 0.d0
  auxvar%d%Ua_pa = 0.d0
  auxvar%d%Ua_T = 0.d0
  auxvar%d%denv = 0.d0
  auxvar%d%dena = 0.d0
  auxvar%d%denv_T = 0.d0
  auxvar%d%dena_T = 0.d0
  auxvar%d%denv_pg = 0.d0
  auxvar%d%dena_pg = 0.d0
  auxvar%d%Hc = 0.d0
    
  auxvar%d%psat_p = 0.d0
  auxvar%d%psat_T = 0.d0
  auxvar%d%pv_p = 0.d0
  auxvar%d%pv_pa = 0.d0
  auxvar%d%pv_T = 0.d0
  auxvar%d%Hc_p = 0.d0
  auxvar%d%Hc_T = 0.d0
  auxvar%d%mobilityl_pl = 0.d0
  auxvar%d%mobilityl_T = 0.d0
  auxvar%d%mobilityl_satg = 0.d0
  auxvar%d%mobilityg_pg = 0.d0
  auxvar%d%mobilityg_T = 0.d0
  auxvar%d%mobilityg_satg = 0.d0
  auxvar%d%mobilityg_pa = 0.d0
  auxvar%d%mug = 0.d0
  auxvar%d%mug_T = 0.d0
  auxvar%d%mug_pg = 0.d0
  auxvar%d%xmol_p(:,:) = 0.d0
  auxvar%d%xmol_T(:,:) = 0.d0

end subroutine HydrateZeroAnalyticalDeriv

! ************************************************************************** !

subroutine HydrateAuxVarSingleDestroy(auxvar)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  implicit none

  type(hydrate_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call HydrateAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)  

end subroutine HydrateAuxVarSingleDestroy

! ************************************************************************** !

subroutine HydrateAuxVarArray1Destroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  implicit none

  type(hydrate_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call HydrateAuxVarStrip(auxvars(iaux))
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)  

end subroutine HydrateAuxVarArray1Destroy

! ************************************************************************** !

subroutine HydrateAuxVarArray2Destroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  implicit none

  type(hydrate_auxvar_type), pointer :: auxvars(:,:)
  
  PetscInt :: iaux, idof
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars,2)
      do idof = 1, size(auxvars,1)
        call HydrateAuxVarStrip(auxvars(idof-1,iaux))
      enddo
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)  

end subroutine HydrateAuxVarArray2Destroy

! ************************************************************************** !

subroutine HydrateAuxVarStrip(auxvar)
  ! 
  ! HydrateAuxVarDestroy: Deallocates a hydrate auxiliary object
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  type(hydrate_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%pres)  
  call DeallocateArray(auxvar%sat)  
  call DeallocateArray(auxvar%den)  
  call DeallocateArray(auxvar%den_kg)  
  call DeallocateArray(auxvar%xmol)  
  call DeallocateArray(auxvar%H)  
  call DeallocateArray(auxvar%U)  
  call DeallocateArray(auxvar%mobility)  
  if (associated(auxvar%d)) then
    deallocate(auxvar%d)
    nullify(auxvar%d)
  endif
  
end subroutine HydrateAuxVarStrip

! ************************************************************************** !

subroutine HydrateAuxDestroy(aux)
  ! 
  ! Deallocates a hydrate auxiliary object
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  type(hydrate_type), pointer :: aux
  
  if (.not.associated(aux)) return
  
  call HydrateAuxVarDestroy(aux%auxvars)
  call HydrateAuxVarDestroy(aux%auxvars_bc)
  call HydrateAuxVarDestroy(aux%auxvars_ss)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

  if (associated(aux%hydrate_parameter)) then
    call DeallocateArray(aux%hydrate_parameter%diffusion_coefficient)
    deallocate(aux%hydrate_parameter)
  endif
  nullify(aux%hydrate_parameter)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine HydrateAuxDestroy

! ************************************************************************** !

subroutine HydrateCompositeThermalCond(phi,sat,kdry,kwet,keff)

  implicit none

  PetscReal :: phi, kdry, kwet
  PetscReal, pointer :: sat(:)

  PetscReal :: keff
  PetscReal :: k_h20,k_ch4,k_hyd,k_ice
  PetscInt :: lid, gid, hid, iid

  lid = 1
  gid = 2
  hid = 3
  iid = 4

  k_h20 = 0.59d0 !W/m-K
  k_ch4 = 30.d-3 !W/m-K
  k_hyd = 0.58d0 !W/m-K
  k_ice = 2.2d0   !W/m-K

  ! kdry = (1-phi)*k_rock + phi*k_air, k_air << k_rock
  keff = kdry + phi * (sat(lid)*k_h20 + sat(hid)*k_hyd + sat(iid)*k_ice &
          + sat(gid)*k_ch4)

  ! IGHCC2 function (seems odd if phi = 1 and sat(lid) = 1)
  !keff = kdry + phi * (sat(lid)*kwet + sat(hid)*k_hyd + sat(iid) * k_ice &
  !        + sat(gid)*k_ch4)

end subroutine HydrateCompositeThermalCond

! ************************************************************************** !

subroutine HydratePE(T,sat,xmol,PE,dP,characteristic_curves,material_auxvar, &
                     option, hydrate_guest)

  !This subroutine calculates the 3-phase equilibrium pressure of methane
  !hydrate in pure water, from polynomial fit (Moridis, 2003)
  !
  !Author: Michael Nole
  !Date: 01/22/19
  !
  !Author: David Fukuyama (CO2 equilibrium pressure, Kamath, 1984)
  !Date: 09/23/19

  use Characteristic_Curves_module
  use Option_module
  use Material_Aux_class

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: sat
  PetscReal, intent(in) :: xmol
  PetscInt, intent(in) :: hydrate_guest
  PetscReal, intent(out) :: PE
  PetscReal, intent(out) :: dP

  class(characteristic_curves_type) :: characteristic_curves
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option

  PetscReal :: T_temp, dTf, dTfs

  if (hydrate_with_gibbs_thomson .and. hydrate_salinity) then
    call GibbsThomsonFreezing(1.d0-sat, 54734.d0, HYDRATE_DENSITY, T, dTf, &
         characteristic_curves, material_auxvar, option)
    call SalinityEffects(xmol,dTfs,option)
    dTf = dTf - dTfs
  elseif (hydrate_salinity) then
     call SalinityEffects(xmol,dTf,option)
  elseif (hydrate_with_gibbs_thomson) then
     call GibbsThomsonFreezing(1.d0-sat, 54734.d0, HYDRATE_DENSITY, T, dTf, &
          characteristic_curves, material_auxvar, option)
  else
    dTf = 0.d0
  endif

  T_temp = T + 273.15d0 + dTf
  select case(hydrate_guest)
    case(HYDRATE_FORMER_CO2)
      !Kamath, 1984
      !if (T < TQD) then
      !  PE = exp(18.594 - 3161.41 / T_temp)*1.d-3
      !else
      !  PE = exp(44.580 - 10246.28 / T_temp)*1.d-3
      !endif

      !Sloan compilation fit (Clathrate Hydrates of Natural Gases)
      if (T < TQD) then
        PE = 1.1046 + 0.04449 * (T_temp - 273.15) + 0.000629 * (T_temp - 273.15) ** 2
      else
        PE = 1.2241 + 0.13700 * (T_temp - 273.15) ** 2 - 0.0015018 * (T_temp - &
              273.15) ** 3 + 0.0001733 * (T_temp - 273.15) ** 4
     endif

    case(HYDRATE_FORMER_CH4)

      dP = 0.d0

      if (T < TQD) then
        select case(hydrate_phase_boundary)
          case(1)
            !Kamath, 1984. The reason this is default is because it is
            !              invertible.
            PE = exp(1.4717d1-1.88679d3/T_temp)*1.d-3
          case(2)
            !Moridis, 2003
            PE = exp(-43.8921173434628 + 0.776302133739303 * T_temp &
               - 7.27291427030502d-3 * T_temp**2 + 3.85413985900724d-5 * &
               T_temp**3 - 1.03669656828834d-7 * T_temp**4 + &
               1.09882180475307d-10 * T_temp**5)
        end select
      else
        select case(hydrate_phase_boundary)
          case(1)
            !Kamath, 1984
            PE = exp(3.898d1-8.533d3/T_temp)*1.d-3
            if (hydrate_adjust_ghsz_solubility) then
              dP = PE - exp(3.898d1 - 8.533d3/(T_temp - dTf))* 1.d-3
              dP = dP * 1.d6
            endif
          case(2)
            !Moridis, 2003
            PE = exp(-1.9413850446456d5 + 3.31018213397926d3 * T_temp &
               - 22.5540264493806* T_temp**2 + 0.0767559117787059 * T_temp**3 &
               - 1.30465829788791d-4 * T_temp**4 + 8.86065316687571d-8 * &
               T_temp**5)
            if (hydrate_adjust_ghsz_solubility) then
              dP = PE - exp(-1.9413850446456d5 + 3.31018213397926d3 * &
                 (T_temp-dTf) - 22.5540264493806*(T_temp-dTf)**2 + &
                 0.0767559117787059 * (T_temp-dTf)**3 - &
                 1.30465829788791d-4 * (T_temp-dTf)**4 + 8.86065316687571d-8 * &
                 (T_temp-dTf)**5)
              dP = dP * 1.d6
            endif
        end select
      endif
  end select
  PE = PE * 1.d6

end subroutine HydratePE

! ************************************************************************** !
subroutine HydrateMethanogenesis(z,offset,hydrate_parameter,q_meth)

  ! A simple methanogenesis source parameterized as a function of depth
  ! assuming top of domain is the seafloor.
  ! Author: Michael Nole
  ! Date: 03/05/19
  !

  implicit none

  PetscReal :: z, offset
  type(hydrate_parameter_type), pointer :: hydrate_parameter
  PetscReal :: q_meth

  type(methanogenesis_type), pointer :: methanogenesis
  PetscReal :: alpha, k_alpha, lambda, omega, z_smt

  methanogenesis => hydrate_parameter%methanogenesis

  alpha = methanogenesis%alpha
  k_alpha = methanogenesis%k_alpha
  lambda = methanogenesis%lambda
  omega = methanogenesis%omega
  z_smt = methanogenesis%z_smt

  if (offset - z > z_smt) then
    ! Malinverno, 2011
    q_meth = k_alpha * lambda * alpha * exp(-lambda/omega * (offset - &
                                    z - z_smt))
  else
    q_meth = 0.d0
  endif

  !kg/m^3/s to kmol/s
  q_meth = q_meth / MW_CH4

end subroutine HydrateMethanogenesis
! ************************************************************************** !

subroutine HenrysConstant(T,hydrate_guest,K_H)

  !Calculates the Henry's constant as a function of temperature
  !
  !
  !Author: Michael Nole
  !Date: 01/22/19
  !

  use EOS_Water_module

  implicit none

  PetscReal, intent(in) :: T
  PetscInt, intent(in) :: hydrate_guest
  PetscReal, intent(out) :: K_H

  PetscReal :: T_temp, T_star, Pv
  PetscErrorCode :: ierr
  PetscReal, parameter :: Tc = 647.14d0

  T_temp = T + 273.15d0

  select case(hydrate_guest)
    case(HYDRATE_FORMER_CO2)
      !Harvey, 1996
      T_star=T_temp/Tc
      call EOSWaterSaturationPressure(T,Pv,ierr)
      K_H = exp(log(Pv) -9.4234d0/(T_star) + 4.0087d0*(1.d0-T_star)** &
            (0.355d0)/(T_star) + 10.3199d0 * exp(1.d0-T_star)*T_star**(-0.41d0))
      !K_H = K_H * 1.d6
    case(HYDRATE_FORMER_CH4)
      select case(hydrate_henrys_constant)
        case(1)
          !Carroll & Mather  Units: Pa/mol frac
          K_H = exp(5.1345 + 7837.d0/T_temp - 1.509d6/(T_temp**2) + 2.06d7/ &
                (T_temp**3)) *1.d3
        case(2)
          !Cramer, 1982
          K_H = 1.d5*(24582.4d0 + 6.71091d2*T + 6.87067d0*T**2 - &
                  1.773079d-1*T**3 + 1.09652d-03*T**4 - &
                  3.19599d-6*T**5 + 4.46172d-9*T**6 - &
                  2.40294d-12*T**7)
      end select

    case default
      !Carroll & Mather  Units: Pa/mol frac
      K_H = exp(5.1345 + 7837.d0/T_temp - 1.509d6/(T_temp**2) + 2.06d7/ &
            (T_temp**3)) *1.d3
  end select
end subroutine HenrysConstant

subroutine HydrateGHSZSolubilityCorrection(T,P,dP,K_H)

  !Adjusts methane solubility within the hydrate stabilty zone, following
  !Davie et al., 2004
  !
  !Author: Michael Nole
  !

  implicit none

  PetscReal, intent(in) :: T, P
  PetscReal :: K_H, dP

  PetscReal, parameter :: C3_0 = 156.36d0 !mM
  PetscReal, parameter :: T_0 = 292.d0 !K
  PetscReal, parameter :: P_0 = 20.d0 !MPa
  PetscReal, parameter :: dC3_dT = 6.34d0 !mM
  PetscReal, parameter :: dC3_dP = 1.11d0 !mM
  PetscReal, parameter :: alpha = 14.4d0 !C
  PetscReal :: logP

  PetscReal :: T3


  ! Inverting the Moridis equation
  if (T > TQD) then
    select case (hydrate_phase_boundary)
      case(1)
        !Kamath
        T3 = -8.533d3/(log((P-dP)/1.d6*1.d3)-3.898d1)
      case(2)
        !Moridis
        !Lower-order
        T3 = 9.0622d0 * log((P-dP)*1.d-6) + 264.66d0

        !Higher-order
        !logP = log((P-dP)*1.d-6)
        !T3 = -1.09874018d-2*logP**6 + 1.73301550d-1*logP**5 - &
        !     9.67897401d-1*logP**4 + 2.34919362d0*logP**3  - &
        !     2.77146625d0*logP**2 + 1.13889445d1*logP + 2.63495959d2

    end select
  else
    T3 = T + 273.15d0
  endif

  K_H = K_H / exp((T+273.15d0-T3)/alpha)

end subroutine HydrateGHSZSolubilityCorrection

! ************************************************************************** !

subroutine HydrateScalePermeability(perm_base,permeability,solid_sat_eff,option)
  !
  ! Compute permeability as a function of solid saturation.
  !
  ! Author: Michael Nole
  ! Date: 08/26/20
  !

  use Material_Aux_class
  use Option_module

  implicit none

  PetscReal :: perm_base
  PetscReal :: permeability(3)
  PetscReal :: solid_sat_eff
  type(option_type) :: option

  select case (hydrate_perm_scaling_function)
    case(1) ! Dai and Seol, 2014
      if (perm_base < -999.d0) then
        perm_base = permeability(1)
      endif
      permeability(:) = perm_base * &
                        (1.d0-solid_sat_eff)**3/(1.d0+2.d0*solid_sat_eff)**2
    case default
      option%io_buffer = 'Permeability scaling function not recognized.'
      call PrintErrMsgByRank(option)
  end select

end subroutine HydrateScalePermeability

! ************************************************************************** !

subroutine GibbsThomsonFreezing(sat,Hf,rho,Tb,dTf,characteristic_curves,&
                                material_auxvar,option)

  !This subroutine ties the capillary pressure function to a Gibbs-Thomson
  !subcooling required to precipitate a solid in pores.
  !
  !Author: Michael Nole
  !Date: 04/04/19
  !

  use Characteristic_Curves_module
  use Option_module
  use Material_Aux_class

  implicit none

  PetscReal, intent(in) :: sat
  PetscReal, intent(in) :: Hf
  PetscReal, intent(in) :: rho
  PetscReal, intent(in) :: Tb
  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  class(material_auxvar_type) :: material_auxvar
  PetscReal, intent(out) :: dTf

  PetscReal :: Pc,dpc_dsatl,sigma,theta

  sigma = 0.073d0
  theta = 0.d0

  !if (material_auxvar%pore_size < 0.d0) then
    call characteristic_curves%saturation_function% &
             CapillaryPressure(sat,Pc,dpc_dsatl,option)
    dTf = (Tb+273.15)*Pc/(Hf * rho * 1000.d0)
  !else
  !  dTf = (Tb+273.15)*2*sigma*cos(theta)/(Hf * rho * 1000.d0 * &
  !          material_auxvar%pore_size)
  !endif


  !MAN debugging
  !dTf = 0.d0
end subroutine GibbsThomsonFreezing

! ************************************************************************** !
subroutine SalinityEffects(xmol,dTd,option)
  !
  ! This ties salinity to the subcooling required to precipitate
  ! hydrate in pores, similar to the GibbsThomsonFreezing subroutine
  ! Based on the TOUGH method of PE-salinity behavior, using reference
  ! parameters from Dickens and Quinby-Hunt, 1997
  !
  ! Author: David Fukuyama
  ! Date: 5/1/2020
  !

  use Characteristic_Curves_module
  use Option_module
  use Material_Aux_class

  implicit none

  PetscReal, intent(in) :: xmol
  type(option_type) :: option
  PetscReal, intent(out) :: dTd

  PetscReal :: dTr = -0.37
  PetscReal :: xmol_ref = 0.03311

  dTd = dTr * log(1-xmol)/log(1-xmol_ref)
end subroutine salinityEffects

! ************************************************************************** !


subroutine EOSIceEnergy(T,U)

  !Internal energy of ice as f(Temperature) (Fukusako and Yamamoto, 1993)
  !
  !Author: Michael Nole
  !Date: 04/04/19
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: U
  PetscReal, parameter :: Lw = -6017.1d0 !Latent heat of fusion,  J/mol
  PetscReal :: T_temp

  T_temp = T + 273.15d0

  if (T_temp >= 90.d0) then
    U = Lw + 185.d0 * (T_temp-273.15d0) + 3.445 * (T_temp**2 - 273.15d0**2)
  else
    U = Lw + 4.475 * (T_temp**2 - 273.15d0**2)
  endif

  ! J/mol to MJ/kmol
  U = U / 1.d3

end subroutine EOSIceEnergy

! ************************************************************************** !

subroutine EOSHydrateEnthalpy(T,H,hydrate_type)

  !Enthalpy of gas hydrate as f(Temperature) (Handa, 1986)
  !
  !Author: Michael Nole
  !Date: 01/22/19
  !
  !Author: David Fukuyama (CO2 hydrate enthalpy)
  !Date: 09/23/19

  use Option_module

  implicit none


  PetscReal, intent(in):: T
  PetscInt, intent(in):: hydrate_type
  PetscReal, intent(out) :: H

  PetscReal, parameter :: Hh0 = -54734.d0 ! J/mol
  PetscReal, parameter :: Hh0_c = -66800.d0  !J/mol
  PetscReal :: Cph, T_temp
!  PetscReal :: hydrate_type
  T_temp = T + 273.15d0
  select case(hydrate_type)
    case(HYDRATE_FORMER_CO2)
      !Constant Cp - this is wrong- just to test
      !Cph = 1620.d0 * (MW_H2O * Nhyd + MW_CO2)/1.d3
      !H = Cph * (T - TQD) + Hh0 / (Nhyd + 1.d0)
      !H = H / 1.d3

      !Cph = 1.7 * T_temp + 2098
      !H = Cph * (T - TQD)

      Cph = 1620.d0*(MW_H2O*Nhyd + MW_CO2)/1.d3
      H = Cph * (T-TQD) + Hh0_c / (Nhyd + 1.d0)
      H = H / 1.d3

    case(HYDRATE_FORMER_CH4)

      ! Integral of Cph * dT ; Cph from Handa, 1998

      ! Units: J/mol
      !H = Hh0 + 6.6d0 * (T_temp-273.15d0) + 7.269d-1 * (T_temp**2 - 273.15d0**2) - 1.21333d-3 * &
      !      (T_temp**3 - 273.15d0**3)  + 1.578d-6 * (T_temp**4 - 273.15d0**4)
      ! Units: MJ/kmol
      !H = H / 1.d3

      !H = H / (Nhyd+1.d0)

      !Constant Cp
      Cph = 1620.d0*(MW_H2O*Nhyd + MW_CH4)/1.d3
      H = Cph * (T-TQD) + Hh0 / (Nhyd + 1.d0)
      H = H / 1.d3
    case default

      ! Integral of Cph * dT ; Cph from Handa, 1998

      ! Units: J/mol
      !H = Hh0 + 6.6d0 * (T_temp-273.15d0) + 7.269d-1 * (T_temp**2 - 273.15d0**2) - 1.21333d-3 * &
      !      (T_temp**3 - 273.15d0**3)  + 1.578d-6 * (T_temp**4 - 273.15d0**4)
      ! Units: MJ/kmol
      !H = H / 1.d3

      !H = H / (Nhyd+1.d0)

      !Constant Cp
      Cph = 1620.d0*(MW_H2O*Nhyd + MW_CH4)/1.d3
      H = Cph * (T-TQD) + Hh0 / (Nhyd + 1.d0)
      H = H / 1.d3
  end select
end subroutine EOSHydrateEnthalpy

! ************************************************************************** !

end module Hydrate_Aux_module
