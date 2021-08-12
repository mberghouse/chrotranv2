module Reaction_Gas_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Aux_module
  use Reactive_Transport_Aux_module  
  use Global_Aux_module
  use Reaction_Gas_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
 
  private

  
  public :: RGasRead, &
            RTotalGas, &
            RTotalSorbGas, &
            RAccumulationSorbGas, &
            RTotalCO2

contains

! ************************************************************************** !

subroutine RGasRead(gas_species_list,gas_type,error_msg,input,option)
  ! 
  ! Reads immobile species
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13/ 08/01/16
  ! 
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(gas_species_type), pointer :: gas_species_list
  PetscInt :: gas_type
  character(len=MAXSTRINGLENGTH) :: error_msg
  type(input_type), pointer :: input
  type(option_type) :: option
  
  type(gas_species_type), pointer :: new_gas_species, &
                                     prev_gas_species

  ! since both active and passive gases are in the same list, skip to the
  ! end of the list if it exists.
  if (associated(gas_species_list)) then
    prev_gas_species => gas_species_list
    do
      if (.not.associated(prev_gas_species%next)) exit
      prev_gas_species => prev_gas_species%next
    enddo
  else
    nullify(prev_gas_species)
  endif
  ! read in new gases
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    new_gas_species => GasSpeciesCreate()
    call InputReadCard(input,option,new_gas_species%name)  
    call InputErrorMsg(input,option,'keyword',error_msg)    
    new_gas_species%itype = gas_type
    if (associated(prev_gas_species)) then
      prev_gas_species%next => new_gas_species
      new_gas_species%id = prev_gas_species%id + 1
    else
      gas_species_list => new_gas_species
      new_gas_species%id = 1
    endif
    prev_gas_species => new_gas_species
    nullify(new_gas_species)
  enddo         
  call InputPopBlock(input,option)
                                          
end subroutine RGasRead

! ************************************************************************** !

subroutine RTotalGas(rt_auxvar,global_auxvar,reaction,option)
  ! 
  ! Computes the total component concentrations and derivative with
  ! respect to free-ion
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/01/16
  ! 

  use Option_module

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option
  
  PetscInt, parameter :: iphase = 2
  PetscInt :: i, j, igas, icomp, jcomp, ncomp
  PetscReal :: ln_conc(reaction%naqcomp)
  PetscReal :: ln_act(reaction%naqcomp)
  PetscReal :: lnQK, tempreal
  PetscReal :: RT
  PetscReal :: gas_concentration
  type(gas_type), pointer :: gas
  
  rt_auxvar%total(:,iphase) = 0.d0 !debugging 
  
  gas => reaction%gas
  ! units of ideal gas constant = J/mol-K = kPa-L/mol-K
  ! units of RT = Pa-L/mol
  RT = IDEAL_GAS_CONSTANT*(global_auxvar%temp+273.15d0)*1.d3
  
  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)
  
  ! initialize derivatives
  rt_auxvar%aqueous%dtotal(:,:,iphase) = 0.d0
   
  do igas = 1, gas%nactive_gas ! for each secondary species
    ! compute secondary species concentration
    lnQK = -gas%acteqlogK(igas)*LOG_TO_LN

    ! activity of water
    if (gas%acteqh2oid(igas) > 0) then
      lnQK = lnQK + gas%acteqh2ostoich(igas)*rt_auxvar%ln_act_h2o
    endif

    ncomp = gas%acteqspecid(0,igas)
    do i = 1, ncomp
      icomp = gas%acteqspecid(i,igas)
      lnQK = lnQK + gas%acteqstoich(i,igas)*ln_act(icomp)
    enddo
    ! units = bars
    rt_auxvar%gas_pp(igas) = exp(lnQK)
    ! unit = mol/L gas
    gas_concentration = rt_auxvar%gas_pp(igas) * 1.d5 / RT

    ! add contribution to primary totals
    ! units of total = mol/L gas
    do i = 1, ncomp
      icomp = gas%acteqspecid(i,igas)
      rt_auxvar%total(icomp,iphase) = rt_auxvar%total(icomp,iphase) + &
                                      gas%acteqstoich(i,igas)* &
                                      gas_concentration
    enddo
    
    ! add contribution to derivatives of total with respect to free
    ! units of dtotal = kg water / L gas
    do j = 1, ncomp
      jcomp = gas%acteqspecid(j,igas)
      tempreal = gas%acteqstoich(j,igas)*exp(lnQK-ln_conc(jcomp))*1.d5/RT
      do i = 1, ncomp
        icomp = gas%acteqspecid(i,igas)
        rt_auxvar%aqueous%dtotal(icomp,jcomp,iphase) = &
          rt_auxvar%aqueous%dtotal(icomp,jcomp,iphase) + &
          gas%acteqstoich(i,igas)*tempreal
      enddo
    enddo
  enddo

end subroutine RTotalGas

! ************************************************************************** !

subroutine RTotalSorbGas(rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         isotherm_rxn,option)
  !
  ! Computes the total sorbed component partial pressures in the gas phase
  ! 
  ! Author: Michael Nole
  ! Date: 06/23/21
  !

  use Option_module
  use Material_Aux_class
  use Reaction_Isotherm_Aux_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(reaction_rt_type) :: reaction
  type(isotherm_rxn_type) :: isotherm_rxn
  type(option_type) :: option
  type(gas_type), pointer :: gas

  gas => reaction%gas

  if (gas%isotherm%neqkdrxn > 0) then
      call RTotalSorbGasKD(rt_auxvar,global_auxvar,material_auxvar,gas, &
                           gas%isotherm%isotherm_rxn,option)
  endif

end subroutine RTotalSorbGas

! ************************************************************************** !

subroutine RTotalSorbGasKD(rt_auxvar,global_auxvar,material_auxvar,gas, &
                        isotherm_rxn,option)
  ! 
  ! Computes the total sorbed component concentrations
  ! for the K_D model
  ! 
  ! Author: Michael Nole
  ! Date: 06/23/21
  ! 

  use Option_module

  use Global_Aux_module
  use Material_Aux_class
  use Reaction_Isotherm_module

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(gas_type) :: gas
  type(isotherm_rxn_type) :: isotherm_rxn
  type(option_type) :: option

  type(isotherm_type), pointer :: isotherm
  PetscInt :: irxn,i
  PetscInt :: icomp
  PetscReal :: res
  PetscReal :: dres_dc
  PetscReal :: molality
  PetscReal :: tempreal
  PetscReal :: one_over_n
  PetscReal :: molality_one_over_n
  PetscReal :: kd_kgw_m3b, kd_sorb_gas, rf
  PetscReal :: partial_pres,pp_one_over_n
  PetscReal :: gas_concentration
  PetscReal :: RT
  PetscReal :: L_pore, L_gas
  PetscInt :: ncomp

  PetscInt, parameter :: iphase = 2

  isotherm => gas%isotherm

  ! units of ideal gas constant = J/mol-K = kPa-L/mol-K
  ! units of RT = Pa-L/mol
  RT = IDEAL_GAS_CONSTANT*(global_auxvar%temp+273.15d0)*1.d3
  ! Concentration units = mol/L gas
  gas_concentration = partial_pres * 1.d5 / RT
  ! Pore volume (L)
  L_pore = material_auxvar%porosity*material_auxvar%volume*1.d3
  ! Gas volume (L)
  L_gas = L_pore*global_auxvar%sat(option%gas_phase)
  
  do irxn = 1, isotherm%neqkdrxn
    partial_pres =  rt_auxvar%gas_pp(irxn)
    ncomp = gas%acteqspecid(0,irxn)
    if (isotherm%ikd_units == KD_UNIT_MLW_GSOIL) then
                   !KD units [mL water/g soil]
      kd_kgw_m3b = isotherm_rxn%eqisothermcoeff(irxn) * &
                   global_auxvar%den_kg(iphase) * &
                   (1.d0-material_auxvar%porosity) * &
                   material_auxvar%soil_particle_density * &
                   1.d-3 ! convert mL water/g soil to m^3 water/kg soil

    else
      ! kd_unit = KD_UNIT_KG_M3_BULK
      kd_kgw_m3b = isotherm_rxn%eqisothermcoeff(irxn)
    endif
    if (isotherm%eqkdmineral(irxn) > 0) then
      ! NOTE: mineral volume fraction here is solely a scaling factor.  It has 
      ! nothing to do with the soil volume; that is calculated through as a 
      ! function of porosity.
      kd_kgw_m3b = isotherm_rxn%eqisothermcoeff(irxn) * &
                   (rt_auxvar%mnrl_volfrac(isotherm%eqkdmineral(irxn)))
    endif
    select case(isotherm%eqisothermtype(irxn))
      case(SORPTION_LINEAR)
        ! Csorb = Kd*Pg
        res = kd_kgw_m3b*partial_pres
        dres_dc = kd_kgw_m3b
      case(SORPTION_LANGMUIR)
        ! Csorb = K*Caq*b/(1+K*Pg)
        tempreal = kd_kgw_m3b*partial_pres
        res = tempreal*isotherm_rxn%eqisothermlangmuirb(irxn) / &
              (1.d0 + tempreal)
        dres_dc = res/partial_pres - &
                  res / (1.d0 + tempreal) * tempreal / partial_pres
      case(SORPTION_FREUNDLICH)
        ! Csorb = Kd*Pg**(1/n)
        one_over_n = 1.d0/isotherm_rxn%eqisothermfreundlichn(irxn)
        pp_one_over_n = partial_pres**one_over_n
        res = kd_kgw_m3b*partial_pres**one_over_n
        dres_dc = res/partial_pres*one_over_n
      case(SORPTION_RETENTION_FACTOR)
        ! Retention factor: moles sorbed/moles gas
        ! (1/(1+1/rf)) = moles sorbed / (moles sorbed + moles in gas phase)
        ! units = mole solute/L gas * L gas * sorbed/total / m^3 bulk
        !       = mole solute sorbed / m^3 bulk
        rf = isotherm_rxn%eqisothermretentionfactor(irxn)
        res = gas_concentration * L_gas * (1/(1+1/rf)) / material_auxvar%volume
        !res = gas_concentration*L_gas*kd_sorb_gas / material_auxvar%volume
        dres_dc = res/gas_concentration
      case default
        res = 0.d0
        dres_dc = 0.d0
    end select

    rt_auxvar%total_sorb_eq_gas(irxn) = rt_auxvar%total_sorb_eq_gas(irxn) + res
    rt_auxvar%dtotal_sorb_eq_gas(irxn,irxn) = &
    rt_auxvar%dtotal_sorb_eq_gas(irxn,irxn) + dres_dc

  enddo

end subroutine RTotalSorbGasKD

! ************************************************************************** !

subroutine RAccumulationSorbGas(rt_auxvar,global_auxvar,material_auxvar, &
                             reaction,option,Res)
  ! 
  ! Computes sorbed portion of the accumulation term in
  ! residual function, from the gas phase
  ! 
  ! Author: Michael Nole
  ! Date: 06/23/21
  ! 

  use Option_module
  use Material_Aux_class

  implicit none

  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscReal :: Res(reaction%ncomp)
  PetscInt :: irxn, icomp

  ! units = (mol solute/m^3 bulk)*(m^3 bulk)/(sec) = mol/sec
  ! all residual entries should be in mol/sec
!  v_t = material_auxvar%volume/option%tran_dt
  do irxn = 1, reaction%gas%isotherm%neqkdrxn
    icomp = reaction%gas%isotherm%eqkdspecid(irxn)
    Res(icomp) = Res(icomp) + rt_auxvar%total_sorb_eq_gas(irxn)* &
                              material_auxvar%volume
  enddo

end subroutine RAccumulationSorbGas

! ************************************************************************** !

subroutine RTotalCO2(rt_auxvar,global_auxvar,reaction,option)
  ! 
  ! Computes the total component concentrations and derivative with
  ! respect to free-ion for CO2 modes; this is legacy cod3
  ! 
  ! Author: Glenn Hammond, but originally by Chuan Lu
  ! Date: 08/01/16
  ! 

  use Option_module
  use EOS_Water_module
  use co2eos_module, only: Henry_duan_sun

  implicit none
  
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscErrorCode :: ierr
  PetscInt :: iphase
  PetscInt :: icomp
  PetscReal :: tempreal
  PetscReal :: dg,dddt,dddp,fg,dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,&
               yco2,pco2,sat_pressure,lngamco2
  PetscReal :: den_kg_per_L
  PetscReal :: den
  PetscReal :: lnQK
  PetscReal :: m_cl, m_na, muco2, xmass
  PetscReal :: pressure, temperature, xphico2
  PetscInt :: iactgas

! *********** Add SC phase and gas contributions ***********************  
  ! CO2-specific
  iphase = 2

  rt_auxvar%total(:,iphase) = 0.D0
  rt_auxvar%aqueous%dtotal(:,:,iphase) = 0.D0

  if (associated(global_auxvar%xmass)) xmass = global_auxvar%xmass(iphase)
  den_kg_per_L = global_auxvar%den_kg(iphase)*xmass*1.d-3

  if (global_auxvar%sat(iphase) > 1.D-20) then
    do iactgas = 1, reaction%gas%nactive_gas ! all gas phase species are secondary

      pressure = global_auxvar%pres(2)
      temperature = global_auxvar%temp
      xphico2 = global_auxvar%fugacoeff(1)
      den = global_auxvar%den(2)
 
      call EOSWaterSaturationPressure(temperature, sat_pressure, ierr)
      pco2 = pressure - sat_pressure
!     call co2_span_wagner(pressure*1.D-6,temperature+273.15D0,dg,dddt,dddp,fg, &
!              dfgdp,dfgdt,eng,hg,dhdt,dhdp,visg,dvdt,dvdp,option%itable)
!
!            fg = fg*1D6
!            xphico2 = fg / pco2
!            global_auxvar%fugacoeff(1) = xphico2


      if (abs(reaction%species_idx%co2_gas_id) == iactgas ) then

        if (reaction%species_idx%na_ion_id /= 0 .and. reaction%species_idx%cl_ion_id /= 0) then
          m_na = rt_auxvar%pri_molal(reaction%species_idx%na_ion_id)
          m_cl = rt_auxvar%pri_molal(reaction%species_idx%cl_ion_id)
          call Henry_duan_sun(temperature,pressure*1D-5,muco2, &
                lngamco2,m_na,m_cl)
        else
          call Henry_duan_sun(temperature,pressure*1D-5,muco2, &
                lngamco2,option%m_nacl,option%m_nacl)
        endif
        !lnQk = - log(muco2) 
        lnQk = - log(muco2)-lngamco2
           
      else   
        lngamco2 = 0.d0
        lnQK = -reaction%gas%acteqlogK(iactgas)*LOG_TO_LN
      endif 
          
      if (reaction%gas%acteqh2oid(iactgas) > 0) then
        lnQK = lnQK + reaction%gas%acteqh2ostoich(iactgas)*rt_auxvar%ln_act_h2o
      endif
   
   ! contribute to %total          
   !     do i = 1, ncomp
   ! removed loop over species, suppose only one primary species is related
      icomp = reaction%gas%acteqspecid(1,iactgas)
      pressure = pressure * 1.D-5
        
!     rt_auxvar%gas_pp(iactgas) = &
!         exp(lnQK+lngamco2)*rt_auxvar%pri_molal(icomp) &
!         /(IDEAL_GAS_CONSTANT*1.d-2*(temperature+273.15D0)*xphico2)

!     This form includes factor Z in pV = ZRT for nonideal gas
      rt_auxvar%gas_pp(iactgas) = &
          exp(lnQK)*rt_auxvar%pri_act_coef(icomp)*rt_auxvar%pri_molal(icomp)* &
          den/pressure/xphico2

      rt_auxvar%total(icomp,iphase) = rt_auxvar%total(icomp,iphase) + &
          reaction%gas%acteqstoich(1,iactgas)* &
          rt_auxvar%gas_pp(iactgas)

!       print *,'RTotal: ',icomp,iactgas,pressure, temperature, xphico2, &
!         global_auxvar%sat(iphase),rt_auxvar%gas_pp(iactgas), &
!         rt_auxvar%pri_act_coef(icomp)*exp(lnQK)*rt_auxvar%pri_molal(icomp) &
!         /pressure/xphico2*den


   ! contribute to %dtotal
   !      tempreal = exp(lnQK+lngamco2)/pressure/xphico2*den
!     tempreal = rt_auxvar%pri_act_coef(icomp)*exp(lnQK) &
!         /pressure/xphico2*den
      tempreal = rt_auxvar%gas_pp(iactgas)/rt_auxvar%pri_molal(icomp)
      rt_auxvar%aqueous%dtotal(icomp,icomp,iphase) = &
          rt_auxvar%aqueous%dtotal(icomp,icomp,iphase) + &
          reaction%gas%acteqstoich(1,iactgas)*tempreal
    enddo
  ! rt_auxvar%total(:,iphase) = rt_auxvar%total(:,iphase)!*den_kg_per_L
  ! units of dtotal = kg water/L water
  ! rt_auxvar%dtotal(:, :,iphase) = rt_auxvar%dtotal(:,:,iphase)!*den_kg_per_L
  endif

end subroutine RTotalCO2
  
end module Reaction_Gas_module
