module Reaction_Microbial_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Microbial_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: MicrobialRead, &
            RMicrobial

contains

! ************************************************************************** !

subroutine MicrobialRead(microbial,input,option)
  !
  ! Reads chemical species
  !
  ! Author: Glenn Hammond
  ! Date: 08/16/12
  !
  ! Edited by: Peishi Jiang
  ! Date: 06/01/13
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  type(microbial_type) :: microbial
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  type(microbial_rxn_type), pointer :: microbial_rxn, cur_microbial_rxn
  type(monod_type), pointer :: monod, prev_monod
  type(microbial_biomass_type), pointer :: microbial_biomass
  type(inhibition_type), pointer :: inhibition, prev_inhibition

  microbial%nrxn = microbial%nrxn + 1

  microbial_rxn => MicrobialRxnCreate()
  nullify(prev_monod)
  nullify(prev_inhibition)
  nullify(microbial_biomass)
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,MICROBIAL_REACTION')
    call StringToUpper(word)

    select case(trim(word))
      case('REACTION')
        ! remainder of string should be the reaction equation
        microbial_rxn%reaction = trim(adjustl(input%buf))
        ! set flag for error message
        if (len_trim(microbial_rxn%reaction) < 2) input%ierr = 1
        call InputErrorMsg(input,option,'reaction string', &
                            'CHEMISTRY,MICROBIAL_REACTION,REACTION')
      case('RATE_CONSTANT')
        call InputReadDouble(input,option,microbial_rxn%rate_constant)
        call InputErrorMsg(input,option,'rate constant', &
                           'CHEMISTRY,MICROBIAL_REACTION')
      case('ACTIVATION_ENERGY')
        call InputReadDouble(input,option,microbial_rxn%activation_energy)
        call InputErrorMsg(input,option,'activation energy', &
                           'CHEMISTRY,MICROBIAL_REACTION')
        call InputReadAndConvertUnits(input,microbial_rxn%activation_energy, &
                     'J/mol', &
                     'CHEMISTRY,MICROBIAL_REACTION,ACTIVATION_ENERGY',option)
      case('MONOD')
        monod => MicrobialMonodCreate()
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,MICROBIAL_REACTION,MONOD')
          call StringToUpper(word)
          select case(trim(word))
            case('SPECIES_NAME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'species name', &
                           'CHEMISTRY,MICROBIAL_REACTION,MONOD')
              monod%species_name = word
            case('HALF_SATURATION_CONSTANT')
              call InputReadDouble(input,option,monod%half_saturation_constant)
              call InputErrorMsg(input,option,'half saturation constant', &
                                 'CHEMISTRY,MICROBIAL_REACTION,MONOD')
            case('THRESHOLD_CONCENTRATION')
              call InputReadDouble(input,option,monod%threshold_concentration)
              call InputErrorMsg(input,option,'threshold concdntration', &
                                 'CHEMISTRY,MICROBIAL_REACTION,MONOD')
            case default
              call InputKeywordUnrecognized(input,word, &
                                            'CHEMISTRY,MICROBIAL_REACTION,MONOD', &
                                            option)
          end select
        enddo
        call InputPopBlock(input,option)
        ! append to list
        if (.not.associated(microbial_rxn%monod)) then
          microbial_rxn%monod => monod
        else
          prev_monod%next => monod
        endif
        prev_monod => monod
        nullify(monod)
      case('INHIBITION')
        inhibition => MicrobialInhibitionCreate()
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,MICROBIAL_REACTION,INHIBITION')
          call StringToUpper(word)
          select case(trim(word))
            case('SPECIES_NAME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'species name', &
                                 'CHEMISTRY,MICROBIAL_REACTION,INHIBITION')
              inhibition%species_name = word
            case('TYPE')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'inhibition type', &
                                 'CHEMISTRY,MICROBIAL_REACTION,INHIBITION')
              call StringToUpper(word)
              select case(word)
                case('MONOD')
                  inhibition%itype = INHIBITION_MONOD
                case('INVERSE_MONOD')
                  inhibition%itype = INHIBITION_INVERSE_MONOD
                case('THRESHOLD')
                  inhibition%itype = INHIBITION_THRESHOLD
                  call InputReadDouble(input,option, &
                                       inhibition%inhibition_constant2)
                  call InputErrorMsg(input,option,'scaling factor', &
                                     'CHEMISTRY,MICROBIAL_REACTION,&
                                     &INHIBITION,THRESHOLD_INHIBITION')
                case('THRESHOLD_SMOOTHSTEP')
                  inhibition%itype = INHIBITION_THRESHOLD_SMOOTHSTEP
                  call InputReadDouble(input,option, &
                                       inhibition%inhibition_constant_smoothstep)
                  call InputErrorMsg(input,option,'the interval of the smoothstep function', &
                                     'CHEMISTRY,MICROBIAL_REACTION,&
                                     &INHIBITION,THRESHOLD_INHIBITION')
                case default
                  call InputKeywordUnrecognized(input,word, &
                         'CHEMISTRY,MICROBIAL_REACTION,INHIBITION,TYPE',option)
              end select
            case('INHIBITION_CONSTANT')
              call InputReadDouble(input,option,inhibition%inhibition_constant)
              call InputErrorMsg(input,option,'inhibition constant', &
                                 'CHEMISTRY,MICROBIAL_REACTION,INHIBITION')
            case default
              call InputKeywordUnrecognized(input,word, &
                      'CHEMISTRY,MICROBIAL_REACTION,INHIBITION',option)
          end select
        enddo
        call InputPopBlock(input,option)
        if (len_trim(inhibition%species_name) < 2 .or. &
            inhibition%itype == 0 .or. &
            Uninitialized(inhibition%inhibition_constant)) then
          option%io_buffer = 'A SPECIES_NAME, TYPE, and INHIBITION_CON' // &
            'STANT must be defined for INHIBITION in MICROBIAL_REACTION ' // &
            'with REACTION "' // trim(microbial_rxn%reaction) // '".'
          call PrintErrMsg(option)
        endif
        ! append to list
        if (.not.associated(microbial_rxn%inhibition)) then
          microbial_rxn%inhibition => inhibition
        else
          prev_inhibition%next => inhibition
        endif
        prev_inhibition => inhibition
        nullify(inhibition)
      case('BIOMASS')
        if (associated(microbial_biomass)) then
          call MicrobialBiomassDestroy(microbial_biomass)
        endif
        microbial_biomass => MicrobialBiomassCreate()
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword', &
                             'CHEMISTRY,MICROBIAL_REACTION,BIOMASS')
          call StringToUpper(word)
          select case(trim(word))
            case('SPECIES_NAME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'species name', &
                                 'CHEMISTRY,MICROBIAL_REACTION,BIOMASS')
              microbial_biomass%species_name = word
            case('YIELD')
              call InputReadDouble(input,option,microbial_biomass%yield)
              call InputErrorMsg(input,option,'yield', &
                                 'CHEMISTRY,MICROBIAL_REACTION,BIOMASS')
            case default
              call InputKeywordUnrecognized(input,word, &
                                      'CHEMISTRY,MICROBIAL_REACTION,BIOMASS', &
                                            option)
          end select
        enddo
        call InputPopBlock(input,option)
      case default
        call InputKeywordUnrecognized(input,word,'CHEMISTRY,MICROBIAL_REACTION', &
                                      option)
    end select
  enddo
  call InputPopBlock(input,option)

  ! add linkage to biomass if exists
  microbial_rxn%biomass => microbial_biomass

  ! add to list
  if (.not.associated(microbial%microbial_rxn_list)) then
    microbial%microbial_rxn_list => microbial_rxn
    microbial_rxn%id = 1
  else
    cur_microbial_rxn => microbial%microbial_rxn_list
    do
      if (.not.associated(cur_microbial_rxn%next)) then
        cur_microbial_rxn%next => microbial_rxn
        microbial_rxn%id = cur_microbial_rxn%id + 1
        exit
      endif
      cur_microbial_rxn => cur_microbial_rxn%next
    enddo
  endif
  nullify(microbial_rxn)

end subroutine MicrobialRead

! ************************************************************************** !

subroutine RMicrobial(Res,Jac,compute_derivative,rt_auxvar, &
                      global_auxvar,material_auxvar,reaction,option)
  !
  ! Computes the microbial reaction
  !
  ! Author: Glenn Hammond
  ! Date: 10/31/12
  !
  ! Edited by: Peishi Jiang
  ! Date: 06/01/13
  !

  use String_module
  use Option_module, only : option_type
  use Reactive_Transport_Aux_module, only : reactive_transport_auxvar_type
  use Global_Aux_module, only : global_auxvar_type
  use Material_Aux_class, only : material_auxvar_type
  use Reaction_Aux_module, only : reaction_rt_type
  use Reaction_Immobile_Aux_module, only : immobile_type

  implicit none

  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: por_sat_vol
  PetscInt :: irxn, i, j, ii, icomp, jcomp, ncomp
  PetscInt :: imonod, iinhibition, ibiomass
  PetscReal :: Im
  PetscReal :: rate_constant
  PetscReal :: activity
  PetscReal :: act_coef
  PetscReal :: monod(10)
  PetscReal :: inhibition(10)
  PetscReal :: biomass_conc, yield
  PetscInt :: immobile_id
  PetscReal :: denominator, dR_dX, dX_dc, dR_dc, dR_dbiomass
  PetscReal :: tempreal
  ! Below are some intermediate variables for calculating the smoothstep function for inhibition
  PetscReal :: log_inhibition, interval, log_activity, lower, upper, z
  type(microbial_type), pointer :: microbial
  type(immobile_type), pointer :: immobile

  microbial => reaction%microbial
  immobile => reaction%immobile

  ! units:
  ! Residual: mol/sec
  ! Jacobian: (mol/sec)*(kg water/mol) = kg water/sec

  if (reaction%logging_verbosity > 40) then
    do i = 1, reaction%ncomp
      print *, "Jacobian: ",i,(j,Jac(i,j),j=1,reaction%ncomp)
    enddo
      print *, "Residual: ",(Res(j),j=1,reaction%ncomp)
    print *, " "
  endif

  do irxn = 1, microbial%nrxn

    ! units:
    !   without biomass: mol/L-sec
    !   with biomass: mol/L-sec * (m^3 bulk / mol biomass)

    ncomp = microbial%specid(0,irxn)
    rate_constant = microbial%rate_constant(irxn)
    Im = rate_constant
    if (associated(microbial%activation_energy)) then
      ! ideal gas constant units: J/mol-K
      Im = Im * exp(microbial%activation_energy(irxn)/IDEAL_GAS_CONSTANT* &
                    (1.d0/298.15d0-1.d0/(global_auxvar%temp+273.15d0)))
    endif
    yield = 0.d0
    biomass_conc = 0.d0

    ! monod expressions
    do ii = 1, microbial%monodid(0,irxn)
      imonod = microbial%monodid(ii,irxn)
      icomp = microbial%monod_specid(imonod)
      activity = rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp)
      monod(ii) = (activity - microbial%monod_Cth(imonod)) / &
                  (microbial%monod_K(imonod) + activity - &
                   microbial%monod_Cth(imonod))
      Im = Im*monod(ii)
    enddo

    ! inhibition expressions
    do ii = 1, microbial%inhibitionid(0,irxn)
      iinhibition = microbial%inhibitionid(ii,irxn)
      icomp = microbial%inhibition_specid(iinhibition)
      activity = rt_auxvar%pri_molal(icomp)*rt_auxvar%pri_act_coef(icomp)
      select case(microbial%inhibition_type(iinhibition))
        case(INHIBITION_MONOD)
          inhibition(ii) = microbial%inhibition_C(iinhibition) / &
                          (microbial%inhibition_C(iinhibition) + activity)
        case(INHIBITION_INVERSE_MONOD)
          inhibition(ii) = activity / &
                          (microbial%inhibition_C(iinhibition) + activity)
        case(INHIBITION_THRESHOLD)
          ! if microbial%inhibition_C2 is negative, inhibition kicks in
          ! above the concentration
          inhibition(ii) = 0.5d0 + &
                           atan((activity - &
                                 microbial%inhibition_C(iinhibition)) * &
                                microbial%inhibition_C2(iinhibition)) / PI
        case(INHIBITION_THRESHOLD_SMOOTHSTEP)
          ! TODO: use the smoothstep function for inhibition
          log_inhibition = LOG10(microbial%inhibition_C(iinhibition))
          interval = microbial%inhibition_C_smoothstep(iinhibition)
          log_activity = LOG10(activity)
          lower = log_inhibition - interval / 2.d0
          upper = log_inhibition + interval / 2.d0
          z = (log_activity - lower) / interval
          if (z < 0.) then
            inhibition(ii) = 0.d0
          else if (z > 1.) then
            inhibition(ii) = 1.d0
          else
            inhibition(ii) = 3.d0 * z ** 2 - 2.d0 * z ** 3
          endif

          if (reaction%logging_verbosity > 40) then
            print *, ' Get the inhibition value ' // &
              trim(StringWrite(inhibition(ii))) // ' with z ' // &
              trim(StringWrite(z)) // ' for concentration' // &
              trim(StringWrite(activity))
          endif
      end select
      Im = Im*inhibition(ii)
    enddo

    ! biomass term
    ibiomass = microbial%biomassid(irxn)
    if (ibiomass > 0) then
      immobile_id = reaction%offset_immobile + ibiomass
      biomass_conc = rt_auxvar%immobile(ibiomass)
      yield = microbial%biomass_yield(irxn)
      Im = Im*biomass_conc
    endif

    ! por_sat_vol units: m^3 water
    por_sat_vol = material_auxvar%porosity*global_auxvar%sat(iphase)* &
                  material_auxvar%volume

    ! Im units (before): mol/L-sec
    Im = Im * 1.d3*por_sat_vol
    ! Im units (after): mol/sec

    do i = 1, ncomp
      icomp = microbial%specid(i,irxn)
      Res(icomp) = Res(icomp) - microbial%stoich(i,irxn)*Im
    enddo

    if (ibiomass > 0) then
      Res(immobile_id) = Res(immobile_id) - yield*Im
    endif

    if (.not. compute_derivative) cycle

    ! monod expressions
    do ii = 1, microbial%monodid(0,irxn)
      imonod = microbial%monodid(ii,irxn)
      jcomp = microbial%monod_specid(imonod)
      act_coef = rt_auxvar%pri_act_coef(jcomp)
      activity = rt_auxvar%pri_molal(jcomp)*act_coef

      dR_dX = Im / monod(ii)

      denominator = microbial%monod_K(imonod) + activity - &
                    microbial%monod_Cth(imonod)

      dX_dc = act_coef / denominator - &
              act_coef * (activity - microbial%monod_Cth(imonod)) / &
              (denominator*denominator)

      dR_dc = -1.d0*dR_dX*dX_dc
      do i = 1, ncomp
        icomp = microbial%specid(i,irxn)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                            microbial%stoich(i,irxn)*dR_dc
      enddo
      if (ibiomass > 0) then
        Jac(immobile_id,jcomp) = Jac(immobile_id,jcomp) + yield*dR_dc
      endif
    enddo

    ! inhibition expressions
    do ii = 1, microbial%inhibitionid(0,irxn)
      iinhibition = microbial%inhibitionid(ii,irxn)
      jcomp = microbial%inhibition_specid(iinhibition)
      act_coef = rt_auxvar%pri_act_coef(jcomp)
      activity = rt_auxvar%pri_molal(jcomp)*act_coef

      if (inhibition(ii) == 0) then
        ! This is used to avoid NaN when inhibition kicks in.
        ! Note that when the inhibition is zero, we expect its gradient is also zero
        ! So the following works.
        dR_dX = 0.d0
      else
        dR_dX = Im / inhibition(ii)
      endif
      

      select case(microbial%inhibition_type(iinhibition))
        case(INHIBITION_MONOD)
          denominator = microbial%inhibition_C(iinhibition) + activity
          dX_dc = -1.d0 * act_coef *microbial%inhibition_C(iinhibition) / &
                  (denominator*denominator)
        case(INHIBITION_INVERSE_MONOD)
          denominator = microbial%inhibition_C(iinhibition) + activity
          dX_dc = act_coef / denominator - &
                  act_coef * activity / &
                  (denominator*denominator)
        case(INHIBITION_THRESHOLD)
          ! derivative of atan(X) = 1 / (1 + X^2) dX
          tempreal = (activity - microbial%inhibition_C(iinhibition)) * &
                     microbial%inhibition_C2(iinhibition)
          dX_dc = (microbial%inhibition_C2(iinhibition) * act_coef / &
                   (1.d0 + tempreal*tempreal)) / PI
        case(INHIBITION_THRESHOLD_SMOOTHSTEP)
          ! TODO 
          log_inhibition = LOG10(microbial%inhibition_C(iinhibition))
          interval = microbial%inhibition_C_smoothstep(iinhibition)
          log_activity = LOG10(activity)
          lower = log_inhibition - interval / 2.d0 
          upper = log_inhibition + interval / 2.d0
          z = (log_activity - lower) / interval
          if (z < 0.) then
            dX_dc = 0.d0
          else if (z > 1.d0) then
            dX_dc = 0.d0
          else
            dX_dc = (6.d0*z - 6.d0*z**2) / (interval*activity*LOG(10.d0))
          endif

          if (reaction%logging_verbosity > 40) then
            print *, ' Get the gradient of the inhibition ' // &
              trim(StringWrite(dX_dc)) // ' with z ' // &
              trim(StringWrite(z)) // ' for concentration' // &
              trim(StringWrite(activity))
          endif
      end select

      dR_dc = -1.d0*dR_dX*dX_dc
      do i = 1, ncomp
        icomp = microbial%specid(i,irxn)
        ! units = (mol/sec)*(kg water/mol) = kg water/sec
        Jac(icomp,jcomp) = Jac(icomp,jcomp) + &
                            microbial%stoich(i,irxn)*dR_dc
      enddo
      if (ibiomass > 0) then
        Jac(immobile_id,jcomp) = Jac(immobile_id,jcomp) + yield*dR_dc
      endif
    enddo

    ! biomass expression
    if (ibiomass > 0) then
      dR_dbiomass = -1.d0*Im / biomass_conc
      do i = 1, ncomp
        icomp = microbial%specid(i,irxn)
        ! units = (mol/sec)*(m^3/mol) = m^3/sec
        Jac(icomp,immobile_id) = Jac(icomp,immobile_id) + &
                            microbial%stoich(i,irxn)*dR_dbiomass
      enddo
      Jac(immobile_id,immobile_id) = Jac(immobile_id,immobile_id) + &
        yield*dR_dbiomass
    endif

    ! if (reaction%logging_verbosity > 40) then
    !   do i = 1, reaction%ncomp
    !     print *, "Jacobian: ",i,(j,Jac(i,j),j=1,reaction%ncomp)
    !   enddo
    !   print *, " "
    ! endif

  enddo

  if (reaction%logging_verbosity > 40) then
    do i = 1, reaction%ncomp
      print *, "Jacobian: ",i,(j,Jac(i,j),j=1,reaction%ncomp)
    enddo
      print *, "Residual: ",(Res(j),j=1,reaction%ncomp)
    print *, " "
  endif

end subroutine RMicrobial

end module Reaction_Microbial_module
