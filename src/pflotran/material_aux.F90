module Material_Aux_class
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private


  PetscInt, parameter, public :: perm_xx_index = 1
  PetscInt, parameter, public :: perm_yy_index = 2
  PetscInt, parameter, public :: perm_zz_index = 3
  PetscInt, parameter, public :: perm_xy_index = 4
  PetscInt, parameter, public :: perm_yz_index = 5
  PetscInt, parameter, public :: perm_xz_index = 6
  
  ! do not use 0 as an index as there is a case statement in material.F90
  ! designed to catch erroneous values outside [1,2].
  PetscInt, parameter, public :: POROSITY_CURRENT = 1
  PetscInt, parameter, public :: POROSITY_BASE = 2
  PetscInt, parameter, public :: POROSITY_INITIAL = 3

  ! Tensor to scalar conversion models
  ! default for structured grids = TENSOR_TO_SCALAR_LINEAR
  ! default for unstructured grids = TENSOR_TO_SCALAR_POTENTIAL
  ! Both are set in discretization.F90:DiscretizationReadRequiredCards() 
  ! immediately after the GRID cards is read with a call to 
  ! MaterialAuxSetPermTensorModel()
  PetscInt, parameter, public :: TENSOR_TO_SCALAR_LINEAR = 1
  PetscInt, parameter, public :: TENSOR_TO_SCALAR_FLOW = 2
  PetscInt, parameter, public :: TENSOR_TO_SCALAR_POTENTIAL = 3

  ! flag to determine which model to use for tensor to scalar conversion 
  ! of permeability
  PetscInt :: perm_tens_to_scal_model = TENSOR_TO_SCALAR_LINEAR
  
!  PetscInt, public :: soil_thermal_conductivity_index
!  PetscInt, public :: soil_heat_capacity_index
  PetscInt, public :: soil_compressibility_index
  PetscInt, public :: soil_reference_pressure_index
  PetscInt, public :: max_material_index
  
  type, public :: material_auxvar_type
    PetscInt :: id
    PetscReal :: volume
    PetscReal :: porosity_0 ! initial porosity as defined in input file or 
                            ! initial conditon
    PetscReal :: porosity_base ! base porosity prescribed by pm outside flow 
                               ! (e.g. geomechanics, mineral precip/diss)
    PetscReal :: porosity ! porosity used in calculation, which may be a 
                          ! function of soil compressibity, etc.
    PetscReal :: dporosity_dp
    PetscReal :: tortuosity
    PetscReal :: soil_particle_density
    PetscReal, pointer :: permeability(:)
    PetscReal, pointer :: sat_func_prop(:)
    PetscReal, pointer :: soil_properties(:) ! den, therm. cond., heat cap.
    type(fracture_auxvar_type), pointer :: fracture
    PetscReal, pointer :: geomechanics_subsurface_prop(:)
    PetscInt :: creep_closure_id

!    procedure(SaturationFunction), nopass, pointer :: SaturationFunction
  contains
    procedure, public :: PermeabilityTensorToScalar => &
                           MaterialDiagPermTensorToScalar
    procedure, public :: PermeabilityTensorToScalarSafe => &
                           MaterialDiagPermTensorToScalarSafe
  end type material_auxvar_type
  
  type, public :: fracture_auxvar_type
    PetscBool :: fracture_is_on
    PetscReal :: initial_pressure
    PetscReal :: properties(4)
    PetscReal :: vector(3) ! < 0. 0. 0. >
    PetscInt :: id
    PetscBool :: unit_test
    character(len=MAXSTRINGLENGTH) :: input_filename
  end type fracture_auxvar_type
 
  type, public :: material_parameter_type
    PetscReal, pointer :: soil_heat_capacity(:) ! MJ/kg rock-K
    PetscReal, pointer :: soil_thermal_conductivity(:,:) ! W/m-K
  end type material_parameter_type  
  
  type, public :: material_type
    PetscReal :: time_t, time_tpdt  
    PetscInt :: num_aux
    type(material_parameter_type), pointer :: material_parameter
    class(material_auxvar_type), pointer :: auxvars(:)
  end type material_type
  
  ! procedure pointer declarations
  procedure(MaterialCompressSoilDummy), pointer :: &
    MaterialCompressSoilPtr => null()
 
  ! interface blocks
  interface
    subroutine MaterialCompressSoilDummy(auxvar,pressure,compressed_porosity, &
                                         dcompressed_porosity_dp)
    import material_auxvar_type
    implicit none
    class(material_auxvar_type), intent(in) :: auxvar
    PetscReal, intent(in) :: pressure
    PetscReal, intent(out) :: compressed_porosity
    PetscReal, intent(out) :: dcompressed_porosity_dp
    end subroutine MaterialCompressSoilDummy
  end interface 
  
  interface MaterialCompressSoil
    procedure MaterialCompressSoilPtr
  end interface
  
  public :: MaterialCompressSoilDummy, &
            MaterialCompressSoilPtr, &
            MaterialCompressSoil, &
            MaterialCompressSoilBRAGFLO, &
            MaterialCompressSoilPoroExp, &
            MaterialCompressSoilLeijnse, &
            MaterialCompressSoilLinear, &
            MaterialCompressSoilQuadratic, &
            MaterialCompressSoilTest, &
            MaterialCompressSoilUnitTest
  
  public :: MaterialAuxCreate, &
            MaterialAuxVarInit, &
            MaterialAuxVarCopy, &
            MaterialAuxVarStrip, &
            MaterialAuxVarGetValue, &
            MaterialAuxVarSetValue, &
            MaterialAuxIndexToPropertyName, &
            MaterialAuxDestroy, &
            MaterialAuxVarFractureStrip, &
            MaterialAuxSetPermTensorModel

  public :: MaterialAuxVarCompute
  
contains

! ************************************************************************** !

function MaterialAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  type(material_type), pointer :: MaterialAuxCreate
  
  type(material_type), pointer :: aux

  allocate(aux)
  nullify(aux%auxvars)
  allocate(aux%material_parameter)
  nullify(aux%material_parameter%soil_heat_capacity)
  nullify(aux%material_parameter%soil_thermal_conductivity)
  aux%num_aux = 0
  aux%time_t = 0.d0
  aux%time_tpdt = 0.d0

  MaterialAuxCreate => aux
  
end function MaterialAuxCreate

! ************************************************************************** !

subroutine MaterialAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  class(material_auxvar_type) :: auxvar
  type(option_type) :: option
  
  auxvar%id = UNINITIALIZED_INTEGER
  auxvar%volume = UNINITIALIZED_DOUBLE
  auxvar%porosity_0 = UNINITIALIZED_DOUBLE
  auxvar%porosity_base = UNINITIALIZED_DOUBLE
  auxvar%porosity = UNINITIALIZED_DOUBLE
  auxvar%dporosity_dp = 0.d0
  auxvar%tortuosity = UNINITIALIZED_DOUBLE
  auxvar%soil_particle_density = UNINITIALIZED_DOUBLE
  if (option%iflowmode /= NULL_MODE) then
    allocate(auxvar%permeability(3))
    auxvar%permeability = UNINITIALIZED_DOUBLE
  else
    nullify(auxvar%permeability)
  endif
  nullify(auxvar%sat_func_prop)
  nullify(auxvar%fracture)
  auxvar%creep_closure_id = 1
  
  if (max_material_index > 0) then
    allocate(auxvar%soil_properties(max_material_index))
    ! initialize these to zero for now
    auxvar%soil_properties = UNINITIALIZED_DOUBLE
  else
    nullify(auxvar%soil_properties)
  endif

  nullify(auxvar%geomechanics_subsurface_prop)
  
end subroutine MaterialAuxVarInit

! ************************************************************************** !

subroutine MaterialAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  class(material_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option
  
  auxvar2%volume = auxvar%volume
  auxvar2%porosity_0 = auxvar%porosity_0
  auxvar2%porosity_base = auxvar%porosity_base
  auxvar2%porosity = auxvar%porosity
  auxvar2%tortuosity = auxvar%tortuosity
  auxvar2%soil_particle_density = auxvar%soil_particle_density
  if (associated(auxvar%permeability)) then
    auxvar2%permeability = auxvar%permeability
  endif
  if (associated(auxvar%sat_func_prop)) then
    auxvar2%sat_func_prop = auxvar%sat_func_prop
  endif
  if (associated(auxvar%soil_properties)) then
    auxvar2%soil_properties = auxvar%soil_properties
  endif
  auxvar2%creep_closure_id = auxvar%creep_closure_id
  
end subroutine MaterialAuxVarCopy

! ************************************************************************** !

subroutine MaterialAuxSetPermTensorModel(model,option)

  use Option_module

  implicit none

  PetscInt :: model
  type(option_type) :: option

  !! simple if longwinded little safety measure here, the calling routine
  !! should also check that model is a sane number but in principle this
  !! routine should protect itself too. 
  !! Please note that if you add a new model type above then you MUST add
  !! it to this little list here too. 
  if (model == TENSOR_TO_SCALAR_LINEAR .OR. &
      model == TENSOR_TO_SCALAR_FLOW .OR. &
      model == TENSOR_TO_SCALAR_POTENTIAL) then
    perm_tens_to_scal_model = model
  else
    option%io_buffer  = 'MaterialDiagPermTensorToScalar: tensor to scalar &
                         &model type is not recognized.'
    call PrintErrMsg(option)
  endif

end subroutine MaterialAuxSetPermTensorModel

! ************************************************************************** !

subroutine MaterialDiagPermTensorToScalar(material_auxvar,dist, &
                                      scalar_permeability)
  ! 
  ! Transforms a diagonal permeability tensor to a scalar through a dot 
  ! product.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 
  use Utility_module, only : Equal

  implicit none
  
  class(material_auxvar_type) :: material_auxvar
  ! -1 = fraction upwind
  ! 0 = magnitude
  ! 1 = unit x-dir
  ! 2 = unit y-dir
  ! 3 = unit z-dir
  PetscReal, intent(in) :: dist(-1:3)
  PetscReal, intent(out) :: scalar_permeability

  PetscReal :: kx, ky, kz

  kx = material_auxvar%permeability(perm_xx_index)
  ky = material_auxvar%permeability(perm_yy_index)
  kz = material_auxvar%permeability(perm_zz_index)
#if 0
  if (Equal(kx,ky) .and. Equal(ky,kz)) then
    scalar_permeability = kx
    return
  endif
#endif

  select case(perm_tens_to_scal_model)
    case(TENSOR_TO_SCALAR_LINEAR)
      scalar_permeability = DiagPermTensorToScalar_Linear(kx,ky,kz,dist)
    case(TENSOR_TO_SCALAR_FLOW)
      scalar_permeability = DiagPermTensorToScalar_Flow(kx,ky,kz,dist)
    case(TENSOR_TO_SCALAR_POTENTIAL)
      scalar_permeability = DiagPermTensortoScalar_Potential(kx,ky,kz,dist)
    case default
      ! as default, just do linear 
      !scalar_permeability = DiagPermTensorToScalar_Linear(kx,ky,kz,dist)
      ! as default, do perm in direction of flow
      !scalar_permeability = DiagPermTensorToScalar_Flow(kx,ky,kz,dist)
      ! as default, do perm in direction of potential gradient
      scalar_permeability = DiagPermTensorToScalar_Potential(kx,ky,kz,dist)
  end select


end subroutine MaterialDiagPermTensorToScalar

! ************************************************************************** !

subroutine MaterialDiagPermTensorToScalarSafe(material_auxvar,dist, &
                                      scalar_permeability)
  !
  ! Transforms a diagonal perm. tensor to a scalar through a dot product.
  ! This version will not generate NaNs for zero permeabilities
  !
  ! Author: Dave Ponting
  ! Date: 03/19/19
  !

  implicit none

  class(material_auxvar_type) :: material_auxvar

  PetscReal, intent(in) :: dist(-1:3)
  PetscReal, intent(out) :: scalar_permeability

  PetscReal :: kx, ky, kz

  kx = material_auxvar%permeability(perm_xx_index)
  ky = material_auxvar%permeability(perm_yy_index)
  kz = material_auxvar%permeability(perm_zz_index)

  select case(perm_tens_to_scal_model)
    case(TENSOR_TO_SCALAR_LINEAR)
      scalar_permeability = DiagPermTensorToScalar_Linear(kx,ky,kz,dist)
    case(TENSOR_TO_SCALAR_FLOW)
      scalar_permeability = DiagPermTensorToScalar_Flow(kx,ky,kz,dist)
    case(TENSOR_TO_SCALAR_POTENTIAL)
      scalar_permeability = DiagPermTensortoScalar_PotentialSafe(kx,ky,kz,dist)
    case default
      scalar_permeability = DiagPermTensorToScalar_PotentialSafe(kx,ky,kz,dist)
  end select

end subroutine MaterialDiagPermTensorToScalarSafe
! ************************************************************************** !

function DiagPermTensorToScalar_Linear(kx,ky,kz,dist)
  implicit none
  PetscReal :: DiagPermTensorToScalar_Linear
  PetscReal, intent(in) :: dist(-1:3)
  PetscReal :: kx,ky,kz

  DiagPermTensorToScalar_Linear = kx*dabs(dist(1))+ky*dabs(dist(2))+&
                                  kz*dabs(dist(3))

end function DiagPermTensorToScalar_Linear

! ************************************************************************** !

function DiagPermTensorToScalar_Flow(kx,ky,kz,dist)
  
  !Permeability in the direction of flow

  implicit none
  PetscReal :: DiagPermTensorToScalar_Flow
  PetscReal, intent(in) :: dist(-1:3)
  PetscReal :: kx,ky,kz

  DiagPermTensorToScalar_Flow = kx*dabs(dist(1))**2.0 + &
                                     ky*dabs(dist(2))**2.0 + &
                                     kz*dabs(dist(3))**2.0

end function DiagPermTensorToScalar_Flow

! ************************************************************************** !

function DiagPermTensorToScalar_Potential(kx,ky,kz,dist)
  
  !Permeability in the direction of the potential gradient

  implicit none
  PetscReal :: DiagPermTensorToScalar_Potential
  PetscReal, intent(in) :: dist(-1:3)
  PetscReal :: kx,ky,kz

  DiagPermTensorToScalar_Potential = 1.d0/(dist(1)*dist(1)/kx + &
                                         dist(2)*dist(2)/ky + &
                                         dist(3)*dist(3)/kz)

end function DiagPermTensorToScalar_Potential

! ************************************************************************** !

function DiagPermTensorToScalar_PotentialSafe(kx,ky,kz,dist)

  ! Permeability in the direction of the potential gradient
  ! This version will not generate NaNs for zero permeabilities
  !
  ! Author: Dave Ponting
  ! Date: 03/19/19
  !

  implicit none
  PetscReal :: DiagPermTensorToScalar_PotentialSafe
  PetscReal, intent(in) :: dist(-1:3)
  PetscReal :: kx, ky, kz, kxi, kyi, kzi, den, deni

  !  Form safe inverse permeabilities

  kxi = 0.0
  kyi = 0.0
  kzi = 0.0

  if (kx>0.0) kxi = 1.0/kx
  if (ky>0.0) kyi = 1.0/ky
  if (kz>0.0) kzi = 1.0/kz

  !  Form denominator

  den = dist(1)*dist(1)*kxi + &
        dist(2)*dist(2)*kyi + &
        dist(3)*dist(3)*kzi

  !  Form safe inverse denominator

  deni = 0.0
  if (den>0.0) deni=1.0/den

  !  Store final value

  DiagPermTensorToScalar_PotentialSafe = deni

end function DiagPermTensorToScalar_PotentialSafe

! ************************************************************************** !

function MaterialAuxVarGetValue(material_auxvar,ivar)
  ! 
  ! Returns the value of an entry in material_auxvar_type based on ivar.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/28/14
  ! 

  use Variables_module
  
  implicit none

  class(material_auxvar_type) :: material_auxvar 
  PetscInt :: ivar

  PetscReal :: MaterialAuxVarGetValue

  MaterialAuxVarGetValue = UNINITIALIZED_DOUBLE
  select case(ivar)
    case(VOLUME)
      MaterialAuxVarGetValue = material_auxvar%volume
    case(INITIAL_POROSITY)
      MaterialAuxVarGetValue = material_auxvar%porosity_0
    case(BASE_POROSITY)
      MaterialAuxVarGetValue = material_auxvar%porosity_base
    case(POROSITY)
      MaterialAuxVarGetValue = material_auxvar%porosity
    case(TORTUOSITY)
      MaterialAuxVarGetValue = material_auxvar%tortuosity
    case(PERMEABILITY_X)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_xx_index)
    case(PERMEABILITY_Y)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_yy_index)
    case(PERMEABILITY_Z)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_zz_index)
    case(PERMEABILITY_XY)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_xy_index)
    case(PERMEABILITY_YZ)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_yz_index)
    case(PERMEABILITY_XZ)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_xz_index)
    case(SOIL_COMPRESSIBILITY)
      MaterialAuxVarGetValue = material_auxvar% &
                                 soil_properties(soil_compressibility_index)
    case(SOIL_REFERENCE_PRESSURE)
      MaterialAuxVarGetValue = material_auxvar% &
                                 soil_properties(soil_reference_pressure_index)
  end select
  
end function MaterialAuxVarGetValue

! ************************************************************************** !

subroutine MaterialAuxVarSetValue(material_auxvar,ivar,value)
  ! 
  ! Sets the value of an entry in material_auxvar_type based on ivar.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/28/14
  ! 

  use Variables_module
  
  implicit none

  class(material_auxvar_type) :: material_auxvar 
  PetscInt :: ivar
  PetscReal :: value

  select case(ivar)
    case(VOLUME)
      material_auxvar%volume = value
    case(INITIAL_POROSITY)
      material_auxvar%porosity_0 = value
    case(BASE_POROSITY)
      material_auxvar%porosity_base = value
    case(POROSITY)
      material_auxvar%porosity = value
    case(TORTUOSITY)
      material_auxvar%tortuosity = value
    case(PERMEABILITY_X)
      material_auxvar%permeability(perm_xx_index) = value
    case(PERMEABILITY_Y)
      material_auxvar%permeability(perm_yy_index) = value
    case(PERMEABILITY_Z)
      material_auxvar%permeability(perm_zz_index) = value
    case(PERMEABILITY_XY)
      material_auxvar%permeability(perm_xy_index) = value
    case(PERMEABILITY_YZ)
      material_auxvar%permeability(perm_yz_index) = value
    case(PERMEABILITY_XZ)
      material_auxvar%permeability(perm_xz_index) = value
    case(SOIL_COMPRESSIBILITY)
      material_auxvar%soil_properties(soil_compressibility_index) = value
    case(SOIL_REFERENCE_PRESSURE)
      material_auxvar%soil_properties(soil_reference_pressure_index) = value
  end select
  
end subroutine MaterialAuxVarSetValue

! ************************************************************************** !

subroutine MaterialAuxVarCompute(auxvar,pressure)
  ! 
  ! Updates secondary material properties that are a function of state 
  ! variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/21/19
  ! 

  implicit none

  class(material_auxvar_type), intent(inout) :: auxvar
  PetscReal, intent(in) :: pressure
  
  auxvar%porosity = auxvar%porosity_base
  auxvar%dporosity_dp = 0.d0
  if (soil_compressibility_index > 0) then
    call MaterialCompressSoil(auxvar,pressure,auxvar%porosity, &
                              auxvar%dporosity_dp)
  endif
  
end subroutine MaterialAuxVarCompute

! ************************************************************************** !

subroutine MaterialCompressSoilLeijnse(auxvar,pressure, &
                                       compressed_porosity, &
                                       dcompressed_porosity_dp)
  ! 
  ! Calculates soil matrix compression based on Leijnse, 1992.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/14
  ! 

  implicit none

  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: compressibility
  PetscReal :: compression
  PetscReal :: tempreal
  
  compressibility = auxvar%soil_properties(soil_compressibility_index)
  compression = &
    exp(-1.d0 * compressibility * &
        (pressure - auxvar%soil_properties(soil_reference_pressure_index)))
  tempreal = (1.d0 - auxvar%porosity_base) * compression
  compressed_porosity = 1.d0 - tempreal
  dcompressed_porosity_dp = tempreal * compressibility
  
end subroutine MaterialCompressSoilLeijnse

! ************************************************************************** !

subroutine MaterialCompressSoilBRAGFLO(auxvar,pressure, &
                                       compressed_porosity, &
                                       dcompressed_porosity_dp)
  ! 
  ! Calculates soil matrix compression based on Eq. 9.6.9 of BRAGFLO
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/14
  ! 

  implicit none

  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: compressibility
  

  ! convert to pore compressiblity by dividing by base porosity
  compressibility = auxvar%soil_properties(soil_compressibility_index) / &
                    auxvar%porosity_base
  compressed_porosity = auxvar%porosity_base * &
    exp(compressibility * &
        (pressure - auxvar%soil_properties(soil_reference_pressure_index)))
  dcompressed_porosity_dp = compressibility * compressed_porosity
  
end subroutine MaterialCompressSoilBRAGFLO

! ************************************************************************** !

subroutine MaterialCompressSoilLinear(auxvar,pressure, &
                                      compressed_porosity, &
                                      dcompressed_porosity_dp)
  ! 
  ! Calculates soil matrix compression for standard constant 
  ! aquifer compressibility
  !
  ! variable 'alpha' is Freeze and Cherry, 1982
  ! 
  ! Author: Danny Birdsell and Satish Karra
  ! Date: 07/26/2016
  ! 

  implicit none

  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: compressibility
  
  compressibility = auxvar%soil_properties(soil_compressibility_index)
  compressed_porosity = auxvar%porosity_base + compressibility * &
            (pressure - auxvar%soil_properties(soil_reference_pressure_index)) 
  dcompressed_porosity_dp = compressibility

end subroutine MaterialCompressSoilLinear

! ************************************************************************** !

subroutine MaterialCompressSoilPoroExp(auxvar,pressure, &
                                       compressed_porosity, &
                                       dcompressed_porosity_dp)
  ! 
  ! Calculates soil matrix compression based on Eq. 9.6.9 of BRAGFLO
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/14
  ! 

  implicit none

  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: compressibility
  
  compressibility = auxvar%soil_properties(soil_compressibility_index)
  compressed_porosity = auxvar%porosity_base * &
    exp(compressibility * &
        (pressure - auxvar%soil_properties(soil_reference_pressure_index)))
  dcompressed_porosity_dp = compressibility * compressed_porosity
  
end subroutine MaterialCompressSoilPoroExp

! ************************************************************************** !

subroutine MaterialCompressSoilQuadratic(auxvar,pressure, &
                                         compressed_porosity, &
                                         dcompressed_porosity_dp)
  ! 
  ! Calculates soil matrix compression based on a quadratic model
  ! This is thedefaul model adopted in ECLIPSE 
  !
  ! Author: Paolo Orsini
  ! Date: 02/27/17
  ! 

  implicit none

  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: compressibility
  PetscReal :: compress_factor

  compressibility = auxvar%soil_properties(soil_compressibility_index)

  compress_factor = compressibility * &
          (pressure - auxvar%soil_properties(soil_reference_pressure_index))

  compressed_porosity = auxvar%porosity_base * &
          ( 1.0 + compress_factor + (compress_factor**2)/2.0 )
  
  dcompressed_porosity_dp = auxvar%porosity_base * &
          ( 1.0 + compress_factor) * compressibility  
  
end subroutine MaterialCompressSoilQuadratic

! ************************************************************************** !

subroutine MaterialCompressSoilTest(liq_pressure,auxvar)
  ! 
  ! Author: Jennifer Frederick
  ! Date: 10/22/2020
  ! 

  implicit none
  
  PetscReal, intent(in) :: liq_pressure
  class(material_auxvar_type), intent(in) :: auxvar

  PetscReal :: compressed_porosity, dcompressed_porosity_dp

  call MaterialCompressSoilBRAGFLO(auxvar,liq_pressure,compressed_porosity, &
                                   dcompressed_porosity_dp)
  print *, 'BRAGFLO :: [material ID], [liq. pressure], [intact porosity], &
           &[altered porosity], [d_por/d_press]'
  print *, auxvar%id, liq_pressure, auxvar%porosity_base, compressed_porosity, &
           dcompressed_porosity_dp

  call MaterialCompressSoilPoroExp(auxvar,liq_pressure,compressed_porosity, &
                                   dcompressed_porosity_dp)
  print *, 'POROEXP :: [material ID], [liq. pressure], [intact porosity], &
           &[altered porosity], [d_por/d_press]'
  print *, auxvar%id, liq_pressure, auxvar%porosity_base, compressed_porosity, &
           dcompressed_porosity_dp
  
end subroutine MaterialCompressSoilTest

! ************************************************************************** !

subroutine MaterialCompressSoilUnitTest(auxvars,grid)
  ! 
  ! Author: Jennifer Frederick
  ! Date: 10/22/2020
  ! 

  use Grid_module

  implicit none
  
  class(material_auxvar_type), pointer :: auxvars(:)
  type(grid_type), pointer :: grid

  class(material_auxvar_type), pointer :: auxvar
  PetscReal, pointer :: liq_pressure(:)                ! [Pa]
  PetscReal, pointer :: temp_liq_pressure(:)           ! [Pa]
  PetscReal, pointer :: material_id(:)                 ! [  ]    
  PetscReal, pointer :: temp_material_id(:)            ! [  ] 
  PetscReal, pointer :: corr_compressed_porosity(:)    ! [  ]
  PetscReal, pointer :: corr_dcompressed_porosity_dp(:)! [1/Pa]
  PetscReal, pointer :: temp_corr_comp_porosity(:)     ! [  ]
  PetscReal, pointer :: temp_corr_dcomp_porosity_dp(:) ! [1/Pa]
  PetscReal :: compressed_porosity                     ! [  ]
  PetscReal :: dcompressed_porosity_dp                 ! [1/Pa]        
  PetscReal, parameter :: tolerance = 1.d-8
  PetscReal :: diff
  PetscInt :: i, j, k
  PetscInt :: local_id, ghosted_id
  PetscInt :: prev_mat_id, mat_id
  character(len=MAXWORDLENGTH) :: pass_fail
  character(len=MAXWORDLENGTH) :: filename_in, filename_out, id
  PetscInt :: rc_in, rc_out, fu_in, fu_out
  PetscBool :: input_file_given
  PetscInt :: values(8)
  character(len=8) :: date
  character(len=5) :: zone
  character(len=10) :: time

  allocate(temp_material_id(99))
  allocate(temp_liq_pressure(99))
  allocate(temp_corr_comp_porosity(99))
  allocate(temp_corr_dcomp_porosity_dp(99))

  ! find the compressibility input filename
  filename_in = ''
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    auxvar => auxvars(ghosted_id)
    ! some code goes here!
  enddo

  !if (len(trim(auxvar%unittest_input_filename)) > 0) then
  if (1>0) then
  !------- an input file was provided -------------------------------
    input_file_given = PETSC_TRUE
    i = 1
    !open(action='read', file=trim(auxvar%unittest_input_filename), iostat=rc_in, &
    open(action='read', file=trim('./compressibility.in'), iostat=rc_in, &
         newunit=fu_in)
    read(fu_in, *) ! skip header line
    do
      read (fu_in, *, iostat=rc_in) temp_material_id(i), &
                                    temp_liq_pressure(i), &
                                    temp_corr_comp_porosity(i), &
                                    temp_corr_dcomp_porosity_dp(i)
      if (rc_in /= 0) exit 
      i = i + 1 
    enddo
    ! read in values while counting how many values with i
    allocate(material_id(i-1))
    allocate(liq_pressure(i-1))
    allocate(corr_compressed_porosity(i-1))
    allocate(corr_dcompressed_porosity_dp(i-1))
    material_id(:) = temp_material_id(1:i-1)
    liq_pressure(:) = temp_liq_pressure(1:i-1)
    corr_compressed_porosity(:) = temp_corr_comp_porosity(1:i-1)
    corr_dcompressed_porosity_dp(:) = temp_corr_dcomp_porosity_dp(1:i-1)
  else
  !------- an input file was not provided ---------------------------
    input_file_given = PETSC_FALSE
  endif

  prev_mat_id = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    mat_id = auxvars(ghosted_id)%id
    if (mat_id == prev_mat_id) cycle
    prev_mat_id = mat_id
    auxvar => auxvars(ghosted_id)

  !----- for creating input files only --------------------!
    call MaterialCompressSoilTest(2.5d5,auxvar)
    call MaterialCompressSoilTest(7.5d5,auxvar)
    call MaterialCompressSoilTest(9.0d5,auxvar)
    call MaterialCompressSoilTest(1.5d6,auxvar)
    call MaterialCompressSoilTest(3.0d6,auxvar)
    call MaterialCompressSoilTest(5.0d6,auxvar)
  !--------------------------------------------------------!

    write(id,'(I2)') auxvar%id
    filename_out = trim('./compressibility_id') // trim(id) // trim('.out')

    if (input_file_given) then
      !------- an input file was provided -------------------------------
      open(action='write', file=filename_out, iostat=rc_out, newunit=fu_out)
      call date_and_time(DATE=date,ZONE=zone,TIME=time,VALUES=values)
      write(fu_out,*) date(1:4),'/',date(5:6),'/',date(7:8),' ',time(1:2),':', &
                      time(3:4),' ',zone(1:3),':',zone(4:5),'UTC'
      write(fu_out,*)
      write(fu_out,'(a)') 'NOTE: The input file provided was:'
      !write(fu_out,'(a,a)') '      ', trim(auxvar%unittest_input_filename)
      write(fu_out,'(a,a)') '      ', trim('./compressibility.in')
      write(fu_out,'(a,d17.10,a)') 'NOTE: The validation test tolerance is ', &
                                   tolerance, '.'
      write(fu_out,*)

      i = 0
      do k=1,size(liq_pressure)
        if (material_id(k) == auxvar%id) then
          write(fu_out,'(a,I2,a)') '||-----------TEST-#',k,'-----------------------&
                                     &--------------||'
          write(fu_out,'(a)') '[in]  liquid pressure [Pa]:'
          write(fu_out,'(d17.10)') liq_pressure(k)

          write(fu_out,'(a)') '[out]  altered porosity [-]:'
          call MaterialCompressSoilBRAGFLO(auxvar,liq_pressure(k),compressed_porosity, &
                                           dcompressed_porosity_dp)
          write(fu_out,'(d17.10)') compressed_porosity
          write(fu_out,'(a)') '[correct] altered porosity [-]:'
          write(fu_out,'(d17.10)') corr_compressed_porosity(k)
          diff = abs(corr_compressed_porosity(k)-compressed_porosity)
          if (diff > (tolerance*corr_compressed_porosity(k))) then
            pass_fail = 'FAIL!'
            i = i + 1
          else
            pass_fail = 'pass'
          endif
          write(fu_out,'(a)') trim(pass_fail)

          write(fu_out,'(a)') '[out] d altered porosity dp [1/Pa]:'
          write(fu_out,'(d17.10)') dcompressed_porosity_dp
          write(fu_out,'(a)') '[correct] d altered porosity dp [1/Pa]:'
          write(fu_out,'(d17.10)') corr_dcompressed_porosity_dp(k)
          diff = abs(corr_dcompressed_porosity_dp(k)-dcompressed_porosity_dp)
          if (diff > (tolerance*corr_dcompressed_porosity_dp(k))) then
            pass_fail = 'FAIL!'
            i = i + 1
          else
            pass_fail = 'pass'
          endif
          write(fu_out,'(a)') trim(pass_fail)

          write(fu_out,*)
        endif
      enddo

      write(fu_out,'(a)') 'TEST SUMMARY:'
      if (i == 0) then
        write(fu_out,'(a)') ' All tests passed!'
      else
        write(fu_out,'(a,I3,a)') ' A total of (', i, ') test(s) failed!'
      endif
      close(fu_out)
    endif

    close(fu_in)

  enddo

  if (input_file_given) then
    deallocate(material_id)
    deallocate(liq_pressure)
    deallocate(corr_compressed_porosity)
    deallocate(corr_dcompressed_porosity_dp)
  endif
  deallocate(temp_material_id)
  deallocate(temp_liq_pressure)
  deallocate(temp_corr_comp_porosity)
  deallocate(temp_corr_dcomp_porosity_dp)

end subroutine MaterialCompressSoilUnitTest

! ************************************************************************** !

function MaterialAuxIndexToPropertyName(i)
  ! 
  ! Returns the name of the soil property associated with an index
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/06/16
  ! 
  
  implicit none

  PetscInt :: i

  character(len=MAXWORDLENGTH) :: MaterialAuxIndexToPropertyName

  if (i == soil_compressibility_index) then
    MaterialAuxIndexToPropertyName = 'soil compressibility'
  else if (i == soil_reference_pressure_index) then
    MaterialAuxIndexToPropertyName = 'soil reference pressure'
!  else if (i == soil_thermal_conductivity_index) then
!    MaterialAuxIndexToPropertyName = 'soil thermal conductivity'
!  else if (i == soil_heat_capacity_index) then
!    MaterialAuxIndexToPropertyName = 'soil heat capacity'
  else
    MaterialAuxIndexToPropertyName = 'unknown property'
  end if

end function MaterialAuxIndexToPropertyName

! ************************************************************************** !

subroutine MaterialAuxVarFractureStrip(fracture)
  ! 
  ! Deallocates a fracture auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/14/17
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none

  type(fracture_auxvar_type), pointer :: fracture

  if (.not.associated(fracture)) return

  ! properties and vector are now static arrays.
  deallocate(fracture)
  nullify(fracture)
  
end subroutine MaterialAuxVarFractureStrip

! ************************************************************************** !

subroutine MaterialAuxVarStrip(auxvar)
  ! 
  ! Deallocates a material auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none

  class(material_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%permeability)
  call DeallocateArray(auxvar%sat_func_prop)
  call DeallocateArray(auxvar%soil_properties)
  call MaterialAuxVarFractureStrip(auxvar%fracture)
  if (associated(auxvar%geomechanics_subsurface_prop)) then
    call DeallocateArray(auxvar%geomechanics_subsurface_prop)
  endif
  
end subroutine MaterialAuxVarStrip

! ************************************************************************** !

subroutine MaterialAuxDestroy(aux)
  ! 
  ! Deallocates a material auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/11
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  type(material_type), pointer :: aux
  
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call MaterialAuxVarStrip(aux%auxvars(iaux))
    enddo  
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)
    
  if (associated(aux%material_parameter)) then
    call DeallocateArray(aux%material_parameter%soil_heat_capacity)
    call DeallocateArray(aux%material_parameter%soil_thermal_conductivity)
  endif
  deallocate(aux%material_parameter)
  nullify(aux%material_parameter)
  
  deallocate(aux)
  nullify(aux)

end subroutine MaterialAuxDestroy

end module Material_Aux_class
