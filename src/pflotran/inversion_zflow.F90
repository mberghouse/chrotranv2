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

    ! For Wm
    PetscInt :: num_constraints_local    ! Number of constraints
    PetscInt :: num_constraints_total    ! Total number of constraints
    PetscInt, pointer :: rblock(:,:)     ! array stores info about reg.
    PetscReal, pointer :: Wm(:)          ! Regularization matrix

    !type(constrained_block_type), pointer :: constrained_block

  end type inversion_zflow_type




end module Inversion_ZFlow_class