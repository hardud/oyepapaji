//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IsotropicPlasticityStressUpdate.h"
#include "MooseMesh.h"

/**
 * This class uses the Discrete material in a radial return isotropic plasticity
 * model.  This class is one of the basic radial return constitutive models;
 * more complex constitutive models combine creep and plasticity.
 *
 * This class models power law hardening by using the relation
 * \f$ \sigma = \sigma_y + K \epsilon^n \f$
 * where \f$ \sigma_y \f$ is the yield stress. This class solves for the yield
 * stress as the intersection of the power law relation curve and Hooke's law:
 * \f$ \epsilon_y = \frac{\sigma_y}{E} = \left( \frac{\sigma_y}{K} \right)^n \f$
 * where \f$epsilon_y \f$ is the total strain at the yield point and the stress
 * \f$ \sigma_y \f$ is the von Mises stress.
 * Parameters from the parent class, IsotropicPlasticityStressUpdate, are
 * suppressed to enable this class to solve for yield stress:
 * \f$ \sigma_y = \left( \frac{E^n}{K} \right)^{1/(n-1)} \f$
 */
template <bool is_ad>
class ZerothKinematicHardeningTempl : public IsotropicPlasticityStressUpdateTempl<is_ad>
{
public:
  static InputParameters validParams();

  ZerothKinematicHardeningTempl(const InputParameters & parameters);

  using Material::_qp;
  using RadialReturnStressUpdateTempl<is_ad>::_three_shear_modulus;
  using RadialReturnStressUpdateTempl<is_ad>::_base_name;

protected:
  virtual void
  computeStressInitialize(const GenericReal<is_ad> & effective_trial_stress,
                          const GenericRankFourTensor<is_ad> & elasticity_tensor) override;
  virtual void computeYieldStress(const GenericRankFourTensor<is_ad> & elasticity_tensor) override;
  virtual void computeBackStress(const GenericReal<is_ad> & plastic_strain_increment); // ADDED

  virtual GenericReal<is_ad> computeHardeningDerivative(const GenericReal<is_ad> & scalar) override;
  // virtual void initQpStatefulProperties() override;

  ///@{ Power law hardening coefficients
  Real _C;                         // Kinematic Hardening Modulus ADDED
  Real _K;                         // strength Coefficient
  Real _strain_hardening_exponent; // variable 'n' ADDED

  ///@}

  /// Elastic constants
  GenericReal<is_ad> _youngs_modulus;

  GenericReal<is_ad> _effective_trial_stress;    // ADDED
  GenericMaterialProperty<Real, is_ad> & _alpha; // ADDED
  const MaterialProperty<Real> & _alpha_old;     // ADDED

  GenericReal<is_ad> getIsotropicLameLambda(const GenericRankFourTensor<is_ad> & elasticity_tensor);
};

typedef ZerothKinematicHardeningTempl<false> ZerothKinematicHardening;
typedef ZerothKinematicHardeningTempl<true> ADZerothKinematicHardening;
