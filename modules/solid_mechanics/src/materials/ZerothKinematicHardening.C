//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RadialReturnStressUpdate.h"
#include "libmesh/libmesh_common.h"
#include "ZerothKinematicHardening.h"
#include "ElasticityTensorTools.h"

registerMooseObject("SolidMechanicsApp", ZerothKinematicHardening);
registerMooseObject("SolidMechanicsApp", ADZerothKinematicHardening);

template <bool is_ad>
InputParameters
ZerothKinematicHardeningTempl<is_ad>::validParams()
{
  InputParameters params = IsotropicPlasticityStressUpdateTempl<is_ad>::validParams();
  params.addClassDescription("This class uses the discrete material in a radial return isotropic "
                             "plasticity power law hardening model, solving for the yield stress "
                             "as the intersection of the power law relation curve and Hooke's law. "
                             " This class can be used in conjunction with other creep and "
                             "plasticity materials for more complex simulations.");

  // Set and Suppress parameters to enable calculation of the yield stress
  params.set<Real>("yield_stress") = 1.0;
  params.set<Real>("hardening_constant") = 1.0;
  params.suppressParameter<Real>("yield_stress");
  params.suppressParameter<Real>("hardening_constant");

  // Power law hardening specific parameters
  params.addRequiredParam<Real>("strength_coefficient",
                                "The strength coefficient (K) for power law hardening");
  params.addRequiredRangeCheckedParam<Real>(
      "strain_hardening_exponent",
      "strain_hardening_exponent>=0.0 & strain_hardening_exponent <=1.0",
      "The strain hardening exponent (n) for power law hardening");

  // Kinematic hardening specific parameters ADDED
  params.addRequiredParam<Real>("kinematic_hardening_modulus",
                                "The kinematic hardening modulus (C) for back stress evolution");

  return params;
}

template <bool is_ad>
ZerothKinematicHardeningTempl<is_ad>::ZerothKinematicHardeningTempl(
    const InputParameters & parameters)
  : IsotropicPlasticityStressUpdateTempl<is_ad>(parameters),
    _C(parameters.get<Real>("kinematic_hardening_modulus")), // ADDED
    _K(parameters.get<Real>("strength_coefficient")),
    _strain_hardening_exponent(parameters.get<Real>("strain_hardening_exponent")),
    _alpha(this->template declareGenericProperty<Real, is_ad>("alpha")), // ADDED
    _alpha_old(this->template getMaterialPropertyOld<Real>("alpha"))     // ADDED
{
}

template <bool is_ad>
void
ZerothKinematicHardeningTempl<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  RadialReturnStressUpdateTempl<is_ad>::computeStressInitialize(effective_trial_stress,
                                                                elasticity_tensor);

  computeYieldStress(elasticity_tensor);

  _effective_trial_stress = effective_trial_stress;
  this->_yield_condition =
      effective_trial_stress - this->_alpha[_qp] - this->_yield_stress; // ADDED

  this->_alpha[_qp] = this->_alpha_old[_qp];
  this->_plastic_strain[_qp] = this->_plastic_strain_old[_qp]; // ADDED
}

template <bool is_ad>
GenericReal<is_ad>
ZerothKinematicHardeningTempl<is_ad>::computeHardeningDerivative(const GenericReal<is_ad> & scalar)
{
  const GenericReal<is_ad> stress_delta = _effective_trial_stress - _three_shear_modulus * scalar;
  GenericReal<is_ad> slope = std::pow(stress_delta, (1.0 / _strain_hardening_exponent - 1.0)) /
                             _strain_hardening_exponent * 1.0 /
                             std::pow(_K, 1.0 / _strain_hardening_exponent);
  slope -= 1.0 / _youngs_modulus;

  return 1.0 / slope;
}

template <bool is_ad>
void
ZerothKinematicHardeningTempl<is_ad>::computeYieldStress(
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  // Pull in the Lam\`{e} lambda, and caculate E
  const GenericReal<is_ad> lambda = getIsotropicLameLambda(elasticity_tensor);
  const GenericReal<is_ad> shear_modulus = _three_shear_modulus / 3.0;

  _youngs_modulus = shear_modulus * (3.0 * lambda + 2 * shear_modulus) / (lambda + shear_modulus);

  // Then solve for yield stress using equation from the header file
  this->_yield_stress = std::pow(_K / std::pow(_youngs_modulus, _strain_hardening_exponent),
                                 1.0 / (1.0 - _strain_hardening_exponent));
  if (this->_yield_stress <= 0.0)
    mooseError("The yield stress must be greater than zero, but during the simulation your yield "
               "stress became less than zero.");
}

template <bool is_ad> // ADDED
void
ZerothKinematicHardeningTempl<is_ad>::computeBackStress(
    const GenericReal<is_ad> & plastic_strain_increment)
{
  this->_alpha[_qp] = this->_alpha_old[_qp] +
                      _C * plastic_strain_increment *
                          std::pow((this->_yield_stress / _K), (1.0 / _strain_hardening_exponent));
}

template <bool is_ad>
GenericReal<is_ad>
ZerothKinematicHardeningTempl<is_ad>::getIsotropicLameLambda(
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  const GenericReal<is_ad> lame_lambda = elasticity_tensor(0, 0, 1, 1);

  if (this->_mesh.dimension() == 3 &&
      MetaPhysicL::raw_value(lame_lambda) != MetaPhysicL::raw_value(elasticity_tensor(1, 1, 2, 2)))
    mooseError(
        "Check to ensure that your Elasticity Tensor is truly Isotropic: dixfferent lambda values");
  return lame_lambda;
}

template class ZerothKinematicHardeningTempl<false>;
template class ZerothKinematicHardeningTempl<true>;
