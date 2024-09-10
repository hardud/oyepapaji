//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ZerothKinematicPlasticityStressUpdate.h"

#include "Function.h"
#include "ElasticityTensorTools.h"
#include <algorithm>

registerMooseObject("SolidMechanicsApp", ADZerothKinematicPlasticityStressUpdate);
registerMooseObject("SolidMechanicsApp", ZerothKinematicPlasticityStressUpdate);

template <bool is_ad>
InputParameters
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnStressUpdateTempl<is_ad>::validParams();
  params.addClassDescription("This class uses the discrete material in a radial return isotropic "
                             "plasticity model.  This class is one of the basic radial return "
                             "constitutive models, yet it can be used in conjunction with other "
                             "creep and plasticity materials for more complex simulations.");
  params.addParam<FunctionName>("yield_stress_function",
                                "Yield stress as a function of temperature");
  params.addParam<Real>("yield_stress", "The point at which plastic strain begins accumulating");

  params.addCoupledVar("temperature", 0.0, "Coupled Temperature");
  params.addDeprecatedParam<std::string>(
      "plastic_prepend",
      "",
      "String that is prepended to the plastic_strain Material Property",
      "This has been replaced by the 'base_name' parameter");
  params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";
  params.addRequiredParam<Real>("kinematic_hardening_modulus",
                                "The kinematic hardening modulus (C) for back stress evolution");
  params.addRequiredParam<Real>(
      "material_constant_gamma",
      "The nonlinear hardening parameter (gamma) for back stress evolution");
  return params;
}

template <bool is_ad>
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::ZerothKinematicPlasticityStressUpdateTempl(
    const InputParameters & parameters)
  : RadialReturnStressUpdateTempl<is_ad>(parameters),

    _effective_inelastic_strain(this->template declareGenericProperty<Real, is_ad>(
        this->_base_name +
        this->template getParam<std::string>("effective_inelastic_strain_name"))),
    _effective_inelastic_strain_old(this->template getMaterialPropertyOld<Real>(
        this->_base_name +
        this->template getParam<std::string>("effective_inelastic_strain_name"))),

    _apply_strain(this->template getParam<bool>("apply_strain")),
    _plastic_prepend(this->template getParam<std::string>("plastic_prepend")),

    _yield_stress_function(this->isParamValid("yield_stress_function")
                               ? &this->getFunction("yield_stress_function")
                               : nullptr),
    _yield_stress(this->isParamValid("yield_stress") ? this->template getParam<Real>("yield_stress")
                                                     : 0),

    _yield_condition(-1.0), // set to a non-physical value to catch uninitalized yield condition
    _hardening_slope(0.0),

    _plastic_strain(this->template declareGenericProperty<RankTwoTensor, is_ad>(
        _base_name + _plastic_prepend + "plastic_strain")),
    _plastic_strain_old(this->template getMaterialPropertyOld<RankTwoTensor>(
        _base_name + _plastic_prepend + "plastic_strain")),

    _backstress(this->template declareGenericProperty<RankTwoTensor, is_ad>("backstress")), // ADDED

    _backstress_old(this->template getMaterialPropertyOld<RankTwoTensor>("backstress")),
    _C(parameters.get<Real>("kinematic_hardening_modulus")), // ADDED // added to remove errors
    _gamma(parameters.get<Real>("material_constant_gamma")),

    _temperature(this->template coupledGenericValue<is_ad>("temperature"))

{
  if (parameters.isParamSetByUser("yield_stress") && _yield_stress <= 0.0)
    mooseError("Yield stress must be greater than zero");

  // Both of these parameters are given default values by derived classes, which makes them valid
  if (_yield_stress_function == nullptr && !this->isParamValid("yield_stress"))
    mooseError("Either yield_stress or yield_stress_function must be given");
}

template <bool is_ad>
void
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::updateState(
    GenericRankTwoTensor<is_ad> & strain_increment,
    GenericRankTwoTensor<is_ad> & inelastic_strain_increment,
    const GenericRankTwoTensor<is_ad> & /*rotation_increment*/,
    GenericRankTwoTensor<is_ad> & stress_new,
    const RankTwoTensor & /*stress_old*/,
    const GenericRankFourTensor<is_ad> & elasticity_tensor,
    const RankTwoTensor & elastic_strain_old,
    bool compute_full_tangent_operator,
    RankFourTensor & tangent_operator)
{

  GenericRankTwoTensor<is_ad> deviatoric_trial_stress = stress_new.deviatoric();

  GenericRankTwoTensor<is_ad> relative_stress = deviatoric_trial_stress - _backstress[_qp];

  GenericReal<is_ad> rel_stress_squared = relative_stress.doubleContraction(relative_stress);
  GenericReal<is_ad> norm_rel_stress = std::sqrt(rel_stress_squared);
  GenericReal<is_ad> _effective_trial_stress2 =
      MetaPhysicL::raw_value(rel_stress_squared) ? std::sqrt(3.0 / 2.0 * rel_stress_squared) : 0.0;

  computeStressInitialize(_effective_trial_stress2, elasticity_tensor);

  // Ensure shear modulus is set
  mooseAssert(
      _three_shear_modulus != 0.0,
      "Shear modulus is zero. Ensure that the base class computeStressInitialize() is called.");

  // Use Newton iteration to determine the scalar effective inelastic strain increment
  _effective_inelastic_strain_increment = 0.0;
  if (!MooseUtils::absoluteFuzzyEqual(_effective_trial_stress2, 0.0))
  {
    this->returnMappingSolve(
        _effective_trial_stress2, _effective_inelastic_strain_increment, this->_console);
    if (_effective_inelastic_strain_increment != 0.0)
    {
      inelastic_strain_increment =
          deviatoric_trial_stress *
          (1.5 * _effective_inelastic_strain_increment / _effective_trial_stress2);

      _scalar_one = _three_shear_modulus * _effective_inelastic_strain_increment / std::sqrt(1.5) /
                    norm_rel_stress;
    }
    else
    {
      inelastic_strain_increment.zero();
      _scalar_one = 0;
    }
  }
  else
  {
    inelastic_strain_increment.zero();
    _scalar_one = 0;
  }

  // Apply strain if applicable
  if (_apply_strain)
  {
    strain_increment -= inelastic_strain_increment;
    updateEffectiveInelasticStrain(_effective_inelastic_strain_increment);

    // Update stress state
    stress_new = elasticity_tensor * (strain_increment + elastic_strain_old);
  }

  // Finalize the stress computation
  computeStressFinalize(inelastic_strain_increment);

  // Compute the tangent operator if necessary
  if constexpr (!is_ad)
  {
    if (compute_full_tangent_operator)
    {
      computeTangentOperator(_effective_trial_stress2, stress_new, tangent_operator);
    }
  }
  else
  {
    libmesh_ignore(compute_full_tangent_operator);
    libmesh_ignore(tangent_operator);
  }
}

template <bool is_ad>
void
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::initQpStatefulProperties()
{
  // _hardening_variable[_qp] = 0.0;
  _backstress[_qp] = 0.0;
  _plastic_strain[_qp].zero();
}

template <bool is_ad>
void
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::propagateQpStatefulProperties()
{
  // _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
  _backstress[_qp] = _backstress_old[_qp];

  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
void
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & _effective_trial_stress2,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  RadialReturnStressUpdateTempl<is_ad>::computeStressInitialize(_effective_trial_stress2,
                                                                elasticity_tensor);

  computeYieldStress(elasticity_tensor);

  _yield_condition = _effective_trial_stress2 /*_backstress[_qp]*/ - _yield_stress;

  // std::cout << "Value of _backstress[_qp] before assignment: " << _backstress[_qp] << std::endl;
  std::cout << "Effective Trial Stress: " << _effective_trial_stress2 << std::endl;
  std::cout << "Yield Condition: " << _yield_condition << std::endl;
  _backstress[_qp] = _backstress_old[_qp];

  // std::cout << "Value of _backstress[_qp] after assignment: " << _backstress[_qp] << std::endl;
  //  // updated backstress values

  _plastic_strain[_qp] = _plastic_strain_old[_qp];
}

template <bool is_ad>
GenericReal<is_ad>
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeResidual(
    const GenericReal<is_ad> & _effective_trial_stress2, const GenericReal<is_ad> & _scalar_one)
{
  mooseAssert(_yield_condition != -1.0,
              "the yield stress was not updated by computeStressInitialize");

  if (_yield_condition > 0.0)
  {

    _backstress[_qp] = _backstress_old[_qp];
    _plastic_strain[_qp] = _plastic_strain_old[_qp];

    GenericReal<is_ad> residual =
        (_effective_trial_stress2 - _yield_stress) / _three_shear_modulus - _scalar_one;

    return residual;
  }

  return 0.0;
}

template <bool is_ad>
GenericReal<is_ad>
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeDerivative(
    const GenericReal<is_ad> & /*_effective_trial_stress2*/,
    const GenericReal<is_ad> & /*_scalar_one*/)
{
  if (_yield_condition > 0.0)
    return -1.0 - _hardening_slope / _three_shear_modulus;

  return 1.0;
}

template <bool is_ad>
void
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::iterationFinalize(
    const GenericReal<is_ad> & /*_scalar_one*/)
{
  if (_yield_condition > 0.0)
  {

    _plastic_strain[_qp] = _plastic_strain_old[_qp];

    _backstress[_qp] = _backstress_old[_qp];
  }
}

template <bool is_ad>
void
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _plastic_strain[_qp] += plastic_strain_increment;
  _backstress[_qp] = _backstress_old[_qp] + (2.0 / 3.0) * _C * plastic_strain_increment -
                     _gamma * _backstress[_qp] * _effective_inelastic_strain_increment;
  std::cout << "Backstress: " << MetaPhysicL::raw_value(_backstress[_qp]) << std::endl;
  std::cout << "Plastic_strain: " << MetaPhysicL::raw_value(_plastic_strain[_qp]) << std::endl;
}

template <bool is_ad>
void
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeYieldStress(
    const GenericRankFourTensor<is_ad> & /*elasticity_tensor*/)
{

  if (_yield_stress_function)
  {
    static const Moose::GenericType<Point, is_ad> p;
    _yield_stress = _yield_stress_function->value(_temperature[_qp], p);

    if (_yield_stress <= 0.0)
      mooseError("In ",
                 this->_name,
                 ": The calculated yield stress (",
                 _yield_stress,
                 ") is less than zero");
  }
}

template class ZerothKinematicPlasticityStressUpdateTempl<false>;
template class ZerothKinematicPlasticityStressUpdateTempl<true>;
