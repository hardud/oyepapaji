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
  // Linear strain hardening parameters
  params.addParam<FunctionName>("yield_stress_function",
                                "Yield stress as a function of temperature");
  params.addParam<Real>("yield_stress", "The point at which plastic strain begins accumulating");
  // params.addParam<FunctionName>("hardening_function",
  //                               "True stress as a function of plastic strain");
  // params.addParam<Real>("hardening_constant", "Hardening slope");
  params.addCoupledVar("temperature", 0.0, "Coupled Temperature");
  params.addDeprecatedParam<std::string>(
      "plastic_prepend",
      "",
      "String that is prepended to the plastic_strain Material Property",
      "This has been replaced by the 'base_name' parameter");
  params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";
  params.addRequiredParam<Real>("kinematic_hardening_modulus",
                                "The kinematic hardening modulus (C) for back stress evolution");

  return params;
}

template <bool is_ad>
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::ZerothKinematicPlasticityStressUpdateTempl(
    const InputParameters & parameters)
  : RadialReturnStressUpdateTempl<is_ad>(parameters),

    // ADDED

    // _hardening_constant(this->isParamValid("hardening_constant")
    //                         ? this->template getParam<Real>("hardening_constant")
    //                         : 0),

    // _hardening_function(this->isParamValid("hardening_function")
    //                         ? &this->getFunction("hardening_function")
    //                         : nullptr),

    _effective_inelastic_strain(this->template declareGenericProperty<Real, is_ad>(
        this->_base_name +
        this->template getParam<std::string>("effective_inelastic_strain_name"))),
    _effective_inelastic_strain_old(this->template getMaterialPropertyOld<Real>(
        this->_base_name +
        this->template getParam<std::string>("effective_inelastic_strain_name"))),
    // _hardening_slope(0.0),
    _apply_strain(this->template getParam<bool>("apply_strain")),
    _plastic_prepend(this->template getParam<std::string>("plastic_prepend")),

    _yield_stress_function(this->isParamValid("yield_stress_function")
                               ? &this->getFunction("yield_stress_function")
                               : nullptr),
    _yield_stress(this->isParamValid("yield_stress") ? this->template getParam<Real>("yield_stress")
                                                     : 0),

    _yield_condition(-1.0), // set to a non-physical value to catch uninitalized yield condition
    _plastic_strain(this->template declareGenericProperty<RankTwoTensor, is_ad>(
        _base_name + _plastic_prepend + "plastic_strain")),
    _plastic_strain_old(this->template getMaterialPropertyOld<RankTwoTensor>(
        _base_name + _plastic_prepend + "plastic_strain")),

    _backstress(this->template declareGenericProperty<RankTwoTensor, is_ad>("backstress")), // ADDED

    _backstress_old(this->template getMaterialPropertyOld<RankTwoTensor>("backstress")),
    _C(parameters.get<Real>("kinematic_hardening_modulus")), // ADDED // added to remove errors

    // _hardening_variable(
    //     this->template declareGenericProperty<Real, is_ad>(_base_name + "hardening_variable")),
    // _hardening_variable_old(
    //     this->template getMaterialPropertyOld<Real>(_base_name + "hardening_variable")),

    _temperature(this->template coupledGenericValue<is_ad>("temperature"))

{
  if (parameters.isParamSetByUser("yield_stress") && _yield_stress <= 0.0)
    mooseError("Yield stress must be greater than zero");

  // Both of these parameters are given default values by derived classes, which makes them valid
  if (_yield_stress_function == nullptr && !this->isParamValid("yield_stress"))
    mooseError("Either yield_stress or yield_stress_function must be given");
  // if (!parameters.isParamValid("hardening_constant") &&
  // !this->isParamValid("hardening_function"))
  //   mooseError("Either hardening_constant or hardening_function must be defined");

  // if (parameters.isParamSetByUser("hardening_constant") &&
  // this->isParamValid("hardening_function"))
  //   mooseError(
  //       "Only the hardening_constant or only the hardening_function can be defined but not
  //       both");
}

// template<bool is_ad>
// class RadialReturnStressUpdateTempl : public StressUpdateBaseTempl<bool is_ad>
// {

// protected:
//     GenericRankTwoTensor<is_ad> _deviatoric_stress;
//     Real _norm_dev_stress;
//     Real _norm_dev_stress_squared;
//     Real _effective_trial_stress;
//     Real _three_shear_modulus;
//     Real _effective_inelastic_strain_increment;
//     Real _deriv;
//     Real _scalar_one;
//     RankTwoTensor _flow_direction;
//     RankFourTensor _flow_direction_dyad;
//     RankFourTensor _deviatoric_projection_four;
// }

// template <bool is_ad>
// void
// ZerothKinematicPlasticityStressUpdateTempl<is_ad>::updateState(
//     GenericRankTwoTensor<is_ad> & strain_increment,
//     GenericRankTwoTensor<is_ad> & inelastic_strain_increment,
//     const GenericRankTwoTensor<is_ad> & /*rotation_increment*/,
//     GenericRankTwoTensor<is_ad> & stress_new,
//     const RankTwoTensor & /*stress_old*/,
//     const GenericRankFourTensor<is_ad> & elasticity_tensor,
//     const RankTwoTensor & elastic_strain_old,
//     bool compute_full_tangent_operator,
//     RankFourTensor & tangent_operator)

// {
//   // GenericRankTwoTensor<is_ad> stress_new; // Added

//   const RankTwoTensor _deviatoric_stress = stress_new.deviatoric();
//   _norm_dev_stress_squared = _deviatoric_stress.doubleContraction(_deviatoric_stress);
//   _norm_dev_stress = std::sqrt(_norm_dev_stress_squared);

//   //  other necessary quantities

//   _flow_direction = _deviatoric_stress / _norm_dev_stress;
//   _flow_direction_dyad = _flow_direction.outerProduct(_flow_direction);

//   _deriv = computeStressDerivative(_effective_trial_stress2,
//   _effective_inelastic_strain_increment); _scalar_one = _three_shear_modulus *
//   _effective_inelastic_strain_increment / std::sqrt(1.5) /
//                 _norm_dev_stress;
//   tangent_operator = _scalar_one * _deviatoric_projection_four +
//                      (_three_shear_modulus * _deriv - _scalar_one) * _flow_direction_dyad;
//   GenericRankTwoTensor<is_ad> stress_new;
//   RadialReturnStressUpdateTempl<is_ad>::updateState(strain_increment,
//                                                     inelastic_strain_increment,
//                                                     /*rotation_increment*/
//                                                     stress_new,
//                                                     /*stress_old*/
//                                                     elasticity_tensor,
//                                                     elastic_strain_old,
//                                                     compute_full_tangent_operator,
//                                                     tangent_operator);

//   GenericRankTwoTensor<is_ad> deviatoric_trial_stress = stress_new.deviatoric();

//   GenericRankTwoTensor<is_ad> computeBackStress(
//       const GenericRankTwoTensor<is_ad> & plastic_strain_increment);

//   GenericRankTwoTensor<is_ad> relative_stress = deviatoric_trial_stress - _backstress[_qp];

//   GenericReal<is_ad> rel_stress_squared = relative_stress.doubleContraction(relative_stress);

//   GenericReal<is_ad> effective_trial_stress2 =
//       MetaPhysicL::raw_value(rel_stress_squared) ? std::sqrt(3.0 / 2.0 * rel_stress_squared) :
//       0.0;

//   computeStressInitialize(effective_trial_stress2, elasticity_tensor);

//   mooseAssert(
//       _three_shear_modulus != 0.0,
//       "Shear modulus is zero. Ensure that the base class computeStressInitialize() is called.");

//   // Use Newton iteration to determine the scalar effective inelastic strain increment
//   _effective_inelastic_strain_increment = 0.0;
//   if (!MooseUtils::absoluteFuzzyEqual(effective_trial_stress2, 0.0))
//   {
//     this->returnMappingSolve(
//         effective_trial_stress2, _effective_inelastic_strain_increment, this->_console);
//     if (_effective_inelastic_strain_increment != 0.0)
//       inelastic_strain_increment =
//           deviatoric_trial_stress *
//           (1.5 * _effective_inelastic_strain_increment / effective_trial_stress2);
//     else
//       inelastic_strain_increment.zero();
//   }
//   else
//     inelastic_strain_increment.zero();

//   if (_apply_strain)
//   {
//     strain_increment -= inelastic_strain_increment;
//     updateEffectiveInelasticStrain(_effective_inelastic_strain_increment);

//     // Use the old elastic strain here because we require tensors used by this class
//     // to be isotropic and this method natively allows for changing in time
//     // elasticity tensors
//     stress_new = elasticity_tensor * (strain_increment + elastic_strain_old);
//   }

//   computeStressFinalize(inelastic_strain_increment);

//   if constexpr (!is_ad)
//   {
//     if (compute_full_tangent_operator)
//       computeTangentOperator(effective_trial_stress2, stress_new, tangent_operator);
//   }
//   else
//   {
//     libmesh_ignore(compute_full_tangent_operator);
//     libmesh_ignore(tangent_operator);
//   }
// }

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
  // Compute the deviatoric trial stress and related quantities
  GenericRankTwoTensor<is_ad> deviatoric_trial_stress = stress_new.deviatoric();
  GenericRankTwoTensor<is_ad> computeBackStress(
      const GenericRankTwoTensor<is_ad> & plastic_strain_increment);
  GenericRankTwoTensor<is_ad> relative_stress = deviatoric_trial_stress - _backstress[_qp];

  GenericReal<is_ad> rel_stress_squared = relative_stress.doubleContraction(relative_stress);

  GenericReal<is_ad> _effective_trial_stress2 =
      MetaPhysicL::raw_value(rel_stress_squared) ? std::sqrt(3.0 / 2.0 * rel_stress_squared) : 0.0;

  // GenericReal<is_ad> dev_trial_stress_squared =
  //     deviatoric_trial_stress.doubleContraction(deviatoric_trial_stress);

  // GenericReal<is_ad> _effective_trial_stress2 =
  //     MetaPhysicL::raw_value(dev_trial_stress_squared)
  //         ? std::sqrt(3.0 / 2.0 * dev_trial_stress_squared)
  //         : 0.0;

  // Initialize stress computation
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
    }
    else
    {
      inelastic_strain_increment.zero();
    }
  }
  else
  {
    inelastic_strain_increment.zero();
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

  _yield_condition = _effective_trial_stress2 - /*_backstress[_qp]*/ -_yield_stress;

  // std::cout << "Value of _backstress[_qp] before assignment: " << _backstress[_qp] << std::endl;

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
    // _hardening_slope = computeHardeningDerivative(scalar);
    // _hardening_variable[_qp] = computeHardeningValue(scalar);

    return (_effective_trial_stress2 - /*_backstress[_qp]*/ -_yield_stress) / _three_shear_modulus -
           _scalar_one;
  }

  return 0.0;
}

template <bool is_ad>
GenericReal<is_ad>
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeDerivative(
    const GenericReal<is_ad> & /*_effective_trial_stress2*/, const GenericReal<is_ad> & /*scalar*/)
{
  // if (_yield_condition > 0.0)
  //   return -1.0 - _hardening_slope / _three_shear_modulus;

  return 1.0;
}

// template <bool is_ad>
// void
// ZerothKinematicPlasticityStressUpdateTempl<is_ad>::iterationFinalize(
//     const GenericReal<is_ad> & scalar)
// {
//   if (_yield_condition > 0.0)
//     _hardening_variable[_qp] = computeHardeningValue(scalar);
// }

template <bool is_ad>
void
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  _plastic_strain[_qp] += plastic_strain_increment;
}

// template <bool is_ad>
// GenericReal<is_ad>
// ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeHardeningValue(
//     const GenericReal<is_ad> & scalar)
// {
//   if (_hardening_function)
//   {
//     const Real strain_old = this->_effective_inelastic_strain_old[_qp];
//     return _hardening_function->value(strain_old + scalar) - _yield_stress;
//   }

//   return _hardening_variable_old[_qp] + _hardening_slope * scalar;
// }

// template <bool is_ad>
// GenericReal<is_ad>
// ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeHardeningDerivative(
//     const GenericReal<is_ad> & /*scalar*/)
// {
//   if (_hardening_function)
//   {
//     const Real strain_old = this->_effective_inelastic_strain_old[_qp];
//     return _hardening_function->timeDerivative(strain_old);
//   }

//   return _hardening_constant;
// }

template <bool is_ad>
void
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeYieldStress(
    const GenericRankFourTensor<is_ad> & /*elasticity_tensor*/)
{
  // std::cout << "Yield Condition: " << _yield_condition << std::endl;

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

template <bool is_ad> // ADDED
GenericRankTwoTensor<is_ad>
ZerothKinematicPlasticityStressUpdateTempl<is_ad>::computeBackStress(
    const GenericRankTwoTensor<is_ad> & plastic_strain_increment)
{
  return this->_backstress[_qp] = this->_backstress_old[_qp] + _C * plastic_strain_increment;
}

template class ZerothKinematicPlasticityStressUpdateTempl<false>;
template class ZerothKinematicPlasticityStressUpdateTempl<true>;
