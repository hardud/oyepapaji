//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RadialReturnStressUpdate.h"
#include "libmesh/libmesh_common.h"
#include "StressUpdateBase.h"
#include "SingleVariableReturnMappingSolution.h"
#include "ADSingleVariableReturnMappingSolution.h"

/**
 * This class uses the Discrete material in a radial return isotropic plasticity
 * model.  This class is one of the basic radial return constitutive models;
 * more complex constitutive models combine creep and plasticity.
 *
 * This class inherits from RadialReturnStressUpdate and must be used
 * in conjunction with ComputeReturnMappingStress.  This class calculates
 * an effective trial stress, an effective scalar plastic strain
 * increment, and the derivative of the scalar effective plastic strain increment;
 * these values are passed to the RadialReturnStressUpdate to compute
 * the radial return stress increment.  This isotropic plasticity class also
 * computes the plastic strain as a stateful material property.
 *
 * This class is based on the implicit integration algorithm in F. Dunne and N.
 * Petrinic's Introduction to Computational Plasticity (2004) Oxford University
 * Press, pg. 146 - 149.
 */

template <bool is_ad>
class ZerothKinematicPlasticityStressUpdateTempl : public RadialReturnStressUpdateTempl<is_ad>
{
public:
  static InputParameters validParams();

  ZerothKinematicPlasticityStressUpdateTempl(const InputParameters & parameters);

  using Material::_qp;
  using RadialReturnStressUpdateTempl<is_ad>::_base_name;
  using RadialReturnStressUpdateTempl<is_ad>::_three_shear_modulus;

  using Material::_current_elem;
  using Material::_dt;
  using Material::_q_point;

  enum class SubsteppingType
  {
    NONE,
    ERROR_BASED,
    INCREMENT_BASED
  };

  virtual void
  computeStressInitialize(const GenericReal<is_ad> & _effective_trial_stress2,
                          const GenericRankFourTensor<is_ad> & elasticity_tensor) override;
  virtual GenericReal<is_ad> computeResidual(const GenericReal<is_ad> & _effective_trial_stress2,
                                             const GenericReal<is_ad> & scalar) override;
  virtual GenericReal<is_ad> computeDerivative(const GenericReal<is_ad> & _effective_trial_stress2,
                                               const GenericReal<is_ad> & scalar) override;

  virtual void computeYieldStress(const GenericRankFourTensor<is_ad> & elasticity_tensor);

  GenericReal<is_ad> yieldCondition() const { return _yield_condition; }

  virtual void
  computeStressFinalize(const GenericRankTwoTensor<is_ad> & plastic_strain_increment) override;

protected:
  virtual void updateState(GenericRankTwoTensor<is_ad> & strain_increment,
                           GenericRankTwoTensor<is_ad> & inelastic_strain_increment,
                           const GenericRankTwoTensor<is_ad> & /*rotation_increment*/,
                           GenericRankTwoTensor<is_ad> & stress_new,
                           const RankTwoTensor & /*stress_old*/,
                           const GenericRankFourTensor<is_ad> & elasticity_tensor,
                           const RankTwoTensor & elastic_strain_old,
                           bool compute_full_tangent_operator,
                           RankFourTensor & tangent_operator) override;
  using RadialReturnStressUpdateTempl<is_ad>::computeTangentOperator;
  using RadialReturnStressUpdateTempl<is_ad>::computeStressDerivative;
  virtual void initQpStatefulProperties() override;

  virtual void propagateQpStatefulProperties() override;
  // virtual Real computeStressDerivative(const Real effective_trial_stress2, const Real scalar);
  // {
  //   return 0.0;
  // }

  // virtual GenericReal<is_ad> computeHardeningValue(const GenericReal<is_ad> & scalar);
  // virtual GenericReal<is_ad> computeHardeningDerivative(const GenericReal<is_ad> & scalar);
  // virtual GenericRankTwoTensor<is_ad>
  // computeBackStress(const GenericRankTwoTensor<is_ad> & plastic_strain_increment); // ADDED
  const RankTwoTensor _deviatoric_stress;
  virtual void iterationFinalize(const GenericReal<is_ad> & scalar) override;

  Real _norm_dev_stress;
  Real _norm_dev_stress_squared;
  Real _effective_trial_stress2;
  Real scalar;
  // Real _three_shear_modulus;
  GenericMaterialProperty<Real, is_ad> & _effective_inelastic_strain;

  const MaterialProperty<Real> & _effective_inelastic_strain_old;
  GenericReal<is_ad> _effective_inelastic_strain_increment;

  Real _deriv;
  GenericReal<is_ad> _scalar_one;
  RankTwoTensor _flow_direction;
  RankFourTensor _flow_direction_dyad;
  RankFourTensor _deviatoric_projection_four;
  const bool _apply_strain;

  /// a string to prepend to the plastic strain Material Property name
  const std::string _plastic_prepend;

  const Function * _yield_stress_function;
  GenericReal<is_ad> _yield_stress;
  // const Real _hardening_constant;
  // const Function * const _hardening_function;

  GenericReal<is_ad> _yield_condition;
  GenericReal<is_ad> _hardening_slope;
  /// Debugging option to enable specifying instead of calculating strain

  /// Current value of scalar inelastic strain
  const GenericReal<is_ad> & effectiveInelasticStrainIncrement() const
  {
    return _effective_inelastic_strain_increment;
  }

  void updateEffectiveInelasticStrainIncrement(const GenericReal<is_ad> & eisi)
  {
    _effective_inelastic_strain_increment = eisi;
  }

  void updateEffectiveInelasticStrain(const GenericReal<is_ad> & increment)
  {
    _effective_inelastic_strain[_qp] = _effective_inelastic_strain_old[_qp] + increment;
  }
  /// plastic strain in this model
  GenericMaterialProperty<RankTwoTensor, is_ad> & _plastic_strain;

  /// old value of plastic strain
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;

  /// back stress in this model ADDED
  GenericMaterialProperty<RankTwoTensor, is_ad> & _backstress;

  /// old value of back stress ADDED
  const MaterialProperty<RankTwoTensor> & _backstress_old;

  ///  Kinematic Hardening Modulus ADDED
  Real _C;

  // GenericMaterialProperty<Real, is_ad> & _hardening_variable;
  // const MaterialProperty<Real> & _hardening_variable_old;
  const GenericVariableValue<is_ad> & _temperature;
};

typedef ZerothKinematicPlasticityStressUpdateTempl<false> ZerothKinematicPlasticityStressUpdate;
typedef ZerothKinematicPlasticityStressUpdateTempl<true> ADZerothKinematicPlasticityStressUpdate;
