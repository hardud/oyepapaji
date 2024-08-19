[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
[]

[Mesh]
  file = plate_hole_3d.e
[]

[AuxVariables]
  [effective_plastic_strain]
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxKernels]
  [effective_plastic_strain]
    type = ADMaterialRealAux
    property = effective_plastic_strain
    variable = effective_plastic_strain
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_xy stress_zz strain_xx strain_yy strain_xy strain_zz plastic_strain_xx plastic_strain_yy plastic_strain_xy plastic_strain_zz vonmises_stress'
    incremental = true
    volumetric_locking_correction = true
    use_automatic_differentiation = true
  []
[]

[BCs]
  [symm_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = 4
    value = 0.0
  []
  [symm_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  []
  [fix_back_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = 6
    value = 0.0
  []
  [pull_right]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = 2
    function = 0.0 #0.001*t
  []
  [pull_top]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = 3
    function = 0.001*t
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.3
  []
#  [stress]
#    type = ComputeFiniteStrainElasticStress
#  []
  [stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'lsh'
  []
  [lsh]
    type = ADIsotropicPlasticityStressUpdate
    yield_stress = 1e3
    hardening_constant = 0.0
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-ksp_gmres_restart -pc_type'
  petsc_options_value = '101                lu'
  line_search = 'none'
  l_max_its = 50
  nl_max_its = 20
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  start_time = 0
  end_time = 1
  dt = 0.025

  [Predictor]
    type = SimplePredictor
    scale = 1.0
  []
[]

[Outputs]
  exodus = true
[]
