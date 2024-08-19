[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 4
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  ny = 4
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  # scale with one over Young's modulus
  [./disp_x]
    scaling = 1e-10
  [../]
  [./disp_y]
    scaling = 1e-10
  [../]
[]

[Kernels]
  [./stress_x]
    type = ADStressDivergenceTensors
    component = 0
    variable = disp_x
    use_displaced_mesh = true
  [../]
  [./stress_y]
    type = ADStressDivergenceTensors
    component = 1
    variable = disp_y
    use_displaced_mesh = true
  [../]
[]

[BCs]
  [./symmx_right]
    type = ADDirichletBC
    variable = disp_y
    boundary = right
    value = 0
  [../]
  [./symmx_left]
    type = ADDirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./tdisp_top]
    type = ADDirichletBC
    variable = disp_y
    boundary = top
    value = 0.1
  [../]
  [./tdisp_bottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.1
  [../]
[]

[Materials]
  [./elasticity]
    type = ADComputeIsotropicElasticityTensor
    poissons_ratio = 0.3
    youngs_modulus = 1e10
  [../]
[]

[Materials]
  [./strain]
    type = ADComputeFiniteStrain
  [../]
  [./stress]
    type = ADComputeFiniteStrainElasticStress
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.5
  solve_type = 'NEWTON'

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomeramg

  dtmin = 0.05
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
