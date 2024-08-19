[Mesh]
    file = mesh-2.e
[]

[GlobalParams]
    displacements = 'disp_x disp_y'
[]

[Variables]
    # scale with one over Young's Modulus
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
    [./tdisp_top]
        type = ADFunctionDirichletBC
        variable = disp_y
        boundary = top
        function = 0.1*t
    [../]
    [./symmx_left]
        type = ADDirichletBC
        variable = disp_x
        boundary = left
        value = 0
    [../]
    [./tdisp_bottom]
        type = ADFunctionDirichletBC
        variable = disp_y
        boundary = bottom
        function = -0.1*t
    [../]
    [./symmy_right]
        type = ADDirichletBC
        variable = disp_x
        boundary = right
        value = 0
    [../]
[]

[Materials]
    [./elasticity]
        type = ADComputeIsotropicElasticityTensor
        poissons_ratio = 0.3
        youngs_modulus = 1e10
    [../]
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
    solve_type = 'Newton'

    petsc_options_iname = -pc_hypre_type
    petsc_options_value = boomeramg

    dtmin = 0.05
    num_steps = 5
[]

[Outputs]
    exodus = true
[]
