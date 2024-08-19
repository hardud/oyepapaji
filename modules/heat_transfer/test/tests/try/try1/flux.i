[Mesh]
file = mesh.e
[]

[Variables]
[Temperature]
order = FIRST
family = LAGRANGE
[]
[]

[Kernels]
[heat_transfer]
type = HeatConduction
variable = Temperature
[]

[]

[BCs]
[bottom]
type = NeumannBC
variable = Temperature
boundary = 'bottom'
[]
[top]
type = NeumannBC
variable = Temperature
boundary = 'top'
[]


[left]
type = DirichletBC
variable = Temperature
boundary = 'left'
value = 1700
[]

[right]
type = ConvectiveHeatFluxBC
variable = Temperature
boundary = 'right'
T_infinity = 273.0
heat_transfer_coefficient = 0
heat_transfer_coefficient_dT = 0
[]

[he]
type = DirichletBC
variable = Temperature
boundary = 'he'
value = 500
[]

[flow_channel]
type = DirichletBC
variable = Temperature
boundary = 'flow_channel'
value = 800
[]

[]
[Materials]
[w]
type = HeatConductionMaterial
block = 'W'
thermal_conductivity = 0.01
specific_heat = 0
[]

[steel]
type = HeatConductionMaterial
block = 'Steel'
thermal_conductivity = 0.02
specific_heat = 2
[]

[sic]
type = HeatConductionMaterial
block = 'SIC'
thermal_conductivity = 0.03
specific_heat = 3
[]

[]

[Executioner]
type = Steady
solve_type = 'PJFNK'
[]
[Outputs]
execute_on = 'timestep_end'
exodus = true
[]
