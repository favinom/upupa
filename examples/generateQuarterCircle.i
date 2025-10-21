[Mesh]
  [ccmg]
    type = ConcentricCircleMeshGenerator
    num_sectors = ${ns}
    radii = '${radius}'
    rings = '5'
    has_outer_square = off
    pitch = 1.42063
    portion = top_right
    preserve_volumes = off
    smoothing_max_it = 200
[]
second_order = true
[]


[Variables]
  [u]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  file_base = quarterCircle
  xda = true
  # exodus = true
[]
