[Mesh]
 [sphere]
  type = SphereMeshGenerator
  radius = 0.2 #${radius}
  nr = 3 #${ns}
  #elem_type = HEX27
 []
 second_order = true

  [./central_boundary]
    type = SideSetsAroundSubdomainGenerator
    input = sphere
    block = 0
    new_boundary = 7
  []


  [SubdomainBoundingBoxBottom]
    type = SubdomainBoundingBoxGenerator
    input = central_boundary
    block_id = 1
    bottom_left = '-99999999 -99999999 -99999999'
    top_right = '0 99999999 99999999'
  []

  [ed0]
    type = BlockDeletionGenerator
    input = SubdomainBoundingBoxBottom
    block = 1
    new_boundary = bottom
  []

  [SubdomainBoundingBoxLeft]
    type = SubdomainBoundingBoxGenerator
    input = ed0
    block_id = 1
    bottom_left = '-99999999 -99999999 -99999999'
    top_right = '99999999  0 99999999'
  []

  [ed1]
    type = BlockDeletionGenerator
    input = SubdomainBoundingBoxLeft
    block = 1
    new_boundary = left
  []

  [SubdomainBoundingBoxBack]
    type = SubdomainBoundingBoxGenerator
    input = ed1
    block_id = 1
    bottom_left = '-99999999 -99999999 -99999999'
    top_right = '99999999 99999999 0'
  []

  [ed2]
    type = BlockDeletionGenerator
    input = SubdomainBoundingBoxBack
    block = 1
    new_boundary = back
  []


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
  file_base = octSphere
  xda = true
#  exodus = true
[]
