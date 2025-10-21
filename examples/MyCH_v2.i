[Mesh]

#[gmg]
 type = GeneratedMesh  # DistributedRectilinearMeshGenerator # GeneratedMeshGenerator
 dim = 2
 nx = 80
 ny = 80
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
 elem_type = QUAD9
#[]
[]

[Variables]
[disp_x] order = SECOND []
[disp_y] order = SECOND []
[p] order = FIRST []
[J_k][]
[b_k][]
[]



[Materials]
[./material] type = NeoHookeanPlasticityGradeZero mu = 1e2 lambda = 2.0e2 disp_x = disp_x disp_y = disp_y []
[./materialms] type = GradeOneGrowth lambda_p = 0 J_k = J_k []
[]

[Kernels]
[./linearelasticityx] type = PoroElastoPlasticity variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y p = p [../]
[./linearelasticityy] type = PoroElastoPlasticity variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y p = p [../]


[./pressure] type = PoroElastoPlasticityMassBalance_v2 variable = p permeability = 0.001 solid_fraction = 0.6 disp_x = disp_x disp_y = disp_y p = p J_k = J_k [../]
[./pressureflux] type = DivergenceChemicalPotential variable = p  m = 4e-5  disp_x = disp_x  disp_y = disp_y  J_k = J_k  b_k = b_k [../]


[./timeevolutionjk] type = StandardTimeDerivative variable = J_k [../]
[./divchempot] type = DivergenceChemicalPotential variable = J_k  m = 1e-4  disp_x = disp_x  disp_y = disp_y  J_k = J_k  b_k = b_k [../]


[./chempotdefinizione] type = ChemicalPotential  variable = b_k  d = 1  kv = 1  J_k = J_k  b_k = b_k [../]
[./growthcoupling] type = CouplingGrowthMechanics  variable = b_k J_k = J_k  b_k = b_k disp_x = disp_x disp_y = disp_y [../]
[]

[ICs]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1-0.5*cos(pi*x)*cos(pi*y)  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 0.9+0.2*exp(-10*(x^2+y^2))  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1+0.3*(x>0.899)*(y>0.899)  [../]
[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1+0.1*(y<0.199)  [../]
[]

[BCs]
[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = FunctionNeumannBC variable = disp_y function = -0.1*t*(t<1) boundary = top  [../]
[./topp] type = DirichletBC variable = p value = 0.0 boundary = top  [../]

#[./leftp] type = DirichletBC variable = p value = 0.0 boundary = left  [../]

#[./rightp] type = DirichletBC variable = p value = 0.0 boundary = right  [../]

[./bottomx] type = DirichletBC variable = disp_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = DirichletBC variable = disp_y value = 0.0 boundary = bottom  [../]


[./bkall] type = NeumannBC variable = b_k value = 0.0 boundary = 'bottom left right top'  [../]
[./Jkall] type = NeumannBC variable = J_k value = 0.0 boundary = 'bottom left right top' [../]
[]


[Preconditioning]
#[./myprec] type = FDP full = true []
[]

[Executioner]
 type = Transient
 solve_type = NEWTON
 #solve_type = LINEAR
 start_time = 0.0
 end_time = 2.0
 dt = 0.1
 dtmin = 0.1
 #nl_rel_tol = 1e-1
 # petsc_options = '-ksp_view_pmat'   # Stampa a schermo la matrice Jacobiana
 line_search = none
 # petsc_options_iname = '-pc_type -pc_hypre_type '
 # petsc_options_value = ' hypre    boomeramg     '
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '


#[Quadrature] type = GRID order = ELEVENTH []

[]

[Outputs]
 file_base = prova
 exodus = true
[]


