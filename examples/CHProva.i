[Mesh]

#[gmg]
 type = GeneratedMesh  # DistributedRectilinearMeshGenerator # GeneratedMeshGenerator
 dim = 2
 nx = 20
 ny = 20
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
 elem_type = QUAD9
#[]
[]

[Variables]
[J_K][]
[b_k][]
[]



[Materials]
[]

[Kernels]
[./timeevolutionjk] type = StandardTimeDerivative variable = J_K [../]
[./pippo] type = DiffusiveTimeDerivative variable = J_K  m = 0.01 kv = 1 J_K = J_K [../]
[./divchempot] type = ChemicalPotentialDiffusion variable = J_K  m = 0.01  J_K = J_K  b_k = b_k [../]

#[./chempotdefinizione] type = ChemicalPotential  variable = b_k  d = 10  kv = 1  J_K = J_K  b_k = b_k [../]
[./chempot] type = ChemicalPotentialReaction  variable = b_k  d = 0.02  J_K = J_K  b_k = b_k [../]
[]

[ICs]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1-0.5*cos(pi*x)*cos(pi*y)  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 0.9+0.2*exp(-10*(x^2+y^2))  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1+0.3*(x>0.899)*(y>0.899)  [../]
[./J_k_ic] type = FunctionIC  variable = 'J_K'  function = 1+0.1*(x<0.099)*(y<0.099)  [../]
[]

[BCs]
[./bkall] type = NeumannBC variable = b_k value = 0.0 boundary = 'bottom left right top'  [../]
[./Jkall] type = NeumannBC variable = J_K value = 0.0 boundary = 'bottom left right top' [../]
[]


[Preconditioning]
#[./myprec] type = FDP full = true []
[]

[Executioner]
 type = Transient
 solve_type = NEWTON
 #solve_type = LINEAR
 start_time = 0.0
 end_time = 1
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


