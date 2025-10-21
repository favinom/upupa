[Mesh]

[gmg]
 type = GeneratedMeshGenerator # DistributedRectilinearMeshGenerator # GeneratedMesh   
 dim = 2
 nx = 2
 ny = 2
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
 elem_type = QUAD9
[]

#[rename] type = RenameBoundaryGenerator input = gmg old_boundary = 'bottom' new_boundary = 'bottom2'
#[]

[]

[Variables]
[disp_x] order = SECOND []
[disp_y] order = SECOND []
[p] order = FIRST []
[J_K] order = FIRST []
[]



[Materials]
[./materialelastic] type = NeoHookeanGrowingRemodeling mu = 1 lambda = 2.0 disp_x = disp_x disp_y = disp_y []
[./materialgrowth] type = GradeOneGrowth lambda_p = 0.0 J_K = J_K[]
[./materialremodeling] type = GradeZeroRemodelingExplicit lambda_p = 0.0[]
[]

[Kernels]
[./elasticityx] type = PoroElastoGrowth variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y p = p J_K = J_K [../]
[./elasticityy] type = PoroElastoGrowth variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y p = p J_K = J_K [../]
[./pressure] type = PoroElastoGrowthMassBalance variable = p disp_x = disp_x disp_y = disp_y p = p J_K = J_K [../]

[./timederivative] type = StandardTimeDerivative variable = J_K [../]
[]

[ICs]
[./J_K_IC] type = FunctionIC  variable = 'J_K'  function = 1  [../]
[]

[BCs]
[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = FunctionNeumannBC variable = disp_y function = -0.3*t*(t<1.5) boundary = top  [../]

[./leftp] type = DirichletBC variable = p value = 0.0 boundary = left  [../]
[./rightp] type = DirichletBC variable = p value = 0.0 boundary = right  [../]

[./bottomx] type = DirichletBC variable = disp_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = DirichletBC variable = disp_y value = 0.0 boundary = bottom  [../]

[]


[Preconditioning]
#[./myprec] type = FDP full = true []
[]

[Executioner]
 type = Transient
 solve_type = NEWTON
 #solve_type = LINEAR
 start_time = 0.0
 end_time = 0.5
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


