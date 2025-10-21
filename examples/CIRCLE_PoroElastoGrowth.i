[Mesh]
  [cyl2d_iga]
    type = FileMeshGenerator
    file = quarterCircle_mesh_refined.xda
  []
  #uniform_refine = 2
[]




[Variables]
[disp_x] order = FIRST []
[disp_y] order = FIRST []
[p] order = FIRST []
[mu_K] order = FIRST []
[J_K] order = FIRST []
[cp] order = FIRST []
[omega] order = FIRST []
[]

[AuxVariables]
 [./residuo_jk]
  order = FIRST
  family = LAGRANGE
 [../]
 
  [./residuo_y]
  order = FIRST
  family = LAGRANGE
 [../]
[]

[Materials]
[./materialelastic] type = NeoHookeanGrowingRemodeling lambda = 2 mu = 1 disp_x = disp_x disp_y = disp_y []
[./materialgrowth] type = GradeOneGrowth J_K = J_K []
[./materialremodeling] type = GradeZeroRemodelingExplicit lambda_p = 0.0 []
[]



[Kernels]

[./elasticityx] type = PoroElastoGrowth variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y p = p J_K = J_K [../]

[./elasticityy] type = PoroElastoGrowth variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y p = p J_K = J_K save_in = residuo_y [../]



[./pressure] type = PoroElastoGrowthMassBalance variable = p permeability = 0.001 disp_x = disp_x disp_y = disp_y p = p J_K = J_K [../]



[./massbalancediffusion] type = ChemicalFluxDivergence variable = mu_K m = 1.0 disp_x = disp_x disp_y = disp_y J_K = J_K mu_K = mu_K [../]

[./offdiagonaltimederivative] type = OffDiagonalTimeDerivative variable = mu_K J_K = J_K [../]

[./growth] type = GrowthSolid variable = mu_K disp_x = disp_x disp_y = disp_y J_K = J_K cp = cp omega = omega [../]



[./potential] type = ChemicalPotential variable = J_K alpha = 0.01 disp_x = disp_x disp_y = disp_y J_K = J_K mu_K = mu_K save_in = residuo_jk [../]

[./coupling] type = EshelbyGrowthCoupling variable = J_K disp_x = disp_x disp_y = disp_y J_K = J_K save_in = residuo_jk [../]



[.proliferating] type = EvolutionProliferatingCell variable = cp m = 1.0 disp_x = disp_x disp_y = disp_y J_K = J_K mu_K = mu_K cp = cp [../]

[./growthcp] type = GrowthProliferantCells variable = cp disp_x = disp_x disp_y = disp_y J_K = J_K cp = cp omega = omega [../]



[.nutrients] type = EvolutionNutrient variable = omega permeability = 0.001 diffusivity = 0.5 disp_x = disp_x disp_y = disp_y p = p J_K = J_K mu_K = mu_K omega = omega [../]

[.nutrientsink] type = NutrientsAbsorbtion variable = omega J_K = J_K cp = cp omega = omega [../]
[]



[ICs]
[./J_K_IC] type = FunctionIC  variable = 'J_K'  function = 1  [../]
[./cp_IC] type = FunctionIC  variable = 'cp'  function = 0.5  [../]
[./omega_IC] type = FunctionIC  variable = 'omega'  function = 0.7 [../]
[]



[BCs]
# DISPLACEMENT
[./bottomy] type = DirichletBC variable = disp_y value = 0.0 boundary = bottom [../]
[./leftx] type = DirichletBC variable = disp_x value = 0.0 boundary = left [../]
[./outerx] type = NeumannPressure variable = disp_x disp_x = disp_x disp_y = disp_y component = 0 pressure = -1.0 function = 0.0 boundary = outer [../]
[./outery] type = NeumannPressure variable = disp_y disp_x = disp_x disp_y = disp_y component = 1 pressure = -1.0 function = 0.0 boundary = outer [../]


# PRESSURE
[./bottomp] type = NeumannBC variable = p value = 0.0 boundary = bottom [../]
[./leftp] type = NeumannBC variable = p value = 0.0 boundary = left [../]
[./outerp] type = DirichletBC variable = p value = 0.0 boundary = outer [../]

# POTENZIALE CHIMICO / TRASFERIMENTI DI MASSA
[./bottommuK] type = NeumannBC variable = mu_K value = 0.0 boundary = bottom [../]
[./leftmuK] type = NeumannBC variable = mu_K value = 0.0 boundary = left [../]

# CRESCITA / TRASFERIMENTI DI POTENZA
[./bottomJK] type = NeumannBC variable = J_K value = 0.0 boundary = bottom [../]
[./leftJK] type = NeumannBC variable = J_K value = 0.0 boundary = left [../]

# CELLULE PROLIFERANTI
[./bottomcp] type = NeumannBC variable = cp value = 0.0 boundary = bottom [../]
[./leftcp] type = NeumannBC variable = cp value = 0.0 boundary = left [../]

# NUTRIENTI
[./bottomomega] type = NeumannBC variable = omega value = 0.0 boundary = bottom [../]
[./leftomega] type = NeumannBC variable = omega value = 0.0 boundary = left [../]
[./outeromega] type = DirichletBC variable = omega value = 0.7 boundary = outer [../]

[]


[Preconditioning]
#[./myprec] type = FDP full = true []
[]

[Executioner]
 type = Transient
 solve_type = NEWTON
 #solve_type = LINEAR
 start_time = 0.0
 end_time = 20.0
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


