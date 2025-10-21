[Mesh]
  [cyl2d_iga]
    type = FileMeshGenerator
    file = octSphere_mesh.xda
  []
  #uniform_refine = 2
[]

[GlobalParams]
lambda = 13330
mu = 19990
solid_fraction = 0.6
lambda_p = 0.0 #0.00003
sigma_y = 10
permeability = 0.001755
alpha = 1e-04
kv = 3.6
m = 1e-04
diffusivity = 11.52
gamma_upt = 1.9253e-08
gamma_nec = 5.4e-07
gamma_lys = 4.14e-08
omega_cr = 2e-06
omega_env = 7e-06
delta1 = 0.7138
delta2 = 1541
mass_density = 1e-06
gamma_n1 = 1.08e-09
gamma_n2 = 1.48e-07
[]

[Variables]
[disp_x] order = FIRST []
[disp_y] order = FIRST []
[disp_z] order = FIRST []
[p] order = FIRST []
[mu_K] order = FIRST []
[J_K] order = FIRST []
[cp] order = FIRST []
[omega] order = FIRST []
[]


[Materials]
[./materialelastic] type = NeoHookeanGrowingRemodeling disp_x = disp_x disp_y = disp_y disp_z = disp_z[]
[./materialgrowth] type = GradeOneGrowth J_K = J_K []
[./materialremodeling] type = GradeZeroRemodelingExplicit []
[]



[Kernels]

[./elasticityx] type = PoroElastoGrowth variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y disp_z = disp_z p = p J_K = J_K [../]
[./elasticityy] type = PoroElastoGrowth variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y disp_z = disp_z p = p J_K = J_K [../]
[./elasticityz] type = PoroElastoGrowth variable = disp_z component = 2 disp_x = disp_x disp_y = disp_y disp_z = disp_z p = p J_K = J_K [../]



[./pressure] type = PoroElastoGrowthMassBalance variable = p disp_x = disp_x disp_y = disp_y disp_z = disp_z p = p J_K = J_K [../]



[./massbalancediffusion] type = ChemicalFluxDivergence variable = mu_K disp_x = disp_x disp_y = disp_y disp_z = disp_z J_K = J_K mu_K = mu_K [../]
[./offdiagonaltimederivative] type = OffDiagonalTimeDerivative variable = mu_K J_K = J_K [../]
[./growth] type = GrowthSolid_v2 variable = mu_K disp_x = disp_x disp_y = disp_y disp_z = disp_z J_K = J_K cp = cp omega = omega [../]



[./potential] type = ChemicalPotential variable = J_K disp_x = disp_x disp_y = disp_y disp_z = disp_z J_K = J_K mu_K = mu_K  [../]
[./coupling] type = EshelbyGrowthCoupling variable = J_K disp_x = disp_x disp_y = disp_y disp_z = disp_z J_K = J_K  [../]



[.proliferating] type = EvolutionProliferatingCell variable = cp disp_x = disp_x disp_y = disp_y disp_z = disp_z J_K = J_K mu_K = mu_K cp = cp [../]
[./growthcp] type = GrowthProliferantCells_v2 variable = cp disp_x = disp_x disp_y = disp_y disp_z = disp_z J_K = J_K cp = cp omega = omega [../]



[.nutrients] type = EvolutionNutrient variable = omega disp_x = disp_x disp_y = disp_y disp_z = disp_z p = p J_K = J_K mu_K = mu_K omega = omega [../]
[.nutrientsink] type = NutrientsAbsorbtion variable = omega J_K = J_K cp = cp omega = omega [../]

[]



[ICs]
[./J_K_IC] type = FunctionIC  variable = 'J_K'  function = 1  [../]
[./cp_IC] type = FunctionIC  variable = 'cp'  function = 1  [../]
[./omega_IC] type = FunctionIC  variable = 'omega'  function = 7e-06 [../]
[]



[BCs]
# "back" 	-> piano XY
# "bottom" 	-> piano YZ
# "left" 	-> piano XZ
# "7" 		-> superficie

# DISPLACEMENT
[./bottomy] type = DirichletBC variable = disp_x value = 0.0 boundary = bottom [../]
[./leftx] type = DirichletBC variable = disp_y value = 0.0 boundary = left [../]
[./backz] type = DirichletBC variable = disp_z value = 0.0 boundary = back [../]
[./outerx] type = NeumannPressure variable = disp_x disp_x = disp_x disp_y = disp_y disp_z = disp_z component = 0 pressure = -1.0 function = 0.0 boundary = 7 [../]
[./outery] type = NeumannPressure variable = disp_y disp_x = disp_x disp_y = disp_y disp_z = disp_z component = 1 pressure = -1.0 function = 0.0 boundary = 7 [../]
[./outerz] type = NeumannPressure variable = disp_z disp_x = disp_x disp_y = disp_y disp_z = disp_z component = 2 pressure = -1.0 function = 0.0 boundary = 7 [../]


# PRESSURE
[./bottomp] type = NeumannBC variable = p value = 0.0 boundary = bottom [../]
[./leftp] type = NeumannBC variable = p value = 0.0 boundary = left [../]
[./backp] type = NeumannBC variable = p value = 0.0 boundary = back [../]
[./outerp] type = DirichletBC variable = p value = 0.0 boundary = 7 [../]

# POTENZIALE CHIMICO / TRASFERIMENTI DI MASSA
[./bottommuK] type = NeumannBC variable = mu_K value = 0.0 boundary = bottom [../]
[./leftmuK] type = NeumannBC variable = mu_K value = 0.0 boundary = left [../]
[./backmuK] type = NeumannBC variable = mu_K value = 0.0 boundary = back [../]

# CRESCITA / TRASFERIMENTI DI POTENZA
[./bottomJK] type = NeumannBC variable = J_K value = 0.0 boundary = bottom [../]
[./leftJK] type = NeumannBC variable = J_K value = 0.0 boundary = left [../]
[./backJK] type = NeumannBC variable = J_K value = 0.0 boundary = back [../]

# CELLULE PROLIFERANTI
[./bottomcp] type = NeumannBC variable = cp value = 0.0 boundary = bottom [../]
[./leftcp] type = NeumannBC variable = cp value = 0.0 boundary = left [../]
[./backcp] type = NeumannBC variable = cp value = 0.0 boundary = back [../]

# NUTRIENTI
[./bottomomega] type = NeumannBC variable = omega value = 0.0 boundary = bottom [../]
[./leftomega] type = NeumannBC variable = omega value = 0.0 boundary = left [../]
[./backomega] type = NeumannBC variable = omega value = 0.0 boundary = back [../]
[./outeromega] type = DirichletBC variable = omega value = 7e-06 boundary = 7 [../]

[]


[Preconditioning]
#[./myprec] type = FDP full = true []
[]

[Executioner]
 type = Transient
 solve_type = NEWTON
 #solve_type = LINEAR
 start_time = 0.0
 end_time = 480
 dt = 24
 dtmin = 0.1
 nl_rel_tol = 1e-04
 # petsc_options = '-ksp_view_pmat'   # Stampa a schermo la matrice Jacobiana
 line_search = none
 # petsc_options_iname = '-pc_type -pc_hypre_type '
 # petsc_options_value = ' hypre    boomeramg     '
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '

[Quadrature] type = TRAP []
#[Quadrature] type = GRID order = ELEVENTH []

[]

[Outputs]
 file_base = prova
 exodus = true
[]


