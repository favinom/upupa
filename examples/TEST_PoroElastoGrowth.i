[Mesh]

[gmg]
 type = GeneratedMeshGenerator # DistributedRectilinearMeshGenerator # GeneratedMesh   
 dim = 2
 nx = 20
 ny = 20
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
  order = SECOND
  family = LAGRANGE
 [../]
[]

[Materials]
[./materialelastic] type = NeoHookeanGrowingRemodeling mu = 1 lambda = 2 disp_x = disp_x disp_y = disp_y []
[./materialgrowth] type = GradeOneGrowth J_K = J_K []
[./materialremodeling] type = GradeZeroRemodelingExplicit lambda_p = 1.0 []
[]

[Kernels]
[./elasticityx] type = PoroElastoGrowth variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y p = p J_K = J_K [../]
[./elasticityy] type = PoroElastoGrowth variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y p = p J_K = J_K save_in = residuo_y [../]

[./pressure] type = PoroElastoGrowthMassBalance variable = p permeability = 0.001 disp_x = disp_x disp_y = disp_y p = p J_K = J_K [../]
#[./additionalfluxp] type = ChemicalFluxDivergence variable = p m = 1.0 solid_fraction = 0.5 sign_chemkernel = -1 disp_x = disp_x disp_y = disp_y J_K = J_K mu_K = mu_K [../]

[./massbalancediffusion] type = ChemicalFluxDivergence variable = mu_K m = 1.0 disp_x = disp_x disp_y = disp_y J_K = J_K mu_K = mu_K [../]
[./offdiagonaltimederivative] type = OffDiagonalTimeDerivative variable = mu_K J_K = J_K [../]

[./potential] type = ChemicalPotential variable = J_K alpha = 0.01 disp_x = disp_x disp_y = disp_y J_K = J_K mu_K = mu_K save_in = residuo_jk [../]
[./coupling] type = EshelbyGrowthCoupling variable = J_K disp_x = disp_x disp_y = disp_y J_K = J_K save_in = residuo_jk [../]
[./growth] type = GrowthSolid variable = J_K disp_x = disp_x disp_y = disp_y J_K = J_K cp = cp omega = omega save_in = residuo_jk [../]

[.proliferating] type = EvolutionProliferatingCell variable = cp m = 1.0 disp_x = disp_x disp_y = disp_y J_K = J_K mu_K = mu_K cp = cp [../]
[./growthcp] type = GrowthProliferantCells variable = cp disp_x = disp_x disp_y = disp_y J_K = J_K cp = cp omega = omega [../]

[.nutrients] type = EvolutionNutrient variable = omega permeability = 0.001 diffusivity = 0.5 disp_x = disp_x disp_y = disp_y p = p J_K = J_K mu_K = mu_K omega = omega [../]
[.nutrientsink] type = NutrientsAbsorbtion variable = omega cp = cp omega = omega [../]
[]

[ICs]
#[./J_K_IC] type = FunctionIC  variable = 'J_K'  function = 1  [../]
#[./J_K_IC] type = FunctionIC  variable = 'J_K'  function = 1+0.1*(x<0.299)*(y<0.299)  [../]
[./J_K_IC] type = FunctionIC  variable = 'J_K'  function = 1+0.3*(y>0.299)*(y<0.601)  [../]
[./cp_IC] type = FunctionIC  variable = 'cp'  function = 0.4+25*x^2*(1-x)^2*y^2*(1-y)^2  [../]
#[./omega_IC] type = FunctionIC  variable = 'omega'  function = 0.3+0.3*x-25*x^2*(1-x)^2*y^2*(1-y)^2 [../]
[./omega_IC] type = FunctionIC  variable = 'omega'  function = 0.7 [../]
[]

[BCs]
[./topx] type = NeumannBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = FunctionNeumannBC variable = disp_y function = -0.0*t*(t<1.5) boundary = top  [../]
#[./topy] type = FunctionDirichletBC variable = disp_y function = -0.3*t*((1+x)*(2-x)-2) boundary = top  [../]
#[./topomega] type = FunctionDirichletBC variable = omega function = 0.3+0.3*x boundary = top [../]

[./leftp] type = DirichletBC variable = p value = 0.0 boundary = left  [../]
[./rightp] type = DirichletBC variable = p value = 0.0 boundary = right  [../]

#[./leftomega] type = DirichletBC variable = omega value = 0.3 boundary = left  [../]
#[./rightomega] type = DirichletBC variable = omega value = 0.6 boundary = right  [../]

[./bottomx] type = DirichletBC variable = disp_x value = 0.0 boundary = bottom  [../]
[./bottomy] type = DirichletBC variable = disp_y value = 0.0 boundary = bottom  [../]
#[./bottomomega] type = FunctionDirichletBC variable = omega function = 0.3+0.3*x boundary = bottom [../]

[./topomega] type = FunctionDirichletBC variable = omega function = 0.7 boundary = top [../]
[./bottomomega] type = FunctionDirichletBC variable = omega function = 0.7 boundary = bottom [../]
[./leftomega] type = DirichletBC variable = omega value = 0.7 boundary = left  [../]
[./rightomega] type = DirichletBC variable = omega value = 0.7 boundary = right  [../]
[]


[Preconditioning]
[./myprec] type = FDP full = true []
[]

[Executioner]
 type = Transient
 solve_type = NEWTON
 #solve_type = LINEAR
 start_time = 0.0
 end_time = 1.0
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


