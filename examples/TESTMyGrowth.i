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
[J_K][]
[]



[Materials]
[./materiale] type = NeoHookeanGrowingRemodeling mu = 1 lambda = 2.0 disp_x = disp_x disp_y = disp_y []
[./materialg] type = GradeOneGrowth lambda_p = 1 J_K = J_K[]
[]

[Kernels]
[./nonlinearelasticityx] type = PoroElastoGrowth variable = disp_x component = 0 disp_x = disp_x disp_y = disp_y p = p J_K = J_K [../]
[./nonlinearelasticityy] type = PoroElastoGrowth variable = disp_y component = 1 disp_x = disp_x disp_y = disp_y p = p J_K = J_K [../]

[./pressure] type = PoroElastoPlasticityMassBalance_v2 variable = p permeability = 0.0001 solid_fraction = 0.6 disp_x = disp_x disp_y = disp_y p = p J_K = J_K [../]

[./growth] type = DiffusiveTimeDerivative variable = J_K  m = 1 kv = 10000 J_K = J_K[../]
[]


[ICs]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1-0.5*cos(pi*x)*cos(pi*y)  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_K'  function = 0.9+0.2*exp(-10*(x^2+y^2))  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1  [../]
#[./J_k_ic] type = FunctionIC  variable = 'J_k'  function = 1+0.3*(x>0.899)*(y>0.899)  [../]
[./J_k_ic] type = FunctionIC  variable = 'J_K'  function = 1+20*(x<0.199)*(y<0.199)  [../]
[]


[BCs]
[./topx] type = DirichletBC variable = disp_x value = 0.0 boundary = top  [../]
[./topy] type = DirichletBC variable = disp_y value = 0.0 boundary = top  [../]

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
 # solve_type = LINEAR
 start_time = 0.0
 end_time = 1.0
 dt = 0.1
 dtmin = 0.01
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


