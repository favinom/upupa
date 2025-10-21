//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CouplingGrowthBC.h"

registerMooseObject("upupaApp", CouplingGrowthBC);

InputParameters
CouplingGrowthBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
        params.addRequiredCoupledVar("disp_x", "Displacement along x");
  	params.addRequiredCoupledVar("disp_y", "Displacement along y");
  	params.addCoupledVar("disp_z","Displacement along z");
  	params.addRequiredCoupledVar("J_K", "determinante di K, parametro d'ordine per questo modello");
  return params;
}

CouplingGrowthBC::CouplingGrowthBC(const InputParameters & parameters) :
  IntegratedBC(parameters), 
_dim(_mesh.dimension()),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_JK(coupled("J_K")),
_disp_x(coupledValue("disp_x")),
_disp_y(coupledValue("disp_y")),
_disp_z(_dim>2 ? coupledValue("disp_z") : _zero),
_J_K(getMaterialProperty<Real>("growthdeterminant")),
_J(getMaterialProperty<Real>("deformationdeterminant")),
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
_P(getMaterialProperty<RealTensorValue>("firstpiolakirchhoff")),
_dPdJK(getMaterialProperty<RealTensorValue>("derivativefirstpiolajk")),
_elastoplasticgrowingMaterial(getMaterialProperty<ElastoPlasticGrowingMaterial *>("elastoplasticgrowingmaterial"))
{
}

Real CouplingGrowthBC::computeQpResidual() {
RealVectorValue disp;
disp(0) = _disp_x[_qp]; disp(1) = _disp_y[_qp]; disp(2) = _disp_z[_qp];

  return -(_P[_qp]*_normals[_qp])*disp*_test[_i][_qp]/(3*_J_K[_qp]);
}

Real CouplingGrowthBC::computeQpJacobian() {
RealVectorValue disp;
disp(0) = _disp_x[_qp]; disp(1) = _disp_y[_qp]; disp(2) = _disp_z[_qp];
  return ((_P[_qp] - _J_K[_qp]*_dPdJK[_qp])*_normals[_qp])*disp*_test[_i][_qp]*_phi[_j][_qp]/std::pow(_J_K[_qp],2);
}


Real CouplingGrowthBC::computeQpOffDiagJacobian(unsigned int jvar) {
RealVectorValue disp;
disp(0) = _disp_x[_qp]; disp(1) = _disp_y[_qp]; disp(2) = _disp_z[_qp];
RealVectorValue h;
RealTensorValue H;
RealTensorValue Jac;
int component_H;

for(int i=0;i<3;++i){
	h(i) = 0;
}

if(jvar == _id_JK){
std::cout << 101 << std::endl;
	return 0.0;
}
else if(jvar==_id_x){
	component_H = 0;
}
else if(jvar==_id_y){
	component_H = 1;
}
else if(jvar==_id_z){
	component_H = 2;
}
else{
	return 0.0;
}
h(component_H) = _phi[_j][_qp];
H(component_H,0) = _grad_phi[_j][_qp](0); H(component_H,1) = _grad_phi[_j][_qp](1); H(component_H,2) = _grad_phi[_j][_qp](2);

Jac = _elastoplasticgrowingMaterial[_qp][0].evaluateJac(H,_qp);

  return -(Jac*_normals[_qp])*disp*_test[_i][_qp]/_J_K[_qp] - (_P[_qp]*_normals[_qp])*h*_test[_i][_qp]/_J_K[_qp];
}

