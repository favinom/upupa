//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PoroElastoGrowth.h"

registerMooseObject("upupaApp", PoroElastoGrowth);

InputParameters
PoroElastoGrowth::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<int>("component", "Component");
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");
  params.addRequiredCoupledVar("p", "pressure");
  params.addRequiredCoupledVar("J_K", "volumetric growth");
  //params.addCoupledVar("density",
  //                     "The name of the temperature variable used in the "
  //                     "ComputeThermalExpansionEigenstrain.  (Not required for "
  //                     "simulations without temperature coupling.)");
  return params;
}

PoroElastoGrowth::PoroElastoGrowth(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_component(getParam<int>("component")),
_p(coupledValue("p")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_p(coupled("p")),
_id_J_K(coupled("J_K")),
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
_J(getMaterialProperty<Real>("deformationdeterminant")),
_J_K(getMaterialProperty<Real>("growthdeterminant")),
_P(getMaterialProperty<RealTensorValue>("firstpiolakirchhoff")),
_dPdJK(getMaterialProperty<RealTensorValue>("derivativefirstpiolajk")),
_elastoplasticgrowingMaterial(getMaterialProperty<ElastoPlasticGrowingMaterial *>("elastoplasticgrowingmaterial"))
{}



Real PoroElastoGrowth::computeQpResidual()
{
RealTensorValue V;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		V(i,j) = 0;
	}
}

V(_component,0) = _grad_test[_i][_qp](0); V(_component,1) = _grad_test[_i][_qp](1); V(_component,2) = _grad_test[_i][_qp](2);

//std::cout << _J_K[_qp] << std::endl;
	return (_P[_qp]-_J[_qp]*_p[_qp]*(_F[_qp].inverse()).transpose()).contract(V);
}





Real PoroElastoGrowth::computeQpJacobian(){
RealTensorValue H;
RealTensorValue Jac;
RealTensorValue V;
RealTensorValue Finv;
RealTensorValue FinvT;


//U(0,0) = _grad_disp_x[_qp](0); U(0,1) = _grad_disp_x[_qp](1); U(0,2) = _grad_disp_x[_qp](2);
//U(1,0) = _grad_disp_y[_qp](0); U(1,1) = _grad_disp_y[_qp](1); U(1,2) = _grad_disp_y[_qp](2);
//U(2,0) = _grad_disp_z[_qp](0); U(2,1) = _grad_disp_z[_qp](1); U(2,2) = _grad_disp_z[_qp](2);

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		V(i,j) = 0;
	}
}

Finv = _F[_qp].inverse();
FinvT = Finv.transpose();


H(_component,0) = _grad_phi[_j][_qp](0); H(_component,1) = _grad_phi[_j][_qp](1); H(_component,2) = _grad_phi[_j][_qp](2);

Jac = _elastoplasticgrowingMaterial[_qp][0].evaluateJac(H,_qp);

V(_component,0) = _grad_test[_i][_qp](0); V(_component,1) = _grad_test[_i][_qp](1); V(_component,2) = _grad_test[_i][_qp](2);

  return (      Jac    -    _J[_qp]*_p[_qp]*(FinvT.contract(H))*FinvT     +     _J[_qp]*_p[_qp]*(Finv*H*Finv).transpose()       ).contract(V);
}




Real PoroElastoGrowth::computeQpOffDiagJacobian(unsigned int jvar) 
{
RealTensorValue H;
RealTensorValue Jac;
RealTensorValue V;
int component_H;

RealTensorValue Finv;
RealTensorValue FinvT;


for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		V(i,j) = 0;
	}
}

V(_component,0) = _grad_test[_i][_qp](0); V(_component,1) = _grad_test[_i][_qp](1); V(_component,2) = _grad_test[_i][_qp](2);

Finv = _F[_qp].inverse();
FinvT = Finv.transpose();

if(jvar == _id_p){
	return -_J[_qp]*_phi[_j][_qp]*(FinvT.contract(V));
}
else if(jvar == _id_J_K){
	return _dPdJK[_qp].contract(V)*_phi[_j][_qp];
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

H(component_H,0) = _grad_phi[_j][_qp](0); H(component_H,1) = _grad_phi[_j][_qp](1); H(component_H,2) = _grad_phi[_j][_qp](2);

Jac = _elastoplasticgrowingMaterial[_qp][0].evaluateJac(H,_qp);

  return (Jac-_J[_qp]*_p[_qp]*(FinvT.contract(H))*FinvT + _J[_qp]*_p[_qp]*(Finv*H*Finv).transpose()).contract(V);
}


