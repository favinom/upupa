//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemicalPotential.h"

registerMooseObject("upupaApp", ChemicalPotential);

InputParameters
ChemicalPotential::validParams()
{
  InputParameters params = Kernel::validParams();
  	params.addParam<Real>("alpha", 1.0, "Diffusion interface parameter");
  	params.addParam<Real>("kv", 1.0, "Dissipative rate");
  	params.addRequiredCoupledVar("disp_x", "Displacement along x");
  	params.addRequiredCoupledVar("disp_y", "Displacement along y");
  	params.addCoupledVar("disp_z","Displacement along z");
  	params.addRequiredCoupledVar("J_K", "determinante di K, parametro d'ordine per questo modello");
        params.addRequiredCoupledVar("mu_K", "parte energetica potenziale chimico");
  return params;
}

ChemicalPotential::ChemicalPotential(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_alpha(getParam<Real>("alpha")),
_kv(getParam<Real>("kv")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_JK(coupled("J_K")),
_id_muk(coupled("mu_K")),
_J_K(coupledValue("J_K")),
_J_K_old(coupledValueOld("J_K")),
_grad_JK(coupledGradient("J_K")),
_mu_K(coupledValue("mu_K")),
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
_J(getMaterialProperty<Real>("deformationdeterminant"))
{
}

Real ChemicalPotential::computeQpResidual(){

// Da discutere M = _J[_qp]/_Jk[_qp]*_F[_qp].inverse()*_M*(_F[_qp].inverse()).transpose();

RealTensorValue Alpha;
RealTensorValue AlphaMat;
for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		Alpha(i,j) = _alpha*(i==j);
	}
}

AlphaMat = _J[_qp]*_F[_qp].inverse()*Alpha*(_F[_qp].inverse()).transpose();


  return _mu_K[_qp]*_test[_i][_qp] + _kv/_dt*(_J_K[_qp]-_J_K_old[_qp])*_test[_i][_qp] + (AlphaMat *_grad_JK[_qp]) * _grad_test[_i][_qp];
}



Real ChemicalPotential::computeQpJacobian(){
RealTensorValue Alpha;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		Alpha(i,j) = _alpha*(i==j);
	}
}
RealTensorValue AlphaMat;
AlphaMat = _J[_qp]*_F[_qp].inverse()*Alpha*(_F[_qp].inverse()).transpose();

  return _kv/(_dt)*_phi[_j][_qp]*_test[_i][_qp]  + (AlphaMat*_grad_phi[_j][_qp]) * _grad_test[_i][_qp] ;
}

//_phi[_j][_qp]*_test[_i][_qp]

Real ChemicalPotential::computeQpOffDiagJacobian(unsigned int jvar){

RealTensorValue Alpha;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		Alpha(i,j) = _alpha*(i==j);
	}
}
RealTensorValue AlphaMat;
AlphaMat = _J[_qp]*_F[_qp].inverse()*Alpha*(_F[_qp].inverse()).transpose();
RealTensorValue H;
RealTensorValue Finv;
int component_H;

Finv = _F[_qp].inverse();

if(jvar == _id_muk){
//std::cout << "ci sono 2" << std::endl;
	return _phi[_j][_qp]*_test[_i][_qp];

} else if(jvar == _id_JK){

	return 0.0;

} else if(jvar==_id_x){
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

RealTensorValue temp = Finv*H*AlphaMat;

  return (Finv.transpose()).contract(H)*(AlphaMat*_grad_JK[_qp])*_grad_test[_i][_qp]-((temp+temp.transpose())*_grad_JK[_qp])*_grad_test[_i][_qp];


}


