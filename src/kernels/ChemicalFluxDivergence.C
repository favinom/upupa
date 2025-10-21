//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemicalFluxDivergence.h"

registerMooseObject("upupaApp", ChemicalFluxDivergence);

InputParameters
ChemicalFluxDivergence::validParams()
{
  InputParameters params = Kernel::validParams();
    	params.addParam<Real>("m", 0.1, "motility coefficient");
	params.addParam<Real>("sign_chemkernel", 1.0, "Sign of the kernel");
    	params.addRequiredCoupledVar("disp_x", "Displacement along x");
  	params.addRequiredCoupledVar("disp_y", "Displacement along y");
  	params.addCoupledVar("disp_z","Displacement along z");
  	params.addRequiredCoupledVar("J_K", "determinante di K, parametro d'ordine per questo modello");
        params.addRequiredCoupledVar("mu_K", "parte energetica potenziale chimico");
  return params;
}

ChemicalFluxDivergence::ChemicalFluxDivergence(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_m(getParam<Real>("m")),
_sign(getParam<Real>("sign_chemkernel")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_JK(coupled("J_K")),
_id_muk(coupled("mu_K")),
_J_K(coupledValue("J_K")),
_mu_K(coupledValue("mu_K")),
_grad_muk(coupledGradient("mu_K")),
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
_J(getMaterialProperty<Real>("deformationdeterminant"))
{
}



Real ChemicalFluxDivergence::computeQpResidual(){
// Da discutere M = _J[_qp]/_Jk[_qp]*_F[_qp].inverse()*_M*(_F[_qp].inverse()).transpose();

RealTensorValue M;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		M(i,j) = _m*(i==j);
	}
}

M = _J[_qp]*_F[_qp].inverse()*M*(_F[_qp].inverse()).transpose();

//std::cout << _var.number() << std::endl; stampa il numero associato alla variabile passata

	return _sign*(M*_grad_muk[_qp])*_grad_test[_i][_qp];
}




Real ChemicalFluxDivergence::computeQpJacobian(){
if (_var.number() == _id_muk){

RealTensorValue M;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		M(i,j) = _m*(i==j);
	}
}

M = _J[_qp]*_F[_qp].inverse()*M*(_F[_qp].inverse()).transpose();

	return _sign*(M*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
}

  return 0.0;
}



Real ChemicalFluxDivergence::computeQpOffDiagJacobian(unsigned int jvar){
RealTensorValue M;
for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		M(i,j) = _m*(i==j);
	}
}
RealTensorValue Mmat;
Mmat = _J[_qp]*_F[_qp].inverse()*M*(_F[_qp].inverse()).transpose();
RealTensorValue H;
RealTensorValue Finv;
int component_H;

Finv = _F[_qp].inverse();

if(jvar == _id_JK){
	return 0.0;
} 
else if(jvar == _id_muk && jvar != _var.number()){
// std::cout << "Off diagonal" << std::endl;
	return _sign*(Mmat*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
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

RealTensorValue temp = Finv*H*Mmat;

	return _sign*((Finv.transpose()).contract(H)*(Mmat*_grad_muk[_qp])*_grad_test[_i][_qp]-((temp+temp.transpose())*_grad_muk[_qp])*_grad_test[_i][_qp]);

}


