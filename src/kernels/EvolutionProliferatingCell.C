//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EvolutionProliferatingCell.h"

registerMooseObject("upupaApp", EvolutionProliferatingCell);

InputParameters
EvolutionProliferatingCell::validParams()
{
  InputParameters params = Kernel::validParams();
      	params.addParam<Real>("m", 0.1, "motility coefficient");
  	params.addRequiredCoupledVar("disp_x", "Displacement along x");
  	params.addRequiredCoupledVar("disp_y", "Displacement along y");
  	params.addCoupledVar("disp_z","Displacement along z");
  	params.addRequiredCoupledVar("J_K", "determinante di K, parametro d'ordine per questo modello");
        params.addRequiredCoupledVar("mu_K", "parte energetica potenziale chimico");
        params.addRequiredCoupledVar("cp", "cellule proliferanti");
  return params;
}

EvolutionProliferatingCell::EvolutionProliferatingCell(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_m(getParam<Real>("m")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_JK(coupled("J_K")),
_id_muk(coupled("mu_K")),
_id_cp(coupled("cp")),
_J_K(coupledValue("J_K")),
_J_K_old(coupledValueOld("J_K")),
_grad_muk(coupledGradient("mu_K")),
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
_J(getMaterialProperty<Real>("deformationdeterminant")),
_cp(coupledValue("cp")),
_cp_old(coupledValueOld("cp"))
{
}

Real EvolutionProliferatingCell::computeQpResidual(){

// Da discutere M = _J[_qp]/_Jk[_qp]*_F[_qp].inverse()*_M*(_F[_qp].inverse()).transpose();

RealTensorValue M;
RealTensorValue MMat;
for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		M(i,j) = _m*(i==j);
	}
}

MMat = _J[_qp]*_F[_qp].inverse()*M*(_F[_qp].inverse()).transpose();


  return -(_J_K[_qp]*_cp[_qp] - _J_K_old[_qp]*_cp_old[_qp] )/_dt*_test[_i][_qp] + _cp[_qp]* (MMat *_grad_muk[_qp]) * _grad_test[_i][_qp];
}



Real EvolutionProliferatingCell::computeQpJacobian(){

RealTensorValue M;
RealTensorValue MMat;
for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		M(i,j) = _m*(i==j);
	}
}

MMat = _J[_qp]*_F[_qp].inverse()*M*(_F[_qp].inverse()).transpose();

  return -_J_K[_qp]*_phi[_j][_qp]*_test[_i][_qp]/_dt  + _phi[_j][_qp]* (MMat *_grad_muk[_qp]) * _grad_test[_i][_qp];
}




Real EvolutionProliferatingCell::computeQpOffDiagJacobian(unsigned int jvar){

RealTensorValue M;
RealTensorValue MMat;
for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		M(i,j) = _m*(i==j);
	}
}

MMat = _J[_qp]*_F[_qp].inverse()*M*(_F[_qp].inverse()).transpose();

RealTensorValue H;
RealTensorValue Finv;
int component_H;

Finv = _F[_qp].inverse();

if(jvar == _id_muk){
//std::cout << "ci sono 2" << std::endl;
	return  _cp[_qp]* (MMat *_grad_phi[_j][_qp]) * _grad_test[_i][_qp];

} else if(jvar == _id_JK){

	return -_phi[_j][_qp]*_cp[_qp]*_test[_i][_qp]/_dt;

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

RealTensorValue temp = Finv*H*MMat;

  return _cp[_qp]*( (Finv.transpose()).contract(H)*(MMat*_grad_muk[_qp])*_grad_test[_i][_qp]-((temp+temp.transpose())*_grad_muk[_qp])*_grad_test[_i][_qp]  );


}


