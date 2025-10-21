//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EshelbyGrowthCoupling.h"

registerMooseObject("upupaApp", EshelbyGrowthCoupling);

InputParameters
EshelbyGrowthCoupling::validParams()
{
  InputParameters params = Kernel::validParams();
  	params.addParam<Real>("coupling", 1.0, "parametro per ridurre o aumentare il coupling");
        params.addRequiredCoupledVar("disp_x", "Displacement along x");
  	params.addRequiredCoupledVar("disp_y", "Displacement along y");
  	params.addCoupledVar("disp_z","Displacement along z");
  	params.addRequiredCoupledVar("J_K", "determinante di K, parametro d'ordine per questo modello");
  return params;
}

EshelbyGrowthCoupling::EshelbyGrowthCoupling(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_control(getParam<Real>("coupling")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_JK(coupled("J_K")),
_disp_x(coupledValue("disp_x")),
_disp_y(coupledValue("disp_y")),
_disp_z(_dim>2 ? coupledValue("disp_z") : _zero),
_J_K(coupledValue("J_K")),
_J(getMaterialProperty<Real>("deformationdeterminant")),
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
_P(getMaterialProperty<RealTensorValue>("firstpiolakirchhoff")),
_dPdJK(getMaterialProperty<RealTensorValue>("derivativefirstpiolajk")),
_Psinu(getMaterialProperty<Real>("naturalenergy")),
_elastoplasticgrowingMaterial(getMaterialProperty<ElastoPlasticGrowingMaterial *>("elastoplasticgrowingmaterial"))
{
}


Real EshelbyGrowthCoupling::computeQpResidual(){
//std::cout << _J_K[_qp] << _F[_qp] << _Psinu[_qp] << std::endl;
//  return _control*( _Psinu[_qp]    -    *_P[_qp].contract(_F[_qp])/(3*_J_K[_qp]) )*_test[_i][_qp];
    return  (_Psinu[_qp] -   _P[_qp].contract(_F[_qp]) / (  3*_J_K[_qp] )       )*_test[_i][_qp];
//  return _control*( -   _P[_qp].contract(_F[_qp])/(3*_J_K[_qp]) )*_test[_i][_qp];
}


Real EshelbyGrowthCoupling::computeQpJacobian(){

 return - _phi[_j][_qp]*_test[_i][_qp]*_dPdJK[_qp].contract(_F[_qp])/(3*_J_K[_qp]);

//  return -_phi[_j][_qp]*_test[_i][_qp]*_P[_qp].contract(_F[_qp])/(3.0*_J_K[_qp]*_J_K[_qp]);
  
//  return -_control*_phi[_j][_qp]*_test[_i][_qp]*_dPdJK[_qp].contract(_F[_qp])/(3*_J_K[_qp]) + _control*_phi[_j][_qp]*_test[_i][_qp]*_P[_qp].contract(_F[_qp])/(3*std::pow(_J_K[_qp],2));
}




Real EshelbyGrowthCoupling::computeQpOffDiagJacobian(unsigned int jvar){
RealTensorValue H;
RealTensorValue Jac;
int component_H;

if(jvar == _id_JK){
std::cout << jvar << std::endl;
//	return -_control*_phi[_j][_qp]*_test[_i][_qp]*_dPdJK[_qp].contract(_F[_qp])/(3*_J_K[_qp]);
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

H(component_H,0) = _grad_phi[_j][_qp](0); H(component_H,1) = _grad_phi[_j][_qp](1); H(component_H,2) = _grad_phi[_j][_qp](2);

Jac = _elastoplasticgrowingMaterial[_qp][0].evaluateJac(H,_qp);

 return (2*_P[_qp].contract(H) - Jac.contract(_F[_qp]))*_test[_i][_qp]/(3*_J_K[_qp]);
//  return _P[_qp].contract(H)*_test[_i][_qp]/(_J_K[_qp]);

}


