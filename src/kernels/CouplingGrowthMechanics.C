//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CouplingGrowthMechanics.h"

registerMooseObject("upupaApp", CouplingGrowthMechanics);

InputParameters
CouplingGrowthMechanics::validParams()
{
  InputParameters params = Kernel::validParams();
        params.addRequiredCoupledVar("disp_x", "Displacement along x");
  	params.addRequiredCoupledVar("disp_y", "Displacement along y");
  	params.addCoupledVar("disp_z","Displacement along z");
  	params.addCoupledVar("p","Pressure field");
  	params.addRequiredCoupledVar("J_K", "determinante di K, parametro d'ordine per questo modello");
  return params;
}

CouplingGrowthMechanics::CouplingGrowthMechanics(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_p(coupled("p")),
_id_JK(coupled("J_K")),
_disp_x(coupledValue("disp_x")),
_disp_y(coupledValue("disp_y")),
_disp_z(_dim>2 ? coupledValue("disp_z") : _zero),
_grad_p(coupledGradient("p")),
_J_K(getMaterialProperty<Real>("growthdeterminant")),
_grad_JK(coupledGradient("J_K")),
_J(getMaterialProperty<Real>("deformationdeterminant")),
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
_P(getMaterialProperty<RealTensorValue>("firstpiolakirchhoff")),
_dPdJK(getMaterialProperty<RealTensorValue>("derivativefirstpiolajk")),
_Psinu(getMaterialProperty<Real>("naturalenergy")),
_elastoplasticgrowingMaterial(getMaterialProperty<ElastoPlasticGrowingMaterial *>("elastoplasticgrowingmaterial"))
{}


Real CouplingGrowthMechanics::computeQpResidual(){
RealVectorValue disp;
disp(0) = _disp_x[_qp]; disp(1) = _disp_y[_qp]; disp(2) = _disp_z[_qp];

Real temp = _Psinu[_qp] - ( _P[_qp].tr() - _J[_qp]* (((_F[_qp].inverse()).transpose())*_grad_p[_qp])*disp  )/(3*_J_K[_qp]) ;
  return temp*_test[_i][_qp] +(_P[_qp]*(_grad_test[_i][_qp]/_J_K[_qp] - _grad_JK[_qp]*_test[_i][_qp]/std::pow(_J_K[_qp],2)))*disp/3;
}


Real CouplingGrowthMechanics::computeQpJacobian(){
RealVectorValue disp;
disp(0) = _disp_x[_qp]; disp(1) = _disp_y[_qp]; disp(2) = _disp_z[_qp];

RealVectorValue moment = _grad_test[_i][_qp]/_J_K[_qp] - _grad_JK[_qp]*_test[_i][_qp]/std::pow(_J_K[_qp],2);
RealVectorValue moment2 = -_grad_test[_i][_qp]*_phi[_j][_qp]/std::pow(_J_K[_qp],2) + 2*_grad_JK[_qp]*_phi[_j][_qp]*_test[_i][_qp]/std::pow(_J_K[_qp],3) - _test[_i][_qp]*_grad_phi[_j][_qp]/std::pow(_J_K[_qp],2);

Real temp1 = -_P[_qp].contract(_F[_qp])/(3*std::pow(_J_K[_qp],2))*_test[_i][_qp]*_phi[_j][_qp];
Real temp2 = -_test[_i][_qp]*_phi[_j][_qp]*_dPdJK[_qp].tr()/(3*_J_K[_qp]);
Real temp3 = (_P[_qp].tr() - _J[_qp]* (((_F[_qp].inverse()).transpose())*_grad_p[_qp])*disp)*_test[_i][_qp]*_phi[_j][_qp]/(3*std::pow(_J_K[_qp],2));
Real temp4 = _phi[_j][_qp]*(_dPdJK[_qp]*moment)*disp/3;
Real temp5 = (_P[_qp]*moment2)*disp/3;
  return temp1 + temp2 + temp3 + temp4 + temp5;
}




Real CouplingGrowthMechanics::computeQpOffDiagJacobian(unsigned int jvar){
RealVectorValue disp;
disp(0) = _disp_x[_qp]; disp(1) = _disp_y[_qp]; disp(2) = _disp_z[_qp];
RealTensorValue H;
RealTensorValue Jac;
int component_H;

if(jvar == _id_JK){
	return 0.0;
} 
else if(jvar == _id_p){
	return -(_J[_qp]* (((_F[_qp].inverse()).transpose())*_grad_phi[_j][_qp])*disp  )/(3*_J_K[_qp])*_test[_i][_qp];
} 
else if(jvar==_id_x){
	component_H = 0;
}
else if(jvar==_id_y){
	component_H = 1;
}
else if(jvar==_id_z){
	component_H = 2;
} else{
	return 0.0;
}

RealVectorValue h;
h(0) = 0.0; h(1) = 0.0; h(2) = 0.0; h(component_H) = _phi[_j][_qp];

RealTensorValue Finv;
RealTensorValue FinvT;
Finv = _F[_qp].inverse();
FinvT = Finv.transpose();

H(component_H,0) = _grad_phi[_j][_qp](0); H(component_H,1) = _grad_phi[_j][_qp](1); H(component_H,2) = _grad_phi[_j][_qp](2);

Jac = _elastoplasticgrowingMaterial[_qp][0].evaluateJac(H,_qp);


Real temp1 = _P[_qp].contract(H)*_test[_i][_qp]/_J_K[_qp];
Real temp2 = -(Jac.tr()/3)*_test[_i][_qp]/_J_K[_qp];
Real temp3 = (_J[_qp]*(FinvT.contract(H))*(FinvT*_grad_p[_qp])*disp - _J[_qp]*((Finv*H*Finv).transpose()*_grad_p[_qp])*disp + _J[_qp]* (((_F[_qp].inverse()).transpose())*_grad_p[_qp])*h)*_test[_i][_qp]/(3*_J_K[_qp]);

RealVectorValue moment = _grad_test[_i][_qp]/_J_K[_qp] - _grad_JK[_qp]*_test[_i][_qp]/std::pow(_J_K[_qp],2);

Real temp4 = _test[_i][_qp]*((Jac*moment)*disp + (_P[_qp]*moment)*h)/3;

  return temp1 + temp2 + temp3 + temp4;

}


