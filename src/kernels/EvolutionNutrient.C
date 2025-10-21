//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full resictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "EvolutionNutrient.h"

registerMooseObject("upupaApp", EvolutionNutrient);

InputParameters
EvolutionNutrient::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addParam<Real>("permeability", 1.0, "Permeability");
  params.addParam<Real>("m", 0.1, "motility coefficient");
  params.addParam<Real>("solid_fraction", 0.5, "Solid volumetric fraction in the natural state");
  params.addParam<Real>("diffusivity", 0.5, "Diffusivity of the nutrient");
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");
  params.addRequiredCoupledVar("p", "pressure");
  params.addRequiredCoupledVar("J_K", "determinant of the growth tensor");
  params.addRequiredCoupledVar("mu_K", "parte energetica potenziale chimico");
  params.addRequiredCoupledVar("omega", "nutrient concentration");
  return params;
}

EvolutionNutrient::EvolutionNutrient(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_k(getParam<Real>("permeability")),
_m(getParam<Real>("m")),
_Phisv(getParam<Real>("solid_fraction")),
_d(getParam<Real>("diffusivity")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_p(coupled("p")),
_id_Jk(coupled("J_K")),
_id_muk(coupled("mu_K")),
_id_omega(coupled("omega")),
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
_J(getMaterialProperty<Real>("deformationdeterminant")),
_J_old(getMaterialPropertyOld<Real>("deformationdeterminant")),
_P(getMaterialProperty<RealTensorValue>("firstpiolakirchhoff")),
_grad_p(coupledGradient("p")),
_J_K(getMaterialProperty<Real>("growthdeterminant")),
_J_K_old(getMaterialPropertyOld<Real>("growthdeterminant")),
_mu_K(coupledValue("mu_K")),
_grad_muk(coupledGradient("mu_K")),
_omega(coupledValue("omega")),
_grad_omega(coupledGradient("omega")),
_omega_old(coupledValueOld("omega"))
 {
 }





Real EvolutionNutrient::computeQpResidual(){
RealTensorValue K;
RealTensorValue D;
RealTensorValue M;
RealTensorValue Kmat;
RealTensorValue Dmat;
RealTensorValue Mmat;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		M(i,j) = _m*(i==j);
		K(i,j) = _k*(i==j);
		D(i,j) = _d*(i==j);
	}
}

Kmat = (_J[_qp]-_Phisv*_J_K[_qp])*_F[_qp].inverse()*K*(_F[_qp].inverse()).transpose();
Dmat = (_J[_qp]-_Phisv*_J_K[_qp])*_F[_qp].inverse()*D*(_F[_qp].inverse()).transpose();
Mmat = _J[_qp]*_F[_qp].inverse()*M*(_F[_qp].inverse()).transpose();
//std::cout << _J[_qp]-_Phisv*_J_K[_qp] << std::endl;

Real temp = ( (_J[_qp] - _Phisv*_J_K[_qp])*_omega[_qp] - (_J_old[_qp]-_Phisv*_J_K_old[_qp])*_omega_old[_qp] )/_dt;

return temp*_test[_i][_qp] + ( _omega[_qp]*(Kmat*_grad_p[_qp] + _Phisv*Mmat*_grad_muk[_qp]) + Dmat*_grad_omega[_qp] )*_grad_test[_i][_qp];
}



Real EvolutionNutrient::computeQpJacobian()
{
RealTensorValue K;
RealTensorValue D;
RealTensorValue M;
RealTensorValue Kmat;
RealTensorValue Dmat;
RealTensorValue Mmat;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		M(i,j) = _m*(i==j);
		K(i,j) = _k*(i==j);
		D(i,j) = _d*(i==j);
	}
}
Kmat = (_J[_qp]-_Phisv*_J_K[_qp])*_F[_qp].inverse()*K*(_F[_qp].inverse()).transpose();
Dmat = (_J[_qp]-_Phisv*_J_K[_qp])*_F[_qp].inverse()*D*(_F[_qp].inverse()).transpose();
Mmat = _J[_qp]*_F[_qp].inverse()*M*(_F[_qp].inverse()).transpose();



//Real temp = ( (_J[_qp] - _Phisv*_J_K[_qp])*_omega[_qp] - (_J_old[_qp] - _Phisv*_J_K_old[_qp])*_omega_old[_qp] )/_dt;
//return temp*_test[_i][_qp] + ( _omega[_qp]*(Kmat*_grad_p[_qp]) + Dmat*_grad_omega[_qp]   )*_grad_test[_i][_qp] ;
return (_J[_qp] - _Phisv*_J_K[_qp])/_dt*_test[_i][_qp]*_phi[_j][_qp] + (_phi[_j][_qp]*(Kmat*_grad_p[_qp]+_Phisv*Mmat*_grad_muk[_qp]) + Dmat*_grad_phi[_j][_qp] )*_grad_test[_i][_qp];
}




Real EvolutionNutrient::computeQpOffDiagJacobian(unsigned int jvar) 
{
RealTensorValue Finv;
RealTensorValue FinvT;
Finv = _F[_qp].inverse();
FinvT = Finv.transpose();

RealTensorValue K;
RealTensorValue D;
RealTensorValue M;
RealTensorValue Kmat;
RealTensorValue Dmat;
RealTensorValue Mmat;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		M(i,j) = _m*(i==j);
		K(i,j) = _k*(i==j);
		D(i,j) = _d*(i==j);
	}
}
Kmat = (_J[_qp]-_Phisv*_J_K[_qp])*Finv*K*FinvT;
Dmat = (_J[_qp]-_Phisv*_J_K[_qp])*Finv*D*FinvT;
Mmat = _J[_qp]*_F[_qp].inverse()*M*(_F[_qp].inverse()).transpose();

RealTensorValue H;
int component_H;

if(jvar == _id_Jk){
	return -_Phisv*_test[_i][_qp]*_phi[_j][_qp]*_omega[_qp]/_dt -_phi[_j][_qp]*_Phisv*( _omega[_qp]*(Finv*K*FinvT *_grad_p[_qp]) + Finv*D*FinvT*_grad_omega[_qp] ) *_grad_test[_i][_qp];
} 
else if(jvar==_id_p){
	return _omega[_qp]*(Kmat*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
}
else if(jvar== _id_muk){
	return _omega[_qp]*(_Phisv*Mmat*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
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
RealTensorValue tempK = Finv*H*Kmat;
RealTensorValue tempD = Finv*H*Dmat;
RealTensorValue tempM = Finv*H*Mmat;

Real ret1 = _J[_qp]*FinvT.contract(H)*_omega[_qp]/_dt*_test[_i][_qp];
Real ret2 = _J[_qp]*(Finv.transpose()).contract(H)*(_omega[_qp]*(Finv*K*FinvT*_grad_p[_qp]+_Phisv*Finv*M*FinvT*_grad_muk[_qp]) + Finv*D*FinvT*_grad_omega[_qp])*_grad_test[_i][_qp];
return ret1 + ret2 - (_omega[_qp]*((tempK+tempK.transpose())*_grad_p[_qp]+_Phisv*(tempM+tempM.transpose())*_grad_muk[_qp])+ (tempD+tempD.transpose())*_grad_omega[_qp])*_grad_test[_i][_qp];
}


