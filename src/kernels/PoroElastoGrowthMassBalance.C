//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full resictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "PoroElastoGrowthMassBalance.h"

registerMooseObject("upupaApp", PoroElastoGrowthMassBalance);

InputParameters
PoroElastoGrowthMassBalance::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addParam<Real>("permeability", 1.0, "Permeability");
  params.addParam<Real>("solid_fraction", 0.5, "Solid volumetric fraction in the natural state");
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");
  params.addRequiredCoupledVar("p", "pressure");
  params.addRequiredCoupledVar("J_K", "determinant of the growth tensor");
  return params;
}

PoroElastoGrowthMassBalance::PoroElastoGrowthMassBalance(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_k(getParam<Real>("permeability")),
_Phisv(getParam<Real>("solid_fraction")),
_grad_p(coupledGradient("p")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_p(coupled("p")),
_id_Jk(coupled("J_K")),
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
_J(getMaterialProperty<Real>("deformationdeterminant")),
_J_old(getMaterialPropertyOld<Real>("deformationdeterminant")),
_J_K(getMaterialProperty<Real>("growthdeterminant")),
_P(getMaterialProperty<RealTensorValue>("firstpiolakirchhoff"))
//_stress(getMaterialProperty<RealTensorValue>("stress")),
//_elasticmaterial(getMaterialProperty<ElasticMaterial *>("elasticmaterial"))
 {
 //std::cout << _dim << std::endl;
 }





Real PoroElastoGrowthMassBalance::computeQpResidual(){
RealTensorValue K;


for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		K(i,j) = _k*(i==j);
	}
}

K = (_J[_qp]-_Phisv*_J_K[_qp])*_F[_qp].inverse()*K*(_F[_qp].inverse()).transpose();

//std::cout << _J[_qp]-_Phisv*_J_K[_qp] << std::endl;

return (K*_grad_p[_qp])*_grad_test[_i][_qp] + (_J[_qp]-_J_old[_qp])/_dt*_test[_i][_qp];
}



Real PoroElastoGrowthMassBalance::computeQpJacobian()
{
RealTensorValue K;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		K(i,j) = _k*(i==j);
	}
}
K = (_J[_qp]-_Phisv*_J_K[_qp])*_F[_qp].inverse()*K*(_F[_qp].inverse()).transpose();

return (K*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
}




Real PoroElastoGrowthMassBalance::computeQpOffDiagJacobian(unsigned int jvar) 
{
RealTensorValue H;
RealTensorValue Finv;
int component_H;
RealTensorValue K;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		K(i,j) = _k*(i==j);
	}
}

Finv = _F[_qp].inverse();

if(jvar == _id_Jk){
	return -(_Phisv*_phi[_j][_qp]*(Finv*K*(Finv).transpose())*_grad_p[_qp])*_grad_test[_i][_qp];
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

//RealTensorValue temp = H*Finv*K;
//Real ret1=-_test[_i][_qp]/_dt*_J[_qp]*((Finv).transpose()).contract(H);
//Real ret2=-_J[_qp]*(Finv.transpose()).contract(H)*((Finv*K*Finv.transpose())*_grad_p[_qp])*_grad_test[_i][_qp];
//return ret1+ret2+(_J[_qp]-_Phisv*_Jk[_qp])*(Finv*(temp+temp.transpose())*Finv.transpose()*_grad_p[_qp])*_grad_test[_i][_qp];

RealTensorValue Kperm = (_J[_qp]-_Phisv*_J_K[_qp])*_F[_qp].inverse()*K*(_F[_qp].inverse()).transpose();
RealTensorValue temp = Finv*H*Kperm;

Real ret1=_test[_i][_qp]/_dt*_J[_qp]*((Finv).transpose()).contract(H);
Real ret2=_J[_qp]*(Finv.transpose()).contract(H)*((Finv*K*Finv.transpose())*_grad_p[_qp])*_grad_test[_i][_qp];
return ret1+ret2-((temp+temp.transpose())*_grad_p[_qp])*_grad_test[_i][_qp];
}


