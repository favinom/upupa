//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NeoHookeanGrowingRemodeling.h"

registerMooseObject("upupaApp", NeoHookeanGrowingRemodeling);

InputParameters
NeoHookeanGrowingRemodeling::validParams()
{
  InputParameters params = ElastoPlasticGrowingMaterial::validParams();
  params.addRequiredParam<Real>("lambda", "First Lamé constant");
  params.addRequiredParam<Real>("mu", "Second Lamé constant");
  params.addParam<Real>("solid_fraction", 1.0, "Solid volumetric fraction in the natural state");
  return params;
}

NeoHookeanGrowingRemodeling::NeoHookeanGrowingRemodeling(const InputParameters & parameters) : 
ElastoPlasticGrowingMaterial(parameters),
_lambda(getParam<Real>("lambda")*getParam<Real>("solid_fraction")),
_mu(getParam<Real>("mu")*getParam<Real>("solid_fraction")),
_Phisnu(getParam<Real>("solid_fraction"))
{
	//_lambda = _lambda*_Phisnu;
	//_mu = _mu*_Phisnu;
}



RealTensorValue NeoHookeanGrowingRemodeling::computeP(RealTensorValue const & F, RealTensorValue const & B_K, Real const & J_K){
// Quantità dipendenti da F
RealTensorValue FinvT;
Real J;

J = F.det();

FinvT = F.inverse().transpose();

//std::cout << B_K << std::endl;
	return _mu*std::pow(J_K,1.0/3.0)*F*B_K + J_K*(_lambda*log(J/J_K) - _mu)*FinvT;

}


Real NeoHookeanGrowingRemodeling::computePsi(RealTensorValue const & F, RealTensorValue const & B_K, Real const & J_K){
// Quantità dipendenti da F
RealTensorValue C;
Real J;

C = F.transpose()*F;
J = F.det();

return _mu/2.0*(  ( std::pow(J_K,-2.0/3.0) * B_K * C ).tr()  - 3.0 ) - _mu*log(J/J_K) + _lambda/2.0*log(J/J_K)*log(J/J_K);
}




RealTensorValue NeoHookeanGrowingRemodeling::evaluateJac(RealTensorValue const & H, int const & qp){
if(_B_Kmat){
	_B_K = _B_Kmat[qp];
} 
else{
	_B_K = _I;
}
if(_J_Kmat){
	_J_K = _J_Kmat[qp];
} 
else{
	_J_K = 1;
}

// Quantità dipendenti da F
RealTensorValue Finv;
RealTensorValue FinvT;
Real J;

Finv = _F[qp].inverse();

FinvT = Finv.transpose();

J = _F[qp].det();
	
	return _mu*std::pow(_J_K,1.0/3.0)*H*_B_K + _lambda*_J_K*(FinvT.contract(H))*FinvT - _J_K*(_lambda*log(J/_J_K) - _mu )*((Finv*H*Finv).transpose());
}



RealTensorValue NeoHookeanGrowingRemodeling::computedJK(RealTensorValue const & F, RealTensorValue const & B_K, Real const & J_K){
RealTensorValue FinvT;
Real J;

FinvT = (F.inverse()).transpose();

J = F.det();

	return 1.0/3.0*_mu*F*B_K/(std::pow(J_K,2.0/3.0)) + (_lambda*log(J/J_K)-_mu-_lambda )*FinvT;
}




