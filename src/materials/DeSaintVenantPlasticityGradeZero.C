//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DeSaintVenantPlasticityGradeZero.h"

registerMooseObject("upupaApp", DeSaintVenantPlasticityGradeZero);

InputParameters
DeSaintVenantPlasticityGradeZero::validParams()
{
  InputParameters params = ElastoPlasticMaterial::validParams();
  params.addRequiredParam<Real>("mu", "First Lamé constant");
  params.addRequiredParam<Real>("lambda", "Second Lamé constant");
  return params;
}

DeSaintVenantPlasticityGradeZero::DeSaintVenantPlasticityGradeZero(const InputParameters & parameters) : 
ElastoPlasticMaterial(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda"))
{
}



RealTensorValue DeSaintVenantPlasticityGradeZero::computeP(RealTensorValue const & F, RealTensorValue const & B_K){
// Quantità dipendenti da B_K
Real J_K;

J_K = 1.0/std::sqrt(B_K.det());

// Quantità dipendenti da F
RealTensorValue C;

C = F.transpose()*F;

return J_K*F*(2*_mu*(B_K*C*B_K - B_K) +_lambda*((B_K*C).tr()-3)*B_K);

}




Real DeSaintVenantPlasticityGradeZero::computePsi(RealTensorValue const & F, RealTensorValue const & B_K){
Real J_K;

J_K = 1.0/std::sqrt(B_K.det());

// Quantità dipendenti da F
RealTensorValue C;

C = F.transpose()*F;

return J_K*(_lambda/8.0*std::pow( (B_K*C).tr()-3, 2) +_mu/4.0*( (B_K*C*B_K*C).tr() - 2*(B_K*C).tr() -3 ) );

}


RealTensorValue DeSaintVenantPlasticityGradeZero::evaluateJac(RealTensorValue const & H, int const & qp){
if(_B_Kmat){
	_B_K = _B_Kmat[qp];
} 
else{
	_B_K = _I;
}
// Quantità dipendenti da B_K
Real J_K;

J_K = 1.0/std::sqrt(_B_K.det());

// Quantità dipendenti da F
RealTensorValue C;
RealTensorValue S;

C = _F[qp].transpose()*_F[qp];
S = J_K*(2*_mu*(_B_K*C*_B_K - _B_K) +_lambda*((_B_K*C).tr()-3)*_B_K);

RealTensorValue C_H;
RealTensorValue S_H;

C_H = _F[qp].transpose()*H + H.transpose()*_F[qp];

S_H = 2*J_K*(_mu*_B_K*C_H*_B_K+_lambda/2.0*(C_H*_B_K).tr()*_B_K);

	return H*S + _F[qp]*S_H;
}









