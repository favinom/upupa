//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NeoHookeanPlasticityGradeZero.h"

registerMooseObject("upupaApp", NeoHookeanPlasticityGradeZero);

InputParameters
NeoHookeanPlasticityGradeZero::validParams()
{
  InputParameters params = ElastoPlasticMaterial::validParams();
  params.addRequiredParam<Real>("mu", "First Lamé constant");
  params.addRequiredParam<Real>("lambda", "Second Lamé constant");
  return params;
}

NeoHookeanPlasticityGradeZero::NeoHookeanPlasticityGradeZero(const InputParameters & parameters) : 
ElastoPlasticMaterial(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda"))
{
}



RealTensorValue NeoHookeanPlasticityGradeZero::computeP(RealTensorValue const & F, RealTensorValue const & B_K){
// Quantità dipendenti da B_K
Real J_K;

J_K = 1.0/std::sqrt(B_K.det());

// Quantità dipendenti da F
RealTensorValue FinvT;
Real J;

J = F.det();

FinvT = F.inverse().transpose();

return _mu*J_K*F*B_K + (_lambda/J_K*J*(J-J_K) - _mu*J_K)*FinvT;

}


Real NeoHookeanPlasticityGradeZero::computePsi(RealTensorValue const & F, RealTensorValue const & B_K){
// Quantità dipendenti da B_K
Real J_K;

J_K = 1.0/std::sqrt(B_K.det());

// Quantità dipendenti da F
RealTensorValue C;
Real J;

C = F.transpose()*F;
J = F.det();

	return _mu/2.0*( (B_K*C).tr() -3 -2*log(J/J_K)  ) +_lambda/2.0*std::pow(J/J_K-1,2);



}




RealTensorValue NeoHookeanPlasticityGradeZero::evaluateJac(RealTensorValue const & H, int const & qp){
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
RealTensorValue Finv;
RealTensorValue FinvT;
Real J;

Finv = _F[qp].inverse();

FinvT = Finv.transpose();

J = _F[qp].det();

	return _mu*J_K*H*_B_K + _lambda/J_K*(2*J-J_K)*J*(FinvT.contract(H))*FinvT - (_lambda/J_K*J*(J-J_K) - _mu*J_K)*((Finv*H*Finv).transpose());

	//return 2J_K*_mu*(H*B_K - 2/3*((F*B_K).contract(H))*FinvT + trCB_K/3*segnaposto) + (_lambda+2/3*_mu)/J_K*( J*(2*J-J_K)*(FinvT.contract(H))*FinvT - J*(J-J_K)*segnaposto);



	//return J_K*_mu*((-2*(F*B_K*H).tr()*FinvT + (C*B_K).tr()*segnaposto)/3.0 + H*B_K) +(_lambda+2.0/3.0*_mu)/J_K*(J*(2*J-J_K)*FinvT.contract(H)*FinvT - J*(J-J_K)*segnaposto);
	
	
}



