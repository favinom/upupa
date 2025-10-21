//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HolmesMowPlasticityGradeZero.h"

registerMooseObject("upupaApp", HolmesMowPlasticityGradeZero);

InputParameters
HolmesMowPlasticityGradeZero::validParams()
{
  InputParameters params = ElastoPlasticMaterial::validParams();
  params.addParam<Real>("mu", 1.0, "First Lamé constant");
  params.addParam<Real>("lambda", 1.0, "Second Lamé constant");
  return params;
}

HolmesMowPlasticityGradeZero::HolmesMowPlasticityGradeZero(const InputParameters & parameters) : 
ElastoPlasticMaterial(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda"))
{
}




RealTensorValue HolmesMowPlasticityGradeZero::computeP(RealTensorValue const & F, RealTensorValue const & B_K){
// Quantità dipendenti da B_K
Real J_K;

J_K = 1.0/std::sqrt(B_K.det());

// Quantità dipendenti da F
RealTensorValue FinvT;
RealTensorValue C;
RealTensorValue CB_K;
Real J;

J = F.det();
FinvT = F.inverse().transpose();
C = F.transpose()*F;
CB_K = C*B_K;


// Coefficienti del materiale
Real alpha0;
Real alpha1;
Real alpha2;
Real alpha3;

alpha0 = (2*_mu+_lambda)/4;
alpha1 = (2*_mu-_lambda)/(2*_mu+_lambda);
alpha2 = _lambda/(2*_mu+_lambda);
alpha3 = 1;


// Invarianti
Real I_1e;
Real I_2e;
Real I_3e;
Real Psi;

I_1e = CB_K.tr();
I_2e = 0.5*(I_1e*I_1e - (CB_K*CB_K).tr());
I_3e = (J*J)/(J_K*J_K);
Psi = exp(alpha1*(I_1e -3) + alpha2*(I_2e-3) - alpha3*log(I_3e));

	return 2*alpha0*J_K*Psi*( (alpha1+alpha2*I_1e)*F*B_K - alpha2*F*B_K*CB_K - alpha3*FinvT);

}



Real HolmesMowPlasticityGradeZero::computePsi(RealTensorValue const & F, RealTensorValue const & B_K){
// Quantità dipendenti da B_K
Real J_K;

J_K = 1.0/std::sqrt(B_K.det());

// Quantità dipendenti da F
RealTensorValue FinvT;
RealTensorValue C;
RealTensorValue CB_K;
Real J;

J = F.det();
FinvT = F.inverse().transpose();
C = F.transpose()*F;
CB_K = C*B_K;


// Coefficienti del materiale
Real alpha0;
Real alpha1;
Real alpha2;
Real alpha3;

alpha0 = (2*_mu+_lambda)/4;
alpha1 = (2*_mu-_lambda)/(2*_mu+_lambda);
alpha2 = _lambda/(2*_mu+_lambda);
alpha3 = 1;


// Invarianti
Real I_1e;
Real I_2e;
Real I_3e;

I_1e = CB_K.tr();
I_2e = 0.5*(I_1e*I_1e - (CB_K*CB_K).tr());
I_3e = (J*J)/(J_K*J_K);
	
	return alpha0*( exp(alpha1*(I_1e -3) + alpha2*(I_2e-3) - alpha3*log(I_3e))   -1);
}

/*
void HolmesMowPlasticityGradeZero::computeQpProperties()
{
if(_B_Kmat){
	_B_K = _B_Kmat[_qp];
} 
else{
	_B_K = _I;
}
// Quantità dipendenti da K
Real J_K;

J_K = std::sqrt(_B_K.det());

// Proprietà di F
RealTensorValue U;
RealTensorValue F;
RealTensorValue FinvT;
RealTensorValue C;
RealTensorValue CB_K;
Real trC;
Real J;

U(0,0) = _grad_disp_x[_qp](0); U(0,1) = _grad_disp_x[_qp](1); U(0,2) = _grad_disp_x[_qp](2);
U(1,0) = _grad_disp_y[_qp](0); U(1,1) = _grad_disp_y[_qp](1); U(1,2) = _grad_disp_y[_qp](2);
U(2,0) = _grad_disp_z[_qp](0); U(2,1) = _grad_disp_z[_qp](1); U(2,2) = _grad_disp_z[_qp](2);

F = _I+U;

J = F.det();

FinvT = F.inverse().transpose();

C = F.transpose()*F;

trC = C.tr();

CB_K = C*B_K;

// Coefficienti del materiale
Real alpha0;
Real alpha1;
Real alpha2;
Real alpha3;

alpha0 = (2*_mu+_lambda)/4;
alpha1 = (2*_mu-_lambda)/(2*_mu+_lambda);
alpha2 = _lambda/(2*_mu+_lambda);
alpha3 = 1;


// Invarianti
Real I_1e;
Real I_2e;
Real I_3e;
Real Psi;

I_1e = CB_K.tr();
I_2e = 0.5*(I_1e*I_1e - (CB_K*CB_K).tr());
I_3e = (J*J)/(J_K*J_K);

Psi = exp(alpha1*(I_1e -3) + alpha2*(I_2e-3) - alpha3*log(I_3e));


	_stress[_qp] = 2*alpha0*J_K*Psi*( (alpha1+alpha2*I_1e)*F*B_K - alpha2*F*B_K*CB_K - alpha3*FinvT);
	
	_sigma[_qp] = _stress[_qp]*(F.transpose())/J;
	
	_Sigma[_qp] = (F.transpose())*_stress[_qp];

// std::cout << (alpha1+alpha2*I_1e)*F*B_K - alpha2*F*B_K*CB_K;

}*/




RealTensorValue HolmesMowPlasticityGradeZero::evaluateJac(RealTensorValue const & H, int const & qp){
if(_B_Kmat){
	_B_K = _B_Kmat[qp];
} 
else{
	_B_K = _I;
}
// Quantità dipendenti da K
Real J_K;

J_K = 1.0/std::sqrt(_B_K.det());

// Proprietà di F
RealTensorValue U;
RealTensorValue F;
RealTensorValue FinvT;
RealTensorValue C;
RealTensorValue Cinv;
RealTensorValue CB_K;
RealTensorValue FB_K;
RealTensorValue D_FTH;
Real trC;
Real J;

U(0,0) = _grad_disp_x[qp](0); U(0,1) = _grad_disp_x[qp](1); U(0,2) = _grad_disp_x[qp](2);
U(1,0) = _grad_disp_y[qp](0); U(1,1) = _grad_disp_y[qp](1); U(1,2) = _grad_disp_y[qp](2);
U(2,0) = _grad_disp_z[qp](0); U(2,1) = _grad_disp_z[qp](1); U(2,2) = _grad_disp_z[qp](2);

F = _I+U;

J = F.det();

FinvT = (F.inverse()).transpose();

C = F.transpose()*F;

Cinv = C.inverse();

trC = C.tr();

CB_K = C*_B_K;

FB_K = F*_B_K;

D_FTH = ((F.transpose()*H) + H.transpose()*F)/2.0;

// Coefficienti del materiale
Real alpha0;
Real alpha1;
Real alpha2;
Real alpha3;

alpha0 = (2*_mu+_lambda)/4;
alpha1 = (2*_mu-_lambda)/(2*_mu+_lambda);
alpha2 = _lambda/(2*_mu+_lambda);
alpha3 = 1;


// Invarianti
Real I_1e;
Real I_2e;
Real I_3e;
Real ePsi;
Real Psi;

I_1e = CB_K.tr();
I_2e = 0.5*(I_1e*I_1e - (CB_K*CB_K).tr());
I_3e = (J*J)/(J_K*J_K);
Psi = alpha1*(I_1e -3) + alpha2*(I_2e-3) - alpha3*log(I_3e);
ePsi = exp(Psi);

// Secondo Piola-Kirchhoff
RealTensorValue S;
Real DerivataPsi;
S = 2*alpha0*J_K*ePsi*((alpha1+alpha2*I_1e)*_B_K - alpha2*_B_K*CB_K - alpha3*Cinv);
DerivataPsi = 2*((alpha1+alpha2*I_1e)*FB_K-alpha3*FinvT).contract(H)-2*alpha2*(CB_K*D_FTH).contract(_B_K);
//DerivataPsi = (2*alpha1+2*alpha2*I_1e)*FB_K.contract(H)-2*alpha2*(CB_K*D_FTH*_B_K).tr()-2*alpha3*FinvT.contract(H);
	
	return (H + DerivataPsi*F)*S + 2*alpha0*J_K*ePsi*( 2*alpha2*(FB_K.contract(H))*FB_K - 2*alpha2*FB_K*D_FTH*_B_K + alpha3*(H*F.inverse() + FinvT*H.transpose())*FinvT);


}
