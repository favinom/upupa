//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GradeZeroRemodelingExplicit.h"

registerMooseObject("upupaApp", GradeZeroRemodelingExplicit);

InputParameters
GradeZeroRemodelingExplicit::validParams()
{
  InputParameters params = Microstructure::validParams();

  return params;
}

GradeZeroRemodelingExplicit::GradeZeroRemodelingExplicit(const InputParameters & parameters) : 
Microstructure(parameters),
_detBK(declareProperty<Real>("determinanteBK")),
_B_K(declareProperty<RealTensorValue>("microstructureBK")),
_B_K_old(getMaterialPropertyOld<RealTensorValue>("microstructureBK"))
{
}


void GradeZeroRemodelingExplicit::initQpStatefulProperties(){
_B_K[_qp] = _I;
_detBK[_qp] = 1;
}

void GradeZeroRemodelingExplicit::computeQpProperties()
{
RealTensorValue F;
RealTensorValue P;


// mettere if

F = _F_old[_qp];
P = _P_old[_qp];

_sigma = (P*F.transpose())/(F.det());
_SIGMA = (F.transpose())*P;

RealTensorValue devsigma;
RealTensorValue devSIGMA;

devsigma = _sigma - (_sigma.tr()/3.0)*_I;
devSIGMA = _SIGMA - (_SIGMA.tr()/3.0)*_I;



Real normdevsigma;
normdevsigma = std::sqrt(devsigma.contract(devsigma));

Real gamma_p;
if(normdevsigma > 0.000001){
	gamma_p = 2.0*_lambda_p/normdevsigma*(normdevsigma-std::sqrt(2.0/3.0)*_sigma_y + std::fabs(normdevsigma-std::sqrt(2.0/3.0)*_sigma_y ));
} 
else{
	gamma_p = 0.0;
}

Real detSigmadev = _SIGMA.det()*gamma_p;
Real eig = std::sqrt(-detSigmadev);
if ( normdevsigma > 0.000001 & _dt > 1.0/eig & -detSigmadev>0.00001){
	std::cout <<"eig="<< eig << std::endl;
	std::cout << "detSigmadev="<< detSigmadev << std::endl;
//	std::cout << eig << std::endl;
	exit(1);
}

RealTensorValue A;
A = _I + gamma_p*_dt*devSIGMA;
_B_K[_qp] = _B_K_old[_qp]*A.inverse();


Real toll;
Real J_K;
toll = 0.01;

//if(abs(J_K-1) > toll){std::cout << J_K << std::endl;}




}
