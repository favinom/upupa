//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElastoPlasticMaterial.h"

//registerMooseObject("OwlApp", ElastoPlasticMaterial);

InputParameters
ElastoPlasticMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");
  return params;
}

ElastoPlasticMaterial::ElastoPlasticMaterial(const InputParameters & parameters) : 
Material(parameters),
_dim(_mesh.dimension()),
_grad_disp_x(coupledGradient("disp_x")),
_grad_disp_y(coupledGradient("disp_y")),
_grad_disp_z(_dim>2 ? coupledGradient("disp_z") : _grad_zero),
_F(declareProperty<RealTensorValue>("deformationgradient")),
_J(declareProperty<Real>("deformationdeterminant")),
//_F_old(getMaterialPropertyOld<RealTensorValue>("deformationgradient")),
_P(declareProperty<RealTensorValue>("firstpiolakirchhoff")),
_Psi(declareProperty<Real>("referenceenergy")),
//_P_old(getMaterialPropertyOld<RealTensorValue>("firstpiolakirchhoff")),
//_sigma(declareProperty<RealTensorValue>("cauchystress")),
//_SIGMA(declareProperty<RealTensorValue>("mandelstress")),
_B_Kmat(getOptionalMaterialProperty<RealTensorValue>("microstructureBK")),
_elastoplasticMaterial(declareProperty<ElastoPlasticMaterial *>("elastoplasticmaterial"))
{
for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
	_I(i,j) = (i==j);
	}
}
}


void ElastoPlasticMaterial::initQpStatefulProperties(){
	if(_B_Kmat){
		_B_K = _B_Kmat[_qp];
	} 
	else{
		_B_K = _I;
	}	
	
	RealTensorValue U;
		for  (int j=0; j<3; ++j){
		U(0,j) = _grad_disp_x[_qp](j);
		U(1,j) = _grad_disp_y[_qp](j);
		U(2,j) = _grad_disp_z[_qp](j);
	}

	_F[_qp] = _I+U;
	_J[_qp] = _F[_qp].det();
	_P[_qp] = computeP(_F[_qp],_B_K);
	_Psi[_qp] = computePsi(_F[_qp],_B_K);

}

void ElastoPlasticMaterial::computeQpProperties()
{
	_elastoplasticMaterial[_qp] = this;
	
	if(_B_Kmat){
		_B_K = _B_Kmat[_qp];
	} 
	else{
		_B_K = _I;
	}	
	
	RealTensorValue U;
	
	for  (int j=0; j<3; ++j){
		U(0,j) = _grad_disp_x[_qp](j);
		U(1,j) = _grad_disp_y[_qp](j);
		U(2,j) = _grad_disp_z[_qp](j);
	}

	_F[_qp] = _I+U;
	_J[_qp] = _F[_qp].det();
	_P[_qp] = computeP(_F[_qp],_B_K);
	_Psi[_qp] = computePsi(_F[_qp],_B_K);
	
//	_sigma[_qp] = (_P[_qp]*_F[_qp].transpose())/_F[_qp].det();
//	_SIGMA[_qp] = _F[_qp].transpose()*_P[_qp];
}
