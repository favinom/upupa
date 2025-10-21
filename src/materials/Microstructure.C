//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Microstructure.h"

registerMooseObject("upupaApp", Microstructure);

InputParameters
Microstructure::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<Real>("lambda_p", 0.0, "rate rimodellamento");
  params.addParam<Real>("sigma_y", 0.0, "sforzo di snervamento");
  return params;
}

Microstructure::Microstructure(const InputParameters & parameters) : 
Material(parameters),
_dim(_mesh.dimension()),
_lambda_p(getParam<Real>("lambda_p")),
_sigma_y(getParam<Real>("sigma_y")),

//_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
_F_old(getMaterialPropertyOld<RealTensorValue>("deformationgradient")),

//_P(getMaterialProperty<RealTensorValue>("firstpiolakirchhoff")),
_P_old(getMaterialPropertyOld<RealTensorValue>("firstpiolakirchhoff"))

//_microstructurematerial(declareProperty<Microstructure *>("microstructurematerial"))


{
for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
	_I(i,j) = (i==j);
	}
}
}

