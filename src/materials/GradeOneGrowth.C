//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GradeOneGrowth.h"

registerMooseObject("upupaApp", GradeOneGrowth);

InputParameters
GradeOneGrowth::validParams()
{
  InputParameters params = Microstructure::validParams();
  	params.addRequiredCoupledVar("J_K", "determinante di K, parametro d'ordine per questo modello");
  return params;
}

GradeOneGrowth::GradeOneGrowth(const InputParameters & parameters) : 
Microstructure(parameters),
_J_K_var(coupledValue("J_K")),
_J_K(declareProperty<Real>("growthdeterminant")),
_F_K(declareProperty<RealTensorValue>("growthtensor"))
{
}


void GradeOneGrowth::initQpStatefulProperties(){
_J_K[_qp] = _J_K_var[_qp];
_F_K[_qp] = std::pow(_J_K[_qp],1.0/3.0)*_I;
}

void GradeOneGrowth::computeQpProperties()
{
_J_K[_qp] = _J_K_var[_qp];
_F_K[_qp] = std::pow(_J_K[_qp],1.0/3.0)*_I;

}
