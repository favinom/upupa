//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OffDiagonalTimeDerivative.h"

registerMooseObject("upupaApp", OffDiagonalTimeDerivative);


InputParameters
OffDiagonalTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addRequiredCoupledVar("J_K", "determinante di K, parametro d'ordine per questo modello");
  return params;
}

OffDiagonalTimeDerivative::OffDiagonalTimeDerivative(const InputParameters & parameters) : 
TimeKernel(parameters),
_id_JK(coupled("J_K")),
_J_K(coupledValue("J_K")),
_J_K_old(coupledValueOld("J_K"))
{}

Real OffDiagonalTimeDerivative::computeQpResidual()
{
  return  -(_J_K[_qp]-_J_K_old[_qp])/_dt*_test[_i][_qp];   //_grad_u[_qp] * _grad_test[_i][_qp]+ _u[_qp]*_test[_i][_qp];
}

Real OffDiagonalTimeDerivative::computeQpJacobian()
{
  return 0.0;
}


Real OffDiagonalTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar){
if(jvar == _id_JK){
	return -_phi[_j][_qp]*_test[_i][_qp]/_dt;
}

return 0.0;


}
