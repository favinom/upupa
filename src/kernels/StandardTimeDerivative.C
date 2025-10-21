//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StandardTimeDerivative.h"

registerMooseObject("upupaApp", StandardTimeDerivative);

InputParameters
StandardTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  return params;
}

StandardTimeDerivative::StandardTimeDerivative(const InputParameters & parameters) : 
TimeKernel(parameters)
{}

Real StandardTimeDerivative::computeQpResidual()
{
  return  -_u_dot[_qp]*_test[_i][_qp];   //_grad_u[_qp] * _grad_test[_i][_qp]+ _u[_qp]*_test[_i][_qp];
}

Real StandardTimeDerivative::computeQpJacobian()
{
  return -_phi[_j][_qp]* _test[_i][_qp]* _du_dot_du[_qp];
}


Real StandardTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar){

return 0.0;

}
