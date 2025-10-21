//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NutrientsAbsorbtion.h"

registerMooseObject("upupaApp", NutrientsAbsorbtion);

InputParameters
NutrientsAbsorbtion::validParams()
{
  InputParameters params = Kernel::validParams();
	params.addParam<Real>("solid_fraction", 0.5, "Solid volumetric fraction in the natural state");
	params.addParam<Real>("gamma_n1", 0.5, "Modulation of nutrient absorbtion 1");
  	params.addParam<Real>("gamma_n2", 0.5, "Modulation of nutrient absorbtion 2");
  	params.addParam<Real>("mass_density", 1000.0, "Mass density");
  	params.addRequiredCoupledVar("J_K", "determinante di K, parametro d'ordine per questo modello");
  	params.addRequiredCoupledVar("cp", "determinante di K, parametro d'ordine per questo modello");
  	params.addRequiredCoupledVar("omega", "determinante di K, parametro d'ordine per questo modello");
  return params;
}

NutrientsAbsorbtion::NutrientsAbsorbtion(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_Phisv(getParam<Real>("solid_fraction")),
_gamma_n1(getParam<Real>("gamma_n1")),
_gamma_n2(getParam<Real>("gamma_n2")),
_rhoR(getParam<Real>("mass_density")),
_id_JK(coupled("J_K")),
_id_cp(coupled("cp")),
_id_omega(coupled("omega")),
_J_K(coupledValue("J_K")),
_cp(coupledValue("cp")),
_omega(coupledValue("omega"))
{
}



Real NutrientsAbsorbtion::computeQpResidual(){

	return _Phisv*_gamma_n1/_rhoR*_omega[_qp]/(_omega[_qp] + _gamma_n2)*_cp[_qp]*_J_K[_qp]*_test[_i][_qp];
}




Real NutrientsAbsorbtion::computeQpJacobian(){

	return _Phisv*_gamma_n1/_rhoR*_gamma_n2/((_omega[_qp] + _gamma_n2)*(_omega[_qp] + _gamma_n2))*_cp[_qp]*_J_K[_qp]* _test[_i][_qp]*_phi[_j][_qp];
}



Real NutrientsAbsorbtion::computeQpOffDiagJacobian(unsigned int jvar){

if(jvar == _id_JK){
	return _Phisv*_gamma_n1/_rhoR*_omega[_qp]/(_omega[_qp] + _gamma_n2)*_cp[_qp]*_phi[_j][_qp]*_test[_i][_qp];
}
else if(jvar == _id_cp){
	return _Phisv*_gamma_n1/_rhoR*_omega[_qp]/(_omega[_qp] + _gamma_n2)*_J_K[_qp]*_phi[_j][_qp]*_test[_i][_qp];;
} 

return 0.0;

}


