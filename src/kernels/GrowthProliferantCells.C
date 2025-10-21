//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GrowthProliferantCells.h"

registerMooseObject("upupaApp", GrowthProliferantCells);

InputParameters
GrowthProliferantCells::validParams()
{
  InputParameters params = Kernel::validParams();
  	params.addParam<Real>("solid_fraction", 0.5, "Solid volumetric fraction in the natural state");
  	params.addParam<Real>("mass_density", 1000.0, "Mass density");
  	params.addParam<Real>("gamma_upt", 0.5, "nutrients uptake");
  	params.addParam<Real>("gamma_nec", 0.5, "necrotic cells");
  	params.addParam<Real>("omega_cr", 0.2, "critial value of nutrients");
  	params.addParam<Real>("omega_env", 1.0, "reference environmental value of nutrients");
  	params.addParam<Real>("delta1", 0.5, "mechanotransduction constant 1");
  	params.addParam<Real>("delta2", 0.5, "mechanotransduction constant 2");
        params.addRequiredCoupledVar("disp_x", "Displacement along x");
  	params.addRequiredCoupledVar("disp_y", "Displacement along y");
  	params.addCoupledVar("disp_z","Displacement along z");
  	params.addRequiredCoupledVar("J_K", "determinante di K, parametro d'ordine per questo modello");
  	params.addRequiredCoupledVar("cp", "determinante di K, parametro d'ordine per questo modello");
  	params.addRequiredCoupledVar("omega", "determinante di K, parametro d'ordine per questo modello");
  return params;
}


GrowthProliferantCells::GrowthProliferantCells(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_Phisv(getParam<Real>("solid_fraction")),
_rhoR(getParam<Real>("mass_density")),
_gamma_upt(getParam<Real>("gamma_upt")),
_gamma_nec(getParam<Real>("gamma_nec")),
_omega_cr(getParam<Real>("omega_cr")),
_omega_env(getParam<Real>("omega_env")),
_delta1(getParam<Real>("delta1")),
_delta2(getParam<Real>("delta2")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_JK(coupled("J_K")),
_id_cp(coupled("cp")),
_id_omega(coupled("omega")),
_J_old(getMaterialPropertyOld<Real>("deformationdeterminant")),
_F_old(getMaterialPropertyOld<RealTensorValue>("deformationgradient")),
_P_old(getMaterialPropertyOld<RealTensorValue>("firstpiolakirchhoff")),
_J_K(coupledValue("J_K")),
_cp(coupledValue("cp")),
_omega(coupledValue("omega"))
{
}


Real GrowthProliferantCells::computeQpResidual(){

Real traceStress = -_P_old[_qp].contract(_F_old[_qp])/_J_old[_qp];
Real posStress = ( traceStress > 0.0);
Real posNutrient = (_omega[_qp]-_omega_cr > 0.0);
Real effectStress = 1 - _delta1*posStress*traceStress/(3*_delta2 + traceStress);
Real negNutrient = 1.0 - posNutrient;

    return _J_K[_qp]*_cp[_qp]*( _gamma_upt/_rhoR*posNutrient*(_omega[_qp]-_omega_cr)/(_omega_env-_omega_cr)*effectStress*(_J_old[_qp] - _J_K[_qp]*_Phisv)/(_J_old[_qp]*(1-_Phisv)) - _gamma_nec/_rhoR*negNutrient*(1.0-_omega[_qp]/_omega_cr) )*_test[_i][_qp];
}



Real GrowthProliferantCells::computeQpJacobian(){

Real traceStress = -_P_old[_qp].contract(_F_old[_qp])/_J_old[_qp];
Real posStress = ( traceStress > 0.0);
Real posNutrient = (_omega[_qp]-_omega_cr > 0.0);
Real effectStress = 1 - _delta1*posStress*traceStress/(3*_delta2 + traceStress);
Real negNutrient = 1.0 - posNutrient;


    return _J_K[_qp]*( _gamma_upt/_rhoR*posNutrient*(_omega[_qp]-_omega_cr)/(_omega_env-_omega_cr)*effectStress*(_J_old[_qp] - _J_K[_qp]*_Phisv)/(_J_old[_qp]*(1-_Phisv)) - _gamma_nec/_rhoR*negNutrient*(1.0-_omega[_qp]/_omega_cr) )*_test[_i][_qp]*_phi[_j][_qp];
}




Real GrowthProliferantCells::computeQpOffDiagJacobian(unsigned int jvar){
RealTensorValue H;
RealTensorValue Jac;
int component_H;

if(jvar == _id_JK){
	Real traceStress = -_P_old[_qp].contract(_F_old[_qp])/_J_old[_qp];
	Real posStress = ( traceStress > 0.0);
	Real posNutrient = (_omega[_qp]-_omega_cr > 0.0);
	Real effectStress = 1 - _delta1*posStress*traceStress/(3*_delta2 + traceStress);
	Real negNutrient = 1.0 - posNutrient;
	

    	return ( _cp[_qp]*( _gamma_upt/_rhoR*posNutrient*(_omega[_qp]-_omega_cr)/(_omega_env-_omega_cr)*effectStress*(_J_old[_qp] - _J_K[_qp]*_Phisv)/(_J_old[_qp]*(1-_Phisv)) - _gamma_nec/_rhoR*negNutrient*(1.0-_omega[_qp]/_omega_cr) ) -_J_K[_qp]*_cp[_qp]*( _gamma_upt/_rhoR*posNutrient*(_omega[_qp]-_omega_cr)/(_omega_env-_omega_cr)*effectStress*_Phisv/(_J_old[_qp]*(1-_Phisv)) )   ) *_test[_i][_qp]*_phi[_j][_qp];
}
else if(jvar==_id_cp){
	return 0.0;
}
else if(jvar==_id_omega){
	Real traceStress = -_P_old[_qp].contract(_F_old[_qp])/_J_old[_qp];
	Real posStress = ( traceStress > 0.0);
	Real posNutrient = (_omega[_qp]-_omega_cr > 0.0);
	Real effectStress = 1 - _delta1*posStress*traceStress/(3*_delta2 + traceStress);
	Real negNutrient = 1.0 - posNutrient;
	
    	return _J_K[_qp]*_cp[_qp]*( _gamma_upt/_rhoR*posNutrient/(_omega_env-_omega_cr)*effectStress*(_J_old[_qp] - _J_K[_qp]*_Phisv)/(_J_old[_qp]*(1-_Phisv)) + _gamma_nec/_rhoR*negNutrient/_omega_cr)*_test[_i][_qp]*_phi[_j][_qp];
}
else if(jvar==_id_x){
	component_H = 0;
}
else if(jvar==_id_y){
	component_H = 1;
}
else if(jvar==_id_z){
	component_H = 2;
} 
else{
	return 0.0;
}

 return 0.0;

}


