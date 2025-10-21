//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// keyword per evitare loop negli include

#include "Kernel.h"
//#include "PoroElastoGrowth.h"
//#include "ElastoPlasticGrowingMaterial.h"

class GrowthProliferantCells_v2 : public Kernel
{
public:
  static InputParameters validParams(); //funzione che legge parametri del Kernel

  GrowthProliferantCells_v2(const InputParameters & parameters); // costruttore della classe

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  	int const _dim;
  	Real const _Phisv;
	Real const _rhoR;
  	Real const _gamma_upt;
  	Real const _gamma_nec;
  	Real const _omega_cr;
  	Real const _omega_env;
  	Real const _delta1;
  	Real const _delta2;
  	
    	int const _id_x;
    	int const _id_y;
    	int const _id_z;
	int const _id_JK;
	int const _id_cp;
	int const _id_omega;
  	
    	MaterialProperty<Real> const & _J;
    	MaterialProperty<Real> const & _J_old;
	MaterialProperty<RealTensorValue> const & _F;
	MaterialProperty<RealTensorValue> const & _F_old;
  	MaterialProperty<RealTensorValue> const & _P_old;
  	
  	VariableValue const & _J_K;
  	
  	VariableValue const & _cp;
  	
  	VariableValue const & _omega;
};
