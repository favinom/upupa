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
#include "ElastoPlasticGrowingMaterial.h"

class EshelbyGrowthCoupling : public Kernel
{
public:
  static InputParameters validParams(); //funzione che legge parametri del Kernel

  EshelbyGrowthCoupling(const InputParameters & parameters); // costruttore della classe

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  	int const _dim;
  	Real const _control;
    	int const _id_x;
    	int const _id_y;
    	int const _id_z;
	int const _id_JK;

	VariableValue const & _disp_x;
    	VariableValue const & _disp_y;
    	VariableValue const & _disp_z;
    	
  	VariableValue const & _J_K;
  	
    	MaterialProperty<Real> const & _J;
    	
	MaterialProperty<RealTensorValue> const & _F;
  	MaterialProperty<RealTensorValue> const & _P;
  	MaterialProperty<RealTensorValue> const & _dPdJK;
  	MaterialProperty<Real> const & _Psinu;
  	
  	MaterialProperty<ElastoPlasticGrowingMaterial *> const & _elastoplasticgrowingMaterial;
  
};
