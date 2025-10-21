//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBC.h"
#include "PoroElastoGrowth.h"

class CouplingGrowthBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  CouplingGrowthBC(const InputParameters & parameters);

protected:
//  virtual Real computeQpResidual();
//  virtual Real computeQpJacobian();
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  	int const _dim;
    	int const _id_x;
    	int const _id_y;
    	int const _id_z;
	int const _id_JK;

	VariableValue const & _disp_x;
    	VariableValue const & _disp_y;
    	VariableValue const & _disp_z;
    	
  	MaterialProperty<Real> const & _J_K;

    	MaterialProperty<Real> const & _J;
	MaterialProperty<RealTensorValue> const & _F;
  	MaterialProperty<RealTensorValue> const & _P;
  	MaterialProperty<RealTensorValue> const & _dPdJK;
  	MaterialProperty<ElastoPlasticGrowingMaterial *> const & _elastoplasticgrowingMaterial;
};
