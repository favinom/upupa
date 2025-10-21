//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
//#include "ElastoPlasticGrowingMaterial.h"

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class EvolutionNutrient : public Kernel
{
public:
  static InputParameters validParams();

  EvolutionNutrient(const InputParameters & parameters);

protected:
	int const _dim;
	Real const _k;
	Real const _m;
	Real const _Phisv;
	Real const _d;

    	int const _id_x;
    	int const _id_y;
    	int const _id_z;
    	int const _id_p;
	int const _id_Jk;
    	int const _id_muk;
	int const _id_omega;
    	
    	MaterialProperty<RealTensorValue> const & _F;
    	MaterialProperty<Real> const & _J;
    	MaterialProperty<Real> const & _J_old;
    	MaterialProperty<RealTensorValue> const & _P;
    	
	VariableGradient const & _grad_p;

    	MaterialProperty<Real> const & _J_K;
    	MaterialProperty<Real> const & _J_K_old;
    	
    	VariableValue const & _mu_K;
  	VariableGradient const & _grad_muk;
  	
  	VariableValue const & _omega;
  	VariableGradient const & _grad_omega;
  	VariableValue const & _omega_old;
  
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
};


