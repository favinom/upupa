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
// include della classe madre

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class ChemicalPotential : public Kernel
{
public:
  static InputParameters validParams(); //funzione che legge parametri del Kernel

  ChemicalPotential(const InputParameters & parameters); // costruttore della classe

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
    	
    	int const _dim;
  	Real const _alpha;
  	Real const _kv;
  	
  	int const _id_x;
    	int const _id_y;
    	int const _id_z;  	
  	
	int const _id_JK;
    	int const _id_muk;
  	VariableValue const & _J_K;
  	VariableValue const & _J_K_old;
  	VariableGradient const & _grad_JK;
  	VariableValue const & _mu_K;  
  	
  	MaterialProperty<RealTensorValue> const & _F;
    	MaterialProperty<Real> const & _J;
  
};
