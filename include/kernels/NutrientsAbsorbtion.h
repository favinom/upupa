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
class NutrientsAbsorbtion : public Kernel
{
public:
  static InputParameters validParams(); //funzione che legge parametri del Kernel

  NutrientsAbsorbtion(const InputParameters & parameters); // costruttore della classe

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  	int const _dim;
	Real const _Phisv;
	Real const _gamma_n1;
	Real const _gamma_n2;
	Real const _rhoR;

	int const _id_JK;
	int const _id_cp;
	int const _id_omega;
	
  	VariableValue const & _J_K;
  
  	VariableValue const & _cp;
  	
  	VariableValue const & _omega;
};
