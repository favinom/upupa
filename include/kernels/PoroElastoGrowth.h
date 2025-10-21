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
#include "ElastoPlasticGrowingMaterial.h"
#include "GradeOneGrowth.h"
/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class PoroElastoGrowth : public Kernel
{
public:
  static InputParameters validParams();

  PoroElastoGrowth(const InputParameters & parameters);

protected:
	int const _dim;
	int const _component;
	VariableValue const & _p;
    	int const _id_x;
    	int const _id_y;
    	int const _id_z;
	int const _id_p;
	int const _id_J_K;
    
    
        MaterialProperty<RealTensorValue> const & _F;
    	MaterialProperty<Real> const & _J;
    	MaterialProperty<Real> const & _J_K;
    	MaterialProperty<RealTensorValue> const & _P;
    	MaterialProperty<RealTensorValue> const & _dPdJK;
     
  MaterialProperty<ElastoPlasticGrowingMaterial *> const & _elastoplasticgrowingMaterial;
  
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
};


