//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

class ElastoPlasticMaterial : public Material
{
public:
  ElastoPlasticMaterial(const InputParameters & parameters);

  static InputParameters validParams();
  
  virtual RealTensorValue evaluateJac(RealTensorValue const & H, int const & qp) = 0;
  
  virtual RealTensorValue computeP(RealTensorValue const & F, RealTensorValue const & B_K) = 0;
  
  virtual Real computePsi(RealTensorValue const & F, RealTensorValue const & B_K) = 0;

protected:
  	virtual void initQpStatefulProperties() override;
        virtual void computeQpProperties() override;
        
    	int const _dim;
    	
    	VariableGradient const & _grad_disp_x;
    	VariableGradient const & _grad_disp_y;
    	VariableGradient const & _grad_disp_z;
    	
  	RealTensorValue _I;
  	
  	RealTensorValue _B_K;
  
    	MaterialProperty<RealTensorValue> & _F;
    	MaterialProperty<Real> & _J;
    	
//    	MaterialProperty<RealTensorValue> const & _F_old;
  
  	MaterialProperty<RealTensorValue> & _P;
  	
  	MaterialProperty<Real> & _Psi;
  	
//  	MaterialProperty<RealTensorValue> const & _P_old;
  	
//  	MaterialProperty<RealTensorValue> & _sigma;
  	
//  	MaterialProperty<RealTensorValue> & _SIGMA;
  	
  	OptionalMaterialProperty<RealTensorValue> const & _B_Kmat;
  	
  	MaterialProperty<ElastoPlasticMaterial *> & _elastoplasticMaterial;

	

};
