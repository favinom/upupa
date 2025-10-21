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

class ElastoPlasticGrowingMaterial : public Material
{
public:
  ElastoPlasticGrowingMaterial(const InputParameters & parameters);

  static InputParameters validParams();
  
  virtual RealTensorValue evaluateJac(RealTensorValue const & H, int const & qp) = 0;
  
  virtual RealTensorValue computeP(RealTensorValue const & F, RealTensorValue const & B_K, Real const & J_K) = 0;
  
  virtual Real computePsi(RealTensorValue const & F, RealTensorValue const & B_K, Real const & J_K) = 0;
  
  virtual RealTensorValue computedJK(RealTensorValue const & F, RealTensorValue const & B_K, Real const & J_K) = 0;

protected:
  	virtual void initQpStatefulProperties() override;
        virtual void computeQpProperties() override;
        
    	int const _dim;
    	
    	VariableGradient const & _grad_disp_x;
    	VariableGradient const & _grad_disp_y;
    	VariableGradient const & _grad_disp_z;
    	
  	RealTensorValue _I;
  	
  	RealTensorValue _B_K;
  	
  	Real _J_K;
  
    	MaterialProperty<RealTensorValue> & _F;
    	
    	MaterialProperty<Real> & _J;
  
  	MaterialProperty<RealTensorValue> & _P;
  	
  	MaterialProperty<RealTensorValue> & _dPdJK; //Proprietà materiale per calcolare il termine dello Jacobiano 
  	
  	MaterialProperty<Real> & _Psi;
  	
  	OptionalMaterialProperty<Real> const & _J_Kmat;
  	
  	OptionalMaterialProperty<RealTensorValue> const & _B_Kmat;
  	
  	MaterialProperty<ElastoPlasticGrowingMaterial *> & _elastoplasticgrowingMaterial;
};
