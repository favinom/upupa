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
#include "Function.h"

class Microstructure : public Material
{
public:
  Microstructure(const InputParameters & parameters);

  static InputParameters validParams();
  
  //virtual RealTensorValue evaluateJac(RealTensorValue const & H, int const & qp) = 0;
  
protected:
 
    	int const _dim;
    	Real const _lambda_p;
    	Real const _sigma_y;
    	
    	RealTensorValue _I;
  	
//  	MaterialProperty<RealTensorValue> const & _F;
  	
  	MaterialProperty<RealTensorValue> const & _F_old;
  	
 // 	MaterialProperty<RealTensorValue> const & _P;
  	
  	MaterialProperty<RealTensorValue> const & _P_old;
  	
  	//MaterialProperty<Microstructure *> & _microstructurematerial;
  	

	

};
