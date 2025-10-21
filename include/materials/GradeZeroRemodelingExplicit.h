//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Microstructure.h"
#include "Function.h"

class GradeZeroRemodelingExplicit : public Microstructure
{
public:
  GradeZeroRemodelingExplicit(const InputParameters & parameters);

  static InputParameters validParams();
  

  
 // virtual RealTensorValue evaluateJac(RealTensorValue const & H, int const & qp);
  
protected:
   	virtual void computeQpProperties() override;
  	virtual void initQpStatefulProperties() override;
  	
  	MaterialProperty<Real> & _detBK;
  	MaterialProperty<RealTensorValue> & _B_K;
  	MaterialProperty<RealTensorValue> const & _B_K_old;
  	RealTensorValue _sigma;
  	RealTensorValue _SIGMA;
  /// member variable to hold the computed diffusivity coefficient
  //	MaterialProperty<RealTensorValue> & _B_K;
  	

	

};
