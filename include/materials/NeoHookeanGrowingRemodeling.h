//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElastoPlasticGrowingMaterial.h"

class NeoHookeanGrowingRemodeling : public ElastoPlasticGrowingMaterial
{
public:
  NeoHookeanGrowingRemodeling(const InputParameters & parameters);

  static InputParameters validParams();
  
  virtual RealTensorValue evaluateJac(RealTensorValue const & H, int const & qp);
  
  virtual RealTensorValue computeP(RealTensorValue const & F, RealTensorValue const & B_K, Real const & J_K);
  
  virtual Real computePsi(RealTensorValue const & F, RealTensorValue const & B_K, Real const & J_K);
  
  virtual RealTensorValue computedJK(RealTensorValue const & F, RealTensorValue const & B_K, Real const & J_K);

protected:
 // virtual void computeQpProperties() override;
	Real const _lambda;
	Real const _mu;
	Real const _Phisnu;
};
