//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSMomentumTractionFormRZ.h"

// Forward Declarations

/**
 * This class computes additional momentum equation residual and
 * Jacobian contributions for the incompressible Navier-Stokes
 * momentum equation in RZ (axisymmetric cylindrical) coordinates.
 */
class INSMomentumTractionFormRZnnMu : public INSMomentumTractionFormRZ
{
public:
  static InputParameters validParams();

  INSMomentumTractionFormRZnnMu(const InputParameters & parameters);

  virtual ~INSMomentumTractionFormRZnnMu() {}

protected:
  virtual RealVectorValue dNNmu();
  virtual Real computeQpJacobian() override;
};
