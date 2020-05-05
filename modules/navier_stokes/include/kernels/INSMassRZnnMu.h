//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSMassRZ.h"

// Forward Declarations
class INSMassRZnnMu;

template <>
InputParameters validParams<INSMassRZnnMu>();

/**
 * This class computes the mass equation residual and Jacobian
 * contributions for the incompressible Navier-Stokes momentum
 * equation in RZ coordinates.  Inherits most of its functionality
 * from INSMass, and calls its computeQpXYZ() functions when
 * necessary.
 */
class INSMassRZnnMu : public INSMassRZ
{
public:
  INSMassRZnnMu(const InputParameters & parameters);
  virtual ~INSMassRZnnMu() {}

protected:
  virtual RealVectorValue dNNmu();
  virtual RealVectorValue dStrongViscDUCompTraction(unsigned comp) override;
  virtual Real dTauDUComp(unsigned comp) override;
};
