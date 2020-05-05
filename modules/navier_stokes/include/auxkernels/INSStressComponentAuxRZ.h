//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// Forward Declarations
class INSStressComponentAuxRZ;

template <>
InputParameters validParams<INSStressComponentAuxRZ>();

/**
 * Computes h_min / |u|
 */
class INSStressComponentAuxRZ : public AuxKernel
{
public:
  INSStressComponentAuxRZ(const InputParameters & parameters);

  virtual ~INSStressComponentAuxRZ() {}

protected:
  virtual Real computeValue();

  // Velocity gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableValue & _pressure;
  const unsigned _comp;
  const MaterialProperty<Real> & _mu;
};
