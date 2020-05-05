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
#include "DerivativeMaterialInterface.h"

// Forward Declarations
class NonNewtonianMu;

template <>
InputParameters validParams<NonNewtonianMu>();

/**
 * Calculates resistivity and electrical conductivity as a function of temperature.
 * It is assumed that resistivity varies linearly with temperature.
 */
class NonNewtonianMu : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  NonNewtonianMu(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

private:
  //const Real _length_scale;
  //const Real _ref_resis;
  //const Real _temp_coeff;
  //const Real _ref_temp;
  //const VariableValue & _T;
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;
  const Real _mu;
  const Real _nexp;
  MaterialProperty<Real> & _mu_nn;
  // MaterialProperty<RealVectorValue> & _dmu_nn;
  // MaterialProperty<Real> & _dmu_nn;
  //const std::string _base_name;
  //MaterialProperty<Real> & _electric_conductivity;
//  MaterialProperty<Real> & _delectric_conductivity_dT;
};
