//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSMomentumTractionFormRZnnMu.h"

registerMooseObject("NavierStokesApp", INSMomentumTractionFormRZnnMu);

InputParameters
INSMomentumTractionFormRZnnMu::validParams()
{
  InputParameters params = INSMomentumTractionFormRZ::validParams();
  params.addClassDescription("This class computes additional momentum equation residual and "
                             "Jacobian contributions for the incompressible Navier-Stokes momentum "
                             "equation in RZ (axisymmetric cylindrical) coordinates.");
  return params;
}

INSMomentumTractionFormRZnnMu::INSMomentumTractionFormRZnnMu(const InputParameters & parameters)
  : INSMomentumTractionFormRZ(parameters)
{
}

RealVectorValue
INSMomentumTractionFormRZnnMu::dNNmu()
{
  // This entire calculation is temporary, this should be calculated on the material,
 // but to ensure this part is not a problem, it goes here. Weee.
  const Real & r = _q_point[_qp](0);
  RealVectorValue mu_der(0, 0, 0);
  const Real & m = 781864.16;
  const Real & n = 0.4;

  auto && shearsq = 2. * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0)
                  + 2. * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1)
                  + 2. * (_u_vel[_qp]/r) * (_u_vel[_qp]/r)
                  + (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));

  if (shearsq == 0)
  {
    mu_der(0) = 0.0;
    mu_der(1) = 0.0;
  }
  else
  {
    mu_der(0) = m * 0.5 * (n - 1.) * std::pow(2./3., 0.5 * (n - 1.))
               * std::pow(shearsq, 0.5 * (n - 3.0))
               * (  4.*_grad_u_vel[_qp](0)*_grad_phi[_j][_qp](0)
                  + 4.*_u_vel[_qp]*_phi[_j][_qp] / (r * r)
                  + 2.*_grad_u_vel[_qp](1)*_grad_phi[_j][_qp](1)
                        );
    mu_der(1) = m * 0.5 * (n - 1.) * std::pow(2./3., 0.5 * (n - 1.))
               * std::pow(shearsq, 0.5 * (n - 3.0))
               * (  4.*_grad_v_vel[_qp](1)*_grad_phi[_j][_qp](1)
                  + 2.*_grad_v_vel[_qp](0)*_grad_phi[_j][_qp](0) );
  }
  return mu_der;
}

Real
INSMomentumTractionFormRZnnMu::computeQpJacobian()
{
  // Initial jacobian (mu * A')
  Real jac_base = INSMomentumTractionFormRZ::computeQpJacobian();

  // Extra bit corresponding to mu'
  if (_component == 0)
  {
    jac_base += (INSMomentumTractionFormRZ::computeQpResidual() - _test[_i][_qp] * strongPressureTerm()(0)) / _mu[_qp]
                 * dNNmu()(0);
  }
  else if (_component == 1)
    jac_base += (INSMomentumTractionFormRZ::computeQpResidual() - _test[_i][_qp] * strongPressureTerm()(1)) / _mu[_qp]
                 * dNNmu()(1);
  // Moose::out<<"jacobian base: "<<jac_base<<"\n";
  return jac_base;
}
