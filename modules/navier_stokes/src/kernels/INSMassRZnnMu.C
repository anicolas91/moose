//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSMassRZnnMu.h"

registerMooseObject("NavierStokesApp", INSMassRZnnMu);

template <>
InputParameters
validParams<INSMassRZnnMu>()
{
  InputParameters params = validParams<INSMassRZ>();
  params.addClassDescription("This class computes the mass equation residual and Jacobian "
                             "contributions for the incompressible Navier-Stokes momentum equation "
                             "in RZ coordinates and a variable viscosity");
  return params;
}

INSMassRZnnMu::INSMassRZnnMu(const InputParameters & parameters) : INSMassRZ(parameters) {}


RealVectorValue
INSMassRZnnMu::dNNmu()
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


RealVectorValue
INSMassRZnnMu::dStrongViscDUCompTraction(unsigned comp)
{
  RealVectorValue add_jac(0, 0, 0);

  if (comp == 0)
  {
    add_jac(0) = (INSMassRZ::strongViscousTermTraction()(0)/_mu[_qp] )
                 * dNNmu()(0);
    add_jac(1) = (INSMassRZ::strongViscousTermTraction()(1)/_mu[_qp] )
                 * dNNmu()(0);
  }
  else if (comp == 1)
    add_jac(0) = (INSMassRZ::strongViscousTermTraction()(0)/_mu[_qp] )
                 * dNNmu()(1);
    add_jac(1) = (INSMassRZ::strongViscousTermTraction()(1)/_mu[_qp] )
                 * dNNmu()(1);

  return INSMassRZ::dStrongViscDUCompTraction(comp) + add_jac;
}


Real
INSMassRZnnMu::dTauDUComp(unsigned comp)
{
  Real nu = _mu[_qp] / _rho[_qp];
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real h = _current_elem->hmax();
  Real transient_part = _transient_term ? 4. / (_dt * _dt) : 0.;
  return -_alpha / 2. *
         std::pow(transient_part + (2. * U.norm() / h) * (2. * U.norm() / h) +
                      9. * (4. * nu / (h * h)) * (4. * nu / (h * h)),
                  -1.5) *
         ( 2. * (2. * U.norm() / h) * 2. / h * U(comp) * _phi[_j][_qp]/
           (U.norm() + std::numeric_limits<double>::epsilon())
         + 2. * 9. * (4. * nu / (h * h)) * 4. *dNNmu()(comp) / (_rho[_qp] * h * h)
         );
}
