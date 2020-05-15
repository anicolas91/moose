//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSADNonNewtonianMu.h"
#include "Function.h"
#include "Assembly.h"
#include "NonlinearSystemBase.h"
#include "INSADObjectTracker.h"

registerMooseObject("NavierStokesApp", INSADNonNewtonianMu);

InputParameters
INSADNonNewtonianMu::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addRequiredCoupledVar("velocity", "The velocity");
  params.addParam<Real>("mu_in",1.,"Viscosity factor");
  params.addParam<Real>("nexp",1.,"Exponential factor in the viscosity power law");
  return params;
}

INSADNonNewtonianMu::INSADNonNewtonianMu(const InputParameters & parameters)
  : ADMaterial(parameters),
    // velocities and gradients
    _velocity(adCoupledVectorValue("velocity")),
    _grad_velocity(adCoupledVectorGradient("velocity")),
    // constants and mu
    _mu_in(getParam<Real>("mu_in")),
    _nexp(getParam<Real>("nexp")),
    _mu(declareADProperty<Real>("mu")) //converted from AD
{
}

void
INSADNonNewtonianMu::computeQpProperties()
{
  // In here the strain rate tensor, when used to calculate the shear rate, is
  // multiplied by 2/3 instead of only 2 to accomodate for Von Mises
	// auto && shear = std::sqrt( 2./3. *
  //                 ( 1. * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0)
  //                 + 1. * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1)
  //                 + 1. * _grad_w_vel[_qp](2) * _grad_w_vel[_qp](2)
  //                 + 1./2.0 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0))
  //                 + 1./2.0 * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0)) * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0))
  //                 + 1./2.0 * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1)) * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1))
  //                 ));

  // auto && shear = std::sqrt(4. * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0)
  //                 + 4. * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1)
  //                 + 4. * _grad_w_vel[_qp](2) * _grad_w_vel[_qp](2)
  //                 + 4. * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0))
  //                 + 4. * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0)) * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0))
  //                 + 4. * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1)) * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1))
  //                   );

  const Real r = _q_point[_qp](0);

  auto && shearsq = 2.  * _grad_velocity[_qp](0,0) * _grad_velocity[_qp](0,0)
                  + 2.  * _grad_velocity[_qp](1,1) * _grad_velocity[_qp](1,1)
                  + 2.  * (_velocity[_qp](0)/r) * (_velocity[_qp](0)/r)
                  + 1.0 * (_grad_velocity[_qp](0,1) + _grad_velocity[_qp](1,0)) * (_grad_velocity[_qp](0,1) + _grad_velocity[_qp](1,0));

  auto && shear = std::sqrt( 2.0/3.0 * shearsq);

// Moose::out<<"eff shear: "<<shear<<"\n";
// Moose::out<<"only rr shear: "<<shear<<"\n";
// Moose::out<<shear<<"\n";

  if (shear <= 1e-100)
  {
  	_mu[_qp] = _mu_in;
  }
  else
  {
    _mu[_qp] = _mu_in *  std::pow( shear , _nexp-1.0);
  }

}
