//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSStressComponentAuxRZ.h"
#include "MooseMesh.h"

registerMooseObject("NavierStokesApp", INSStressComponentAuxRZ);

template <>
InputParameters
validParams<INSStressComponentAuxRZ>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addClassDescription("This class computes the stress component based on "
                             "pressure and velocity for incompressible Navier-Stokes");
  params.addCoupledVar("u", 0, "x-velocity (radial)");
  params.addCoupledVar("v", 0, "y-velocity (axial)"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity (hoop)"); // only required in 2D and 3D
  // params.addCoupledVar("velocity", "The velocity component");
  params.addCoupledVar("pressure", 0, "The pressure");
  params.addRangeCheckedParam<unsigned int>("comp", 0, "0<=comp<=2", "The component");
  params.addParam<MaterialPropertyName>("mu_name", "mu", "The viscosity");

  return params;
}

INSStressComponentAuxRZ::INSStressComponentAuxRZ(const InputParameters & parameters)
  : AuxKernel(parameters),
    _grad_u_vel(isCoupled("u") ? coupledGradient("u") : _grad_zero),
    _grad_v_vel(isCoupled("v") ? coupledGradient("v") : _grad_zero),
    _grad_w_vel(isCoupled("w") ? coupledGradient("w") : _grad_zero),
    _u_vel(coupledValue("u")),
    _v_vel(coupledValue("v")),
    _w_vel(coupledValue("w")),
    _pressure(coupledValue("pressure")),
    _comp(getParam<unsigned>("comp")),
    _mu(getMaterialProperty<Real>("mu_name"))
{
}

Real
INSStressComponentAuxRZ::computeValue()
{
  auto && s = -_pressure[_qp] + 2*_mu[_qp] * _grad_u_vel[_qp](_comp); //
  // Moose::out<<"shear strain dv/dr: "<<_grad_u_vel[_qp](_comp)<<"\n";
  // Moose::out<<_grad_u_vel[_qp](_comp)<<"\n";


  if (_comp == 1)
  {
    s = -_pressure[_qp] +2* _mu[_qp] * _grad_v_vel[_qp](_comp);
  }

  // If this is the hoop component there is an extra term for RZ.
  else if (_comp == 2)
  {
    const Real r = _q_point[_qp](0);
    s = -_pressure[_qp] + 2*_mu[_qp] * _u_vel[_qp] / r;
  }

return s;
}
