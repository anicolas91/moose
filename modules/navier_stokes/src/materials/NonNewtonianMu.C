//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NonNewtonianMu.h"
#include "libmesh/quadrature.h"

registerMooseObject("NavierStokesApp", NonNewtonianMu);

defineLegacyParams(NonNewtonianMu);

InputParameters
NonNewtonianMu::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  //params.addParam<MaterialPropertyName>("mu_nn_name", "mu_nn", "The name of the dynamic viscosity NonNewtonian");
  params.addParam<Real>("mu",1.,"Viscosity factor");
  params.addParam<Real>("nexp",1.,"Exponential factor in the viscosity power law");
  //params.addCoupledVar("temp", 300.0, "variable for temperature");
  //params.addParam<std::string>("base_name", "Material property base name");
  //params.addParam<Real>("length_scale", 1.0e-9, "Length scale of model");
  //params.addParam<Real>(
  //    "ref_resistivity",
  //    6.5e-6,
  //    "Electrical resistivity of the material at reference temperature in ohm-m.");
  //params.addParam<Real>(
  //    "temp_coeff",
  //    0.0045,
  //    "Temperature coefficient for calculating dependence of resistivity on temp.");
  //params.addParam<Real>("ref_temp", 300.0, "Reference temperature for Electrical resistivity in K");
  return params;
}

NonNewtonianMu::NonNewtonianMu(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    // velocities
    _u_vel(coupledValue("u")),
    _v_vel(coupledValue("v")),
    _w_vel(coupledValue("w")),
    // Gradients
    _grad_u_vel(coupledGradient("u")),
    _grad_v_vel(coupledGradient("v")),
    _grad_w_vel(coupledGradient("w")),
    // constants and mu
    _mu(getParam<Real>("mu")),
    _nexp(getParam<Real>("nexp")),
    _mu_nn(declareProperty<Real>("mu_nn"))
    // _dmu_nn(declareProperty<RealVectorValue>("dmu_nn"))
    //_length_scale(getParam<Real>("length_scale")),
    //_ref_resis(getParam<Real>("ref_resistivity")),
    //_temp_coeff(getParam<Real>("temp_coeff")),
    //_ref_temp(getParam<Real>("ref_temp")),
    //_T(coupledValue("temp")),
    //_base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    //_electric_conductivity(declareProperty<Real>(_base_name + "electrical_conductivity")),
    //_delectric_conductivity_dT(declarePropertyDerivative<Real>(
    //   _base_name + "electrical_conductivity", getVar("temp", 0)->name()))
{
}

void
NonNewtonianMu::computeQpProperties()
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

  auto && shearsq = 2.  * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0)
                  + 2.  * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1)
                  + 2.  * (_u_vel[_qp]/r) * (_u_vel[_qp]/r)
                  + 1.0 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));

  auto && shear = std::sqrt( 2.0/3.0 * shearsq);

  // auto && shear = std::sqrt( 2.0/3.0 *
  //                 ( 2. * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0)
  //                 + 2. * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1)
  //                 + 2. * (_u_vel[_qp]/r) * (_u_vel[_qp]/r)
  //                 + 2.0/2.0 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0))
  //                 ));

// auto && shear = std::abs(_grad_u_vel[_qp](0));

// Moose::out<<"eff shear: "<<shear<<"\n";
// Moose::out<<"only rr shear: "<<shear<<"\n";
// Moose::out<<shear<<"\n";

// calculate mu Derivative

// RealVectorValue mu_der(0, 0, 0);
//
// if (shearsq == 0)
// {
//   mu_der(0) = 1.0;
//   mu_der(1) = 1.0;
// }
// else
// {
//   mu_der(0) = _mu * 0.5 * (_nexp - 1.) * std::pow(2./3., 0.5 * (_nexp - 1.))
//              * std::pow(shearsq, 0.5 * (_nexp - 3.0))
//              * (  4.*_grad_u_vel[_qp](0)*_grad_phi[_j][_qp](0)
//                 + 2.*_grad_u_vel[_qp](1)*_grad_phi[_j][_qp](1)
//                 + 4.*_u_vel[_qp]*_phi[_j][_qp] / (r * r)      );
//   mu_der(1) = _mu * 0.5 * (_nexp - 1.) * std::pow(2./3., 0.5 * (_nexp - 1.))
//              * std::pow(shearsq, 0.5 * (_nexp - 3.0))
//              * (  4.*_grad_v_vel[_qp](1)*_grad_phi[_j][_qp](1)
//                 + 2.*_grad_v_vel[_qp](0)*_grad_phi[_j][_qp](0) );
// }
//
// _dmu_nn = mu_der;

// if (shear <= 0.5)
  if (shear <= 1e-100)
  {
  	_mu_nn[_qp] = _mu;
  }
  else
  {
    _mu_nn[_qp] = _mu *  std::pow( shear , _nexp-1.0);
  }

// Moose::out<<"viscosity mu final: "<<_mu_nn[_qp]<<"\n";

  //const Real ref_resis = _ref_resis / _length_scale;
  //const Real resistivity = ref_resis * (1.0 + _temp_coeff * (_T[_qp] - _ref_temp));
  //const Real dresistivity_dT = ref_resis * _temp_coeff;
  //_electric_conductivity[_qp] = 1.0 / resistivity;
  //  _delectric_conductivity_dT[_qp] = -1.0 / (resistivity * resistivity) * dresistivity_dT;
}
