//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RankTwoCylindricalComponent.h"
#include "RankTwoScalarTools.h"

registerMooseObject("TensorMechanicsApp", RankTwoCylindricalComponent);

InputParameters
RankTwoCylindricalComponent::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Compute components of a rank-2 tensor in a cylindrical coordinate system");
  params.addRequiredParam<MaterialPropertyName>("rank_two_tensor",
                                                "The rank two material tensor name");
  params.addRequiredParam<std::string>("property_name",
                                       "Name of the material property computed by this model");
  params.addParam<MooseEnum>(
      "cylindrical_component", RankTwoScalarTools::cylindricalOptions(), "Type of scalar output");
  params.addParam<Point>(
      "cylindrical_axis_point1",
      "Start point for determining axis of rotation for cylindrical stress/strain components");
  params.addParam<Point>(
      "cylindrical_axis_point2",
      "End point for determining axis of rotation for cylindrical stress/strain components");
  return params;
}

RankTwoCylindricalComponent::RankTwoCylindricalComponent(const InputParameters & parameters)
  : Material(parameters),
    _tensor(getMaterialProperty<RankTwoTensor>("rank_two_tensor")),
    _property_name(isParamValid("property_name") ? getParam<std::string>("property_name") : ""),
    _property(declareProperty<Real>(_property_name)),
    _cylindrical_component(getParam<MooseEnum>("cylindrical_component")),
    _cylindrical_axis_point1(isParamValid("cylindrical_axis_point1")
                                 ? getParam<Point>("cylindrical_axis_point1")
                                 : Point(0, 0, 0)),
    _cylindrical_axis_point2(isParamValid("cylindrical_axis_point2")
                                 ? getParam<Point>("cylindrical_axis_point2")
                                 : Point(0, 1, 0))
{
}

void
RankTwoCylindricalComponent::initQpStatefulProperties()
{
  _property[_qp] = 0.0;
}

void
RankTwoCylindricalComponent::computeQpProperties()
{
  Point dummy_direction;

  _property[_qp] = RankTwoScalarTools::getCylindricalComponent(_tensor[_qp],
                                                               _cylindrical_component,
                                                               _cylindrical_axis_point1,
                                                               _cylindrical_axis_point2,
                                                               _q_point[_qp],
                                                               dummy_direction);
}
