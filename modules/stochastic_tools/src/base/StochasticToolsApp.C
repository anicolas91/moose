//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StochasticToolsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
StochasticToolsApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

registerKnownLabel("StochasticToolsApp");

StochasticToolsApp::StochasticToolsApp(InputParameters parameters) : MooseApp(parameters)
{
  StochasticToolsApp::registerAll(_factory, _action_factory, _syntax);
}

StochasticToolsApp::~StochasticToolsApp() {}

void
StochasticToolsApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  Registry::registerObjectsTo(f, {"StochasticToolsApp"});
  Registry::registerActionsTo(af, {"StochasticToolsApp"});

  // Adds [Trainers] block
  registerSyntaxTask("AddSurrogateAction", "Trainers/*", "add_trainer");
  registerMooseObjectTask("add_trainer", SurrogateTrainer, false);
  addTaskDependency("add_trainer", "add_user_object");

  // Adds [Surrogates] block
  registerSyntaxTask("AddSurrogateAction", "Surrogates/*", "add_surrogate");
  registerMooseObjectTask("add_surrogate", SurrogateModel, false);
  addTaskDependency("add_surrogate", "add_trainer");

  // Adds action for loading Surrogate data
  registerTask("load_surrogate_data", true);
  addTaskDependency("load_surrogate_data", "add_surrogate");
}

void
StochasticToolsApp::registerApps()
{
  registerApp(StochasticToolsApp);
}

void
StochasticToolsApp::registerObjects(Factory & factory)
{
  mooseDeprecated("use registerAll instead of registerObjects");
  Registry::registerObjectsTo(factory, {"StochasticToolsApp"});
}

void
StochasticToolsApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & action_factory)
{
  mooseDeprecated("use registerAll instead of associateSyntax");
  Registry::registerActionsTo(action_factory, {"StochasticToolsApp"});
}

void
StochasticToolsApp::registerExecFlags(Factory & /*factory*/)
{
  mooseDeprecated("use registerAll instead of registerExecFlags");
}

extern "C" void
StochasticToolsApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  StochasticToolsApp::registerAll(f, af, s);
}
extern "C" void
StochasticToolsApp__registerApps()
{
  StochasticToolsApp::registerApps();
}
