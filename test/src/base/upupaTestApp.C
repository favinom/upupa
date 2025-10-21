//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "upupaTestApp.h"
#include "upupaApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
upupaTestApp::validParams()
{
  InputParameters params = upupaApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

upupaTestApp::upupaTestApp(InputParameters parameters) : MooseApp(parameters)
{
  upupaTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

upupaTestApp::~upupaTestApp() {}

void
upupaTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  upupaApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"upupaTestApp"});
    Registry::registerActionsTo(af, {"upupaTestApp"});
  }
}

void
upupaTestApp::registerApps()
{
  registerApp(upupaApp);
  registerApp(upupaTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
upupaTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  upupaTestApp::registerAll(f, af, s);
}
extern "C" void
upupaTestApp__registerApps()
{
  upupaTestApp::registerApps();
}
