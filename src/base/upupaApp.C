#include "upupaApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
upupaApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

upupaApp::upupaApp(InputParameters parameters) : MooseApp(parameters)
{
  upupaApp::registerAll(_factory, _action_factory, _syntax);
}

upupaApp::~upupaApp() {}

void 
upupaApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<upupaApp>(f, af, s);
  Registry::registerObjectsTo(f, {"upupaApp"});
  Registry::registerActionsTo(af, {"upupaApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
upupaApp::registerApps()
{
  registerApp(upupaApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
upupaApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  upupaApp::registerAll(f, af, s);
}
extern "C" void
upupaApp__registerApps()
{
  upupaApp::registerApps();
}
