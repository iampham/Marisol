#include "crystal_plasticityApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
crystal_plasticityApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

crystal_plasticityApp::crystal_plasticityApp(InputParameters parameters) : MooseApp(parameters)
{
  crystal_plasticityApp::registerAll(_factory, _action_factory, _syntax);
}

crystal_plasticityApp::~crystal_plasticityApp() {}

void 
crystal_plasticityApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<crystal_plasticityApp>(f, af, s);
  Registry::registerObjectsTo(f, {"crystal_plasticityApp"});
  Registry::registerActionsTo(af, {"crystal_plasticityApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
crystal_plasticityApp::registerApps()
{
  registerApp(crystal_plasticityApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
crystal_plasticityApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  crystal_plasticityApp::registerAll(f, af, s);
}
extern "C" void
crystal_plasticityApp__registerApps()
{
  crystal_plasticityApp::registerApps();
}
