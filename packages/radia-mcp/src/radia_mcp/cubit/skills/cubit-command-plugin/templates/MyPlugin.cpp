#include "MyPlugin.hpp"
#include "HelloCommand.hpp"
#include "ExportNodeCountCommand.hpp"

// Required: exports cubit_plugin_instance() and the verification marker the
// loader scans for.  Without it Cubit silently skips the file.
CUBIT_PLUGIN(MyPlugin)

MyPlugin::MyPlugin() {}
MyPlugin::~MyPlugin() {}

std::vector<std::string> MyPlugin::get_keys()
{
  // Keys only pair get_keys() with create_command(); users never see them.
  return {"HelloCommand", "ExportNodeCountCommand"};
}

CubitCommand* MyPlugin::create_command(const std::string& key)
{
  if (key == "HelloCommand")
    return new HelloCommand();
  if (key == "ExportNodeCountCommand")
    return new ExportNodeCountCommand();
  return nullptr;
}
