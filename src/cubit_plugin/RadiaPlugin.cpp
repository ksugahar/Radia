#include "RadiaPlugin.hpp"
#include "ExportNetgenCommand.hpp"
#include "ExportNastranCommand.hpp"
#include "ExportGmshCommand.hpp"
#include "ExportMegCommand.hpp"
#include "ExportVtkCommand.hpp"

// ============================================================
// Python API plugin (CUBIT_PLUGIN_DIR / cubit.init)
// ============================================================
CUBIT_PLUGIN(RadiaPlugin)

RadiaPlugin::RadiaPlugin() {}
RadiaPlugin::~RadiaPlugin() {}

std::vector<std::string> RadiaPlugin::get_keys()
{
  std::vector<std::string> keys;
  keys.push_back("ExportNetgenCommand");
  keys.push_back("ExportNastranCommand");
  keys.push_back("ExportGmshCommand");
  keys.push_back("ExportMegCommand");
  keys.push_back("ExportVtkCommand");
  return keys;
}

CubitCommand* RadiaPlugin::create_command(const std::string &key)
{
  if (key == "ExportNetgenCommand")
    return new ExportNetgenCommand();
  else if (key == "ExportNastranCommand")
    return new ExportNastranCommand();
  else if (key == "ExportGmshCommand")
    return new ExportGmshCommand();
  else if (key == "ExportMegCommand")
    return new ExportMegCommand();
  else if (key == "ExportVtkCommand")
    return new ExportVtkCommand();
  return nullptr;
}

