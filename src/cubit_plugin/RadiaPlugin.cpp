#include "RadiaPlugin.hpp"
#ifdef HAVE_NETGEN
#include "ExportNetgenCommand.hpp"
#endif
#include "ExportNastranCommand.hpp"
#include "ExportGmshCommand.hpp"
#include "ExportVtkCommand.hpp"
#include "CoilCommand.hpp"

// ============================================================
// Python API plugin (CUBIT_PLUGIN_DIR / cubit.init)
// ============================================================
CUBIT_PLUGIN(RadiaPlugin)

static void dbglog(const char* msg) {
  FILE* f = fopen("C:\\compact_netgen_debug.log", "a");
  if (f) { fprintf(f, "%s\n", msg); fclose(f); }
}

RadiaPlugin::RadiaPlugin() {
  dbglog("RadiaPlugin constructor called");
}
RadiaPlugin::~RadiaPlugin() {}

std::vector<std::string> RadiaPlugin::get_keys()
{
  dbglog("get_keys() called");
  std::vector<std::string> keys;
#ifdef HAVE_NETGEN
  keys.push_back("ExportNetgenCommand");
#endif
  keys.push_back("ExportNastranCommand");
  keys.push_back("ExportGmshCommand");
  keys.push_back("ExportVtkCommand");
  keys.push_back("CoilCommand");
  return keys;
}

CubitCommand* RadiaPlugin::create_command(const std::string &key)
{
#ifdef HAVE_NETGEN
  if (key == "ExportNetgenCommand")
    return new ExportNetgenCommand();
  else
#endif
  if (key == "ExportNastranCommand")
    return new ExportNastranCommand();
  else if (key == "ExportGmshCommand")
    return new ExportGmshCommand();
  else if (key == "ExportVtkCommand")
    return new ExportVtkCommand();
  else if (key == "CoilCommand")
    return new CoilCommand();
  return nullptr;
}

