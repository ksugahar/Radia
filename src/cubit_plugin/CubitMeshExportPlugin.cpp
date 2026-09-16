#include "CubitMeshExportPlugin.hpp"
#ifdef HAVE_NETGEN
#include "ExportNetgenCommand.hpp"
#endif
#include "ExportNastranCommand.hpp"
#include "ExportGmshCommand.hpp"
#include "ExportVtkCommand.hpp"
#include "ExportFemeemCommand.hpp"
#include "ExportMegCommand.hpp"

#include <cstdio>
#include <cstdlib>
#ifdef _WIN32
#include <windows.h>  // MAX_PATH
#else
#include <climits>
#ifndef MAX_PATH
#define MAX_PATH 260
#endif
#endif

// ============================================================
// Python API plugin (CUBIT_PLUGIN_DIR / cubit.init)
// ============================================================
CUBIT_PLUGIN(CubitMeshExportPlugin)

static void dbglog(const char* msg) {
  // Per-user filename under C:\temp: shared machine-level log at
  // C:\ root was not writable by non-admin users (2026-04-24).
  const char* user = std::getenv("USERNAME");
  if (!user || !*user) user = "unknown";
  char path[MAX_PATH];
  snprintf(path, sizeof(path), "C:\\temp\\compact_netgen_debug_%s.log", user);
  FILE* f = fopen(path, "a");
  if (f) { fprintf(f, "%s\n", msg); fclose(f); }
}

CubitMeshExportPlugin::CubitMeshExportPlugin() {
  dbglog("CubitMeshExportPlugin constructor called");
}
CubitMeshExportPlugin::~CubitMeshExportPlugin() {}

std::vector<std::string> CubitMeshExportPlugin::get_keys()
{
  dbglog("get_keys() called");
  std::vector<std::string> keys;
#ifdef HAVE_NETGEN
  keys.push_back("ExportNetgenCommand");
#endif
  keys.push_back("ExportNastranCommand");
  keys.push_back("ExportGmshCommand");
  keys.push_back("ExportVtkCommand");
  keys.push_back("ExportFemeemCommand");
  keys.push_back("ExportMegCommand");
  return keys;
}

CubitCommand* CubitMeshExportPlugin::create_command(const std::string &key)
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
  else if (key == "ExportFemeemCommand")
    return new ExportFemeemCommand();
  else if (key == "ExportMegCommand")
    return new ExportMegCommand();
  return nullptr;
}

