#include "CubitToVolExportPlugin.hpp"
#include "ExportVolCommand.hpp"

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
CUBIT_PLUGIN(CubitToVolExportPlugin)

static void dbglog(const char* msg) {
  // Optional plugin-load trace: written only when CUBIT_TO_VOL_EXPORT_DEBUG_LOG
  // names a file; no fixed machine path is used.
  const char* path = std::getenv("CUBIT_TO_VOL_EXPORT_DEBUG_LOG");
  if (!path || !*path) return;
  FILE* f = fopen(path, "a");
  if (f) { fprintf(f, "%s\n", msg); fclose(f); }
}

CubitToVolExportPlugin::CubitToVolExportPlugin() {
  dbglog("CubitToVolExportPlugin constructor called");
}
CubitToVolExportPlugin::~CubitToVolExportPlugin() {}

std::vector<std::string> CubitToVolExportPlugin::get_keys()
{
  dbglog("get_keys() called");
  std::vector<std::string> keys;
  keys.push_back("ExportVolCommand");
  return keys;
}

CubitCommand* CubitToVolExportPlugin::create_command(const std::string &key)
{
  if (key == "ExportVolCommand")
    return new ExportVolCommand();
  return nullptr;
}

