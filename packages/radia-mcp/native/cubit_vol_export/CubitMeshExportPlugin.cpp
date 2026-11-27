#include "CubitMeshExportPlugin.hpp"
#include "ExportNetgenCommand.hpp"

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
  // Optional plugin-load trace: written only when RADIA_CUBIT_EXPORT_DEBUG_LOG
  // names a file; no fixed machine path is used.
  const char* path = std::getenv("RADIA_CUBIT_EXPORT_DEBUG_LOG");
  if (!path || !*path) return;
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
  keys.push_back("ExportNetgenCommand");
  return keys;
}

CubitCommand* CubitMeshExportPlugin::create_command(const std::string &key)
{
  if (key == "ExportNetgenCommand")
    return new ExportNetgenCommand();
  return nullptr;
}

