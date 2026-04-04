#include "ExportNetgenCommand.hpp"
#include "MeshData.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#include <process.h>   // _popen/_pclose
#include <direct.h>
#include <io.h>
#define popen  _popen
#define pclose _pclose
#else
#include <unistd.h>
#endif

ExportNetgenCommand::ExportNetgenCommand() {}
ExportNetgenCommand::~ExportNetgenCommand() {}

std::vector<std::string> ExportNetgenCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "export netgen <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1-5>'>] "
    "[overwrite]"
  );
  return syntax_list;
}

std::vector<std::string> ExportNetgenCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "export netgen \"filename.vol\" [order {1|2|3|4|5}] [overwrite]"
  );
  return help;
}

std::vector<std::string> ExportNetgenCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Export mesh as Netgen .vol with high-order curving and labels.\n"
    "Runs external Python (calc_export_netgen.py) for curving.\n"
    "Produces .vol file + companion .json with CAD reference values.\n"
    "Order 1-5 supported (default: 2)."
  );
  return help;
}

// Find external Python 3.12 (not Cubit's 3.10)
static std::string find_python()
{
  // RADIA_PYTHON env var
  const char *env = std::getenv("RADIA_PYTHON");
  if (env && env[0]) return env;

#ifdef _WIN32
  // Try py -3
  FILE *p = popen("py -3 -c \"print('ok')\" 2>nul", "r");
  if (p) {
    char buf[16] = {};
    if (fgets(buf, sizeof(buf), p) && std::string(buf).find("ok") != std::string::npos) {
      pclose(p);
      return "py -3";
    }
    pclose(p);
  }
#endif
  return "python";
}

// Find calc_export_netgen.py via python -c "import radia.panels.calc_export_netgen..."
static std::string find_script(const std::string &python)
{
  std::string cmd = python + " -c \"import radia.panels.calc_export_netgen as m; "
    "import os; print(os.path.abspath(m.__file__))\" 2>nul";
  FILE *p = popen(cmd.c_str(), "r");
  if (!p) return "";
  char buf[1024] = {};
  std::string result;
  while (fgets(buf, sizeof(buf), p))
    result += buf;
  pclose(p);
  // Trim
  while (!result.empty() && (result.back() == '\n' || result.back() == '\r'))
    result.pop_back();
  return result;
}

bool ExportNetgenCommand::execute(CubitCommandData &data)
{
  std::string filename;
  data.get_string("filename", filename);
  if (filename.empty()) {
    PRINT_ERROR("Filename required.\n");
    return false;
  }

  int order = 2;
  data.get_value("order", order);
  if (order < 1 || order > 5) {
    PRINT_ERROR("Order must be 1-5.\n");
    return false;
  }

  // Save model to temp cub5 (silent)
  std::string cub5;
#ifdef _WIN32
  char tmpdir[MAX_PATH];
  GetTempPathA(MAX_PATH, tmpdir);
  cub5 = std::string(tmpdir) + "radia_netgen_export.cub5";
#else
  cub5 = "/tmp/radia_netgen_export_" + std::to_string(getpid()) + ".cub5";
#endif
  // Replace backslashes
  for (char &c : cub5) if (c == '\\') c = '/';

  CubitInterface::silent_cmd(
    ("save cub5 \"" + cub5 + "\" overwrite").c_str());

  // Find external Python + script
  std::string python = find_python();
  std::string script = find_script(python);
  if (script.empty()) {
    PRINT_ERROR("Cannot find calc_export_netgen.py. Is radia installed?\n");
    std::remove(cub5.c_str());
    return false;
  }

  // Build command
  std::string cmd = python + " \"" + script + "\" --cub5 \"" + cub5
    + "\" --order " + std::to_string(order) + " --vol \"" + filename + "\"";

  PRINT_INFO("Exporting Netgen .vol (order %d)...\n", order);

  // Run subprocess
  FILE *proc = popen(cmd.c_str(), "r");
  if (!proc) {
    PRINT_ERROR("Failed to start Python subprocess.\n");
    std::remove(cub5.c_str());
    return false;
  }

  // Capture output (last line is JSON)
  std::string output;
  char buf[4096];
  while (fgets(buf, sizeof(buf), proc))
    output += buf;
  int rc = pclose(proc);

  // Cleanup temp cub5
  std::remove(cub5.c_str());

  if (rc != 0) {
    PRINT_ERROR("Netgen export failed (exit %d).\n", rc);
    if (!output.empty())
      PRINT_ERROR("%s\n", output.substr(0, 2000).c_str());
    return false;
  }

  // Parse last line for summary
  std::string last_line;
  std::istringstream iss(output);
  std::string line;
  while (std::getline(iss, line))
    if (!line.empty()) last_line = line;

  PRINT_INFO("Exported Netgen Vol (order %d): %s\n", order, filename.c_str());

  return true;
}
