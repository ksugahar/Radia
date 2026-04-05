#include "CoilCommand.hpp"
#include "CubitMessage.hpp"
#include "CubitInterface.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#include <direct.h>  // _getcwd
#include <io.h>
#define pclose _pclose
#define popen  _popen
#endif

CoilCommand::CoilCommand() {}
CoilCommand::~CoilCommand() {}

std::vector<std::string> CoilCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "coil <string:label='script',help='<script.py>'> "
    "[output <string:label='output',help='<path.step>'>] "
    "[noimport]"
  );
  return syntax_list;
}

std::vector<std::string> CoilCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "coil \"script.py\" [output \"path.step\"] [noimport]"
  );
  return help;
}

std::vector<std::string> CoilCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Generate coil geometry from a Python script using CoilBuilder.\n\n"
    "The script must define a build_coil() function that returns a\n"
    "CoilBuilder instance. The coil is exported as STEP and imported\n"
    "into the current Cubit session.\n\n"
    "Requires external Python 3.12 with NGSolve/OCC.\n"
    "Set RADIA_PYTHON environment variable to override Python path.\n\n"
    "Options:\n"
    "  output \"path.step\"  Output STEP file path (default: coil.step in CWD)\n"
    "  noimport             Do not import STEP into Cubit\n\n"
    "Example:\n"
    "  coil \"my_coil.py\"\n"
    "  coil \"my_coil.py\" output \"C:/models/coil.step\"\n"
    "  coil \"my_coil.py\" noimport\n"
  );
  return help;
}

// Find external Python 3.12 (not Cubit's embedded 3.10).
// Priority: RADIA_PYTHON env > py -3 (Windows) > python
static std::string find_python()
{
  const char *env = std::getenv("RADIA_PYTHON");
  if (env && env[0])
    return std::string(env);

#ifdef _WIN32
  // Test py -3 launcher
  FILE *fp = popen("py -3 -c \"print('ok')\" 2>nul", "r");
  if (fp) {
    char buf[64] = {};
    if (fgets(buf, sizeof(buf), fp) && strstr(buf, "ok")) {
      pclose(fp);
      return "py -3";
    }
    pclose(fp);
  }
#endif

  return "python";
}

// Find generate_coil.py via radia package
static std::string find_generate_coil(const std::string &python)
{
  std::string cmd = python
    + " -c \"import radia.panels.generate_coil as m; import os; "
      "print(os.path.abspath(m.__file__))\" 2>nul";
#ifndef _WIN32
  // Replace 2>nul with 2>/dev/null on non-Windows
  cmd = python
    + " -c \"import radia.panels.generate_coil as m; import os; "
      "print(os.path.abspath(m.__file__))\" 2>/dev/null";
#endif

  FILE *fp = popen(cmd.c_str(), "r");
  if (!fp) return "";

  char buf[1024] = {};
  std::string result;
  while (fgets(buf, sizeof(buf), fp))
    result += buf;
  pclose(fp);

  // Trim whitespace
  while (!result.empty() && (result.back() == '\n' || result.back() == '\r'
                              || result.back() == ' '))
    result.pop_back();

  return result;
}

bool CoilCommand::execute(CubitCommandData &data)
{
  // Parse arguments
  std::string script_path;
  if (!data.get_string("script", script_path) || script_path.empty()) {
    PRINT_ERROR("coil: script path is required.\n");
    return false;
  }

  std::string output_path;
  data.get_string("output", output_path);

  bool no_import = data.find_keyword("noimport");

  // Default output: coil.step in CWD
  if (output_path.empty()) {
    char cwd[1024] = {};
    _getcwd(cwd, sizeof(cwd));
    output_path = std::string(cwd) + "/coil.step";
  }

  // Replace backslashes for consistency
  for (auto &c : script_path) if (c == '\\') c = '/';
  for (auto &c : output_path) if (c == '\\') c = '/';

  // Find Python and generate_coil.py
  std::string python = find_python();
  std::string gen_script = find_generate_coil(python);
  if (gen_script.empty()) {
    PRINT_ERROR("coil: Cannot find generate_coil.py. Is radia installed?\n");
    return false;
  }

  // Build subprocess command
  std::string cmd = python
    + " \"" + gen_script + "\""
    + " --script \"" + script_path + "\""
    + " --output \"" + output_path + "\"";

  PRINT_INFO("Generating coil: %s\n", cmd.c_str());

  // Run subprocess
  FILE *fp = popen(cmd.c_str(), "r");
  if (!fp) {
    PRINT_ERROR("coil: Failed to run Python subprocess.\n");
    return false;
  }

  // Read output
  std::string output;
  char buf[1024];
  while (fgets(buf, sizeof(buf), fp))
    output += buf;

  int ret = pclose(fp);
  if (ret != 0) {
    PRINT_ERROR("coil: Python subprocess failed (exit code %d).\n"
                "Output: %s\n", ret, output.c_str());
    return false;
  }

  // Find JSON in last non-empty line
  std::istringstream iss(output);
  std::string line, last_line;
  while (std::getline(iss, line)) {
    if (!line.empty())
      last_line = line;
  }

  // Report result
  if (last_line.find('{') != std::string::npos) {
    PRINT_INFO("Coil generated: %s\n", last_line.c_str());
  }

  // Check STEP file exists
  std::ifstream check(output_path);
  if (!check.good()) {
    PRINT_ERROR("coil: STEP file not created: %s\n", output_path.c_str());
    return false;
  }
  check.close();

  // Import STEP into Cubit
  if (!no_import) {
    std::string import_cmd = "import step \"" + output_path + "\" heal";
    PRINT_INFO("Importing: %s\n", import_cmd.c_str());
    CubitInterface::cmd(import_cmd.c_str());
  }

  PRINT_INFO("coil: Done.\n");
  return true;
}
