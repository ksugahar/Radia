#include "VerifyLauncherCommand.hpp"
#include "LauncherLogic.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"

#include <iostream>
#include <sstream>

VerifyLauncherCommand::VerifyLauncherCommand() {}
VerifyLauncherCommand::~VerifyLauncherCommand() {}

std::vector<std::string> VerifyLauncherCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  // radia_export verify_launcher [require <names>...]
  syntax_list.push_back(
    "radia_export verify_launcher "
    "[require <string:label='require',help='comma-separated required labels'>]"
  );
  return syntax_list;
}

std::vector<std::string> VerifyLauncherCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "radia_export verify_launcher [require \"<label1>,<label2>,...\"]"
  );
  return help;
}

std::vector<std::string> VerifyLauncherCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Headless self-test for the Radia-NGSolve launcher dialog.\n"
    "Lists block + sideset names exactly as the Qt dialog would, then\n"
    "asserts every name in `require` is present. Prints VERIFY_LAUNCHER:\n"
    "key=value lines for parsing by cubit-smoke-test.\n"
    "Use this in CI / deploy gates so a launcher regression cannot ship."
  );
  return help;
}

namespace {

void emit(const std::string& key, const std::string& val)
{
  // PRINT_INFO terminates with newline; structured prefix is parser-friendly.
  PRINT_INFO("VERIFY_LAUNCHER: %s=%s\n", key.c_str(), val.c_str());
}

std::vector<std::string> split_csv(const std::string& s)
{
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
    // trim whitespace
    size_t a = item.find_first_not_of(" \t\r\n");
    size_t b = item.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) continue;
    out.push_back(item.substr(a, b - a + 1));
  }
  return out;
}

}  // namespace

bool VerifyLauncherCommand::execute(CubitCommandData &data)
{
  // 1. List the labels exactly as the launcher dialog sees them.
  auto labels = RadiaLauncherLogic::get_model_labels();

  std::string label_csv;
  for (size_t i = 0; i < labels.size(); ++i) {
    if (i) label_csv += ",";
    label_csv += labels[i];
  }
  emit("labels_count", std::to_string(labels.size()));
  emit("labels", label_csv);

  // 2. Test default_vol_path on the current journal.
  std::string jou = CubitInterface::get_current_journal_file();
  emit("journal", jou);
  std::string voldef = RadiaLauncherLogic::default_vol_path(jou);
  emit("default_vol", voldef);
  if (voldef.empty()) {
    emit("default_vol_status", "FALLBACK_NEEDED");
  } else if (voldef == "/.vol" || voldef == ".vol" || voldef[0] == '/') {
    emit("default_vol_status", "INVALID");
  } else {
    emit("default_vol_status", "OK");
  }

  // 3. Optional: validate against required labels.
  std::string require;
  data.get_string("require", require);
  if (!require.empty()) {
    auto need = split_csv(require);
    auto missing = RadiaLauncherLogic::validate_labels(need, labels);
    std::string missing_csv;
    for (size_t i = 0; i < missing.size(); ++i) {
      if (i) missing_csv += ",";
      missing_csv += missing[i];
    }
    emit("required_count", std::to_string(need.size()));
    emit("missing_count", std::to_string(missing.size()));
    emit("missing", missing_csv);
    if (!missing.empty()) {
      emit("status", "FAIL");
      PRINT_ERROR("verify_launcher: missing required labels: %s\n",
                  missing_csv.c_str());
      return false;
    }
  }

  emit("status", "OK");
  return true;
}
