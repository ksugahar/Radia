#ifndef VERIFY_LAUNCHER_COMMAND_HPP
#define VERIFY_LAUNCHER_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include <string>
#include <vector>

/// APREPRO command: radia_export verify_launcher
///
/// Headless self-test for the Radia-NGSolve launcher dialog. Calls the
/// SAME pure-logic functions (RadiaLauncherLogic::get_model_labels,
/// default_vol_path, validate_labels) that the .ccl Qt dialog uses, so
/// running this from `cubit-smoke-test` catches the same class of bug
/// (e.g. 2026-04-14 sideset name lookup regression) that the dialog
/// would otherwise expose only at user GUI time.
///
/// Output is a structured "VERIFY_LAUNCHER:" key=value list to stdout
/// for easy parsing by the python smoke test wrapper.
class VerifyLauncherCommand : public CubitCommand
{
public:
  VerifyLauncherCommand();
  ~VerifyLauncherCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);
};

#endif
