#ifndef CUBIT_EXPORT_MODEL_CAPTURE_COMMAND_HPP
#define CUBIT_EXPORT_MODEL_CAPTURE_COMMAND_HPP
#include "CubitCommandInterface.hpp"
#include <string>
#include <vector>
class ExportModelCaptureCommand : public CubitCommand {
public:
  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);
};
#endif
