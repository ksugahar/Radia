#ifndef EXPORT_NETGEN_COMMAND_HPP
#define EXPORT_NETGEN_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include <string>
#include <vector>

/// APREPRO command: export netgen "file.vol" [order N]
/// Runs Python subprocess (calc_export_netgen.py) for high-order curving.
/// Produces .vol + companion .json with CAD reference values.
class ExportNetgenCommand : public CubitCommand
{
public:
  ExportNetgenCommand();
  ~ExportNetgenCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);
};

#endif
