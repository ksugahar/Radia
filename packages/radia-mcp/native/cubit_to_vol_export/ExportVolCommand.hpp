#ifndef EXPORT_VOL_COMMAND_HPP
#define EXPORT_VOL_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include <string>
#include <vector>

/// APREPRO command: export vol "file.vol" [order N]
/// Pure C++ export via MeshData + NetgenCurver (no subprocess).
/// Produces .vol + companion .json with CAD reference values.
class ExportVolCommand : public CubitCommand
{
public:
  ExportVolCommand();
  ~ExportVolCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);
};

#endif
