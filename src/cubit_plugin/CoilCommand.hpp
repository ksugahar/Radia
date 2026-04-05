#ifndef COIL_COMMAND_HPP
#define COIL_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include <string>
#include <vector>

/// APREPRO command: coil "script.py" [output "path.step"] [noimport]
///
/// Generates a coil STEP file via external Python (CoilBuilder) and
/// optionally imports it into Cubit.
///
/// The script must define build_coil() -> CoilBuilder.
/// Requires Python 3.12 with NGSolve/OCC (not Cubit's embedded Python).
class CoilCommand : public CubitCommand
{
public:
  CoilCommand();
  ~CoilCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);
};

#endif // COIL_COMMAND_HPP
