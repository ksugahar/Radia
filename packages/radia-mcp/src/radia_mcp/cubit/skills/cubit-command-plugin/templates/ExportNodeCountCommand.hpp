#ifndef EXPORT_NODE_COUNT_COMMAND_HPP
#define EXPORT_NODE_COUNT_COMMAND_HPP

#include "CubitCommandInterface.hpp"

#include <string>
#include <vector>

// Template for a mesh writer that extends Cubit's own `export` verb.
// Replace the body of write() with the real format.
class ExportNodeCountCommand : public CubitCommand
{
public:
  std::vector<std::string> get_syntax() override;
  std::vector<std::string> get_syntax_help() override;
  std::vector<std::string> get_help() override;
  bool execute(CubitCommandData& data) override;
};

#endif
