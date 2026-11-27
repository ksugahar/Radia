#ifndef HELLO_COMMAND_HPP
#define HELLO_COMMAND_HPP

#include "CubitCommandInterface.hpp"

#include <string>
#include <vector>

// Smallest useful command: proves the plugin loaded and the parser routes
// the verb.  Rename the verb before shipping; a verb that matches a built-in
// command masks it, because plugin commands are parsed first.
class HelloCommand : public CubitCommand
{
public:
  std::vector<std::string> get_syntax() override;
  std::vector<std::string> get_syntax_help() override;
  std::vector<std::string> get_help() override;
  bool execute(CubitCommandData& data) override;
};

#endif
