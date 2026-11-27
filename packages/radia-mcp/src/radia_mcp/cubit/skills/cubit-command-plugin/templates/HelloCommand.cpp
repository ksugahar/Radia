#include "HelloCommand.hpp"
#include "CubitMessage.hpp"

std::vector<std::string> HelloCommand::get_syntax()
{
  // Modified BNF: <type:label=...> required, [..] optional (any order),
  // {a|b} choice.  The label is the key used in execute().
  return {
    "myplugin hello <string:label='message',help='<message>'> "
    "[count <value:label='count',help='<count>'>]"
  };
}

std::vector<std::string> HelloCommand::get_syntax_help()
{
  return {"myplugin hello \"message\" [count <n>]"};
}

std::vector<std::string> HelloCommand::get_help()
{
  return {"Print a message; used to prove the plugin is loaded.\n"};
}

bool HelloCommand::execute(CubitCommandData& data)
{
  // The command object persists between calls: keep state local.
  std::string message;
  data.get_string("message", message);
  int count = 1;
  data.get_value("count", count);  // false (and count unchanged) if absent
  if (count < 1 || count > 1000) {
    PRINT_ERROR("count must be between 1 and 1000, got %d\n", count);
    return false;  // false = command failed; Cubit counts an error
  }
  for (int i = 0; i < count; ++i)
    PRINT_INFO("MYPLUGIN_HELLO %s\n", message.c_str());
  return true;
}
