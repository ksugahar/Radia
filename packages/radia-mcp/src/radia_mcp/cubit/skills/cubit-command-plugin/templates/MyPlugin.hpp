#ifndef MY_CUBIT_PLUGIN_HPP
#define MY_CUBIT_PLUGIN_HPP

#include "CubitCommandInterface.hpp"
#include "CubitPluginExport.hpp"

#include <string>
#include <vector>

// One plugin object per binary.  Cubit asks it for command keys, then for
// one CubitCommand per key; the parser owns the returned commands.
class MyPlugin : public CubitCommandInterface
{
public:
  MyPlugin();
  ~MyPlugin() override;

  std::vector<std::string> get_keys() override;
  CubitCommand* create_command(const std::string& key) override;
};

#endif
