#ifndef RADIA_PLUGIN_HPP
#define RADIA_PLUGIN_HPP

#include "CubitCommandInterface.hpp"
#include "CubitPluginExport.hpp"

class RadiaPlugin : public CubitCommandInterface
{
public:
  RadiaPlugin();
  ~RadiaPlugin();

  std::vector<std::string> get_keys();
  CubitCommand* create_command(const std::string &key);
};

#endif
