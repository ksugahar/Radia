#ifndef CUBIT_TO_VOL_EXPORT_PLUGIN_HPP
#define CUBIT_TO_VOL_EXPORT_PLUGIN_HPP

#include "CubitCommandInterface.hpp"
#include "CubitPluginExport.hpp"

class CubitToVolExportPlugin : public CubitCommandInterface
{
public:
  CubitToVolExportPlugin();
  ~CubitToVolExportPlugin();

  std::vector<std::string> get_keys();
  CubitCommand* create_command(const std::string &key);
};

#endif
