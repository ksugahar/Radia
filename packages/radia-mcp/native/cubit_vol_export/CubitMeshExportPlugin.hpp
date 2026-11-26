#ifndef CUBIT_MESH_EXPORT_PLUGIN_HPP
#define CUBIT_MESH_EXPORT_PLUGIN_HPP

#include "CubitCommandInterface.hpp"
#include "CubitPluginExport.hpp"

class CubitMeshExportPlugin : public CubitCommandInterface
{
public:
  CubitMeshExportPlugin();
  ~CubitMeshExportPlugin();

  std::vector<std::string> get_keys();
  CubitCommand* create_command(const std::string &key);
};

#endif
