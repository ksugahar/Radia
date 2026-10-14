#ifndef EXPORT_VTK_COMMAND_HPP
#define EXPORT_VTK_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include <fstream>
#include <string>
#include <vector>

class MeshData;
struct MeshElement;

class ExportVtkCommand : public CubitCommand
{
public:
  ExportVtkCommand();
  ~ExportVtkCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);

  bool write_vtk(const std::string &filename, const std::string &dim, int order);

private:

  // VTK cell type from element
  static int vtk_type(const MeshElement &elem, int order);
};

#endif
