#ifndef EXPORT_GMSH_COMMAND_HPP
#define EXPORT_GMSH_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include <fstream>
#include <string>
#include <vector>

class MeshData;
struct MeshElement;

class ExportGmshCommand : public CubitCommand
{
public:
  ExportGmshCommand();
  ~ExportGmshCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);

  bool write_gmsh(const std::string &filename, const std::string &version,
                  int order);

private:

  // Gmsh v2.2 writer (all orders, blocks + sidesets + nodesets)
  bool write_gmsh_v22(const std::string &filename, const MeshData &mesh);

  // Gmsh v4.1 writer (all orders, blocks + sidesets + nodesets)
  bool write_gmsh_v41(const std::string &filename, const MeshData &mesh);

  // Gmsh element type code
  static int gmsh_type(const MeshElement &elem, int order);

  // Element topological dimension (3=vol, 2=face, 1=edge, 0=point)
  static int elem_dim(const MeshElement &elem);
};

#endif
