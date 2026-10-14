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

  bool write_gmsh(const std::string &filename, int order);

private:

  // Gmsh v4.1 writer (all orders, blocks + sidesets + nodesets).
  // v4.1 is the only supported format (lab-wide standard, 2026-04).
  bool write_gmsh_v41(const std::string &filename, const MeshData &mesh);

  // Gmsh element type code
  static int gmsh_type(const MeshElement &elem, int order);

  // Element topological dimension (3=vol, 2=face, 1=edge, 0=point)
  static int elem_dim(const MeshElement &elem);
};

#endif
