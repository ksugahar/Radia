#include "ExportVtkCommand.hpp"
#include "MeshData.hpp"
#include "CubitMessage.hpp"
#include "utf8_path.hpp"

#include <fstream>
#include <unordered_map>

// ========================================================================
// Reorder HO mid-edge nodes from internal order to VTK order.
// Derived from GMSH source (GModelIO_VTK.cpp reader permutation).
// ========================================================================
static std::vector<int> reorder_for_vtk(const std::vector<int> &conn,
                                         int nv, ElementType type, int order)
{
  if (order < 2 || (int)conn.size() <= nv)
    return conn;

  const int *reorder = nullptr;
  int n_edge_nodes = 0;

  if (nv == 4 && (type == TETRA4 || type == TETRA)) {
    reorder = EdgeTables::tet_vtk_reorder; n_edge_nodes = 6;
  } else if (nv == 8 && (type == HEX8 || type == HEX)) {
    reorder = EdgeTables::hex_vtk_reorder; n_edge_nodes = 12;
  } else if (nv == 6 && (type == WEDGE6 || type == WEDGE)) {
    reorder = EdgeTables::wedge_vtk_reorder; n_edge_nodes = 9;
  } else if (nv == 5 && (type == PYRAMID5 || type == PYRAMID)) {
    reorder = EdgeTables::pyramid_vtk_reorder; n_edge_nodes = 8;
  } else if (nv == 3 && (type == TRI3 || type == CUBIT_TRI
                       || type == TRISHELL || type == TRISHELL3)) {
    reorder = EdgeTables::tri_vtk_reorder; n_edge_nodes = 3;
  } else if (nv == 4 && (type == QUAD4 || type == QUAD
                       || type == SHEL || type == SHELL4)) {
    reorder = EdgeTables::quad_vtk_reorder; n_edge_nodes = 4;
  }

  if (!reorder) return conn;

  int n_mid = (int)conn.size() - nv;
  if (n_mid < n_edge_nodes) return conn;

  std::vector<int> out(conn.size());
  for (int i = 0; i < nv; i++)
    out[i] = conn[i];
  for (int i = 0; i < n_edge_nodes; i++)
    out[nv + reorder[i]] = conn[nv + i];
  for (int i = nv + n_edge_nodes; i < (int)conn.size(); i++)
    out[i] = conn[i];
  return out;
}

ExportVtkCommand::ExportVtkCommand() {}
ExportVtkCommand::~ExportVtkCommand() {}

std::vector<std::string> ExportVtkCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "radia_export vtk <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1 or 2>'>] "
    "[dimension <value:label='dimension',help='<2 or 3>'>] "
    "[overwrite]"
  );
  return syntax_list;
}

std::vector<std::string> ExportVtkCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "radia_export vtk \"filename\" [order {1|2}] [dimension {2|3}] [overwrite]"
  );
  return help;
}

std::vector<std::string> ExportVtkCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Export mesh to VTK Legacy format (.vtk).\n"
    "Block assignment is NOT required.\n"
    "Sidesets exported as surface cells.\n\n"
    "Options:\n"
    "  order 1      1st-order elements (default)\n"
    "  order 2      2nd-order elements with geometry projection\n"
    "  dimension 3  3D solid mesh output (default)\n"
    "  dimension 2  2D shell mesh output\n"
    "  overwrite    Overwrite existing file\n"
  );
  return help;
}

bool ExportVtkCommand::execute(CubitCommandData &data)
{
  std::string filename;
  data.get_string("filename", filename);

  int order = 1;
  data.get_value("order", order);
  if (order < 1) order = 1;
  if (order > 2) {
    PRINT_WARNING("VTK order %d not supported, using order 2.\n", order);
    order = 2;
  }

  int dim_val = 3;
  data.get_value("dimension", dim_val);
  std::string dim = (dim_val == 2) ? "2d" : "3d";

  return write_vtk(filename, dim, order);
}

// ========================================================================
// VTK cell type mapping
// ========================================================================
int ExportVtkCommand::vtk_type(const MeshElement &elem, int order)
{
  int nv = elem.nv;
  bool is_tet  = (nv == 4 && (elem.type == TETRA4 || elem.type == TETRA));
  bool is_hex  = (nv == 8 && (elem.type == HEX8   || elem.type == HEX));
  bool is_wed  = (nv == 6 && (elem.type == WEDGE6 || elem.type == WEDGE));
  bool is_pyr  = (nv == 5 && (elem.type == PYRAMID5 || elem.type == PYRAMID));
  bool is_tri  = (nv == 3 && (elem.type == TRI3   || elem.type == CUBIT_TRI
                         || elem.type == TRISHELL || elem.type == TRISHELL3));
  bool is_quad = (nv == 4 && (elem.type == QUAD4  || elem.type == QUAD
                         || elem.type == SHEL  || elem.type == SHELL4));
  bool is_line = (elem.type == BAR || elem.type == BAR2 || elem.type == BAR3);
  bool is_pt   = (nv == 1);

  if (order >= 2) {
    if (is_tet)  return 24;  // VTK_QUADRATIC_TETRA
    if (is_hex)  return 25;  // VTK_QUADRATIC_HEXAHEDRON
    if (is_wed)  return 26;  // VTK_QUADRATIC_WEDGE
    if (is_pyr)  return 27;  // VTK_QUADRATIC_PYRAMID
    if (is_tri)  return 22;  // VTK_QUADRATIC_TRIANGLE
    if (is_quad) return 23;  // VTK_QUADRATIC_QUAD
  }

  if (is_tet)  return 10;
  if (is_hex)  return 12;
  if (is_wed)  return 13;
  if (is_pyr)  return 14;
  if (is_tri)  return 5;
  if (is_quad) return 9;
  if (is_line) return 3;
  if (is_pt)   return 1;
  return 0;
}

// ========================================================================
// Main VTK writer
// ========================================================================
bool ExportVtkCommand::write_vtk(const std::string &filename,
                                  const std::string &dim, int order)
{
  MeshData mesh;
  if (!mesh.extract(order))
    return false;

  bool is_3d = (dim == "3d");

  // VTK needs 0-based indices — use mesh.node_id_to_index
  auto idx = [&](int nid) -> int { return mesh.node_id_to_index.at(nid); };

  // Collect all elements (blocks + sideset faces)
  struct CellData {
    int vtk_type;
    std::vector<int> conn_idx;  // 0-based
    int group_id;
    int group_type;  // 0=block, 1=sideset
  };
  std::vector<CellData> cells;

  auto add_elem = [&](const MeshElement &elem, int gid, int gtype) {
    int vt = vtk_type(elem, order);
    if (vt == 0) return;
    const auto &raw = (order >= 2 && !elem.ho_conn.empty()) ? elem.ho_conn : elem.conn;
    const auto c = reorder_for_vtk(raw, elem.nv, elem.type, order);
    CellData cd;
    cd.vtk_type = vt;
    cd.group_id = gid;
    cd.group_type = gtype;
    cd.conn_idx.reserve(c.size());
    for (int nid : c) cd.conn_idx.push_back(idx(nid));
    cells.push_back(std::move(cd));
  };

  for (auto &elem : mesh.elements)
    add_elem(elem, elem.group_id, 0);
  for (auto &sg : mesh.sidesets)
    for (auto &face : sg.faces)
      add_elem(face, sg.id, 1);

  // --- Write VTK file ---
  std::ofstream fid(u8_string_to_path(filename));
  if (!fid.is_open()) {
    PRINT_ERROR("Cannot open file: %s\n", filename.c_str());
    return false;
  }

  int num_nodes = mesh.total_node_count();
  int num_cells = (int)cells.size();

  fid << "# vtk DataFile Version 3.0\n";
  fid << "Radia Cubit Plugin";
  if (order >= 2) fid << " (order " << order << ")";
  fid << "\n";
  fid << "ASCII\n";
  fid << "DATASET UNSTRUCTURED_GRID\n";

  // Points
  fid << "POINTS " << num_nodes << " double\n";
  for (auto &nd : mesh.nodes) {
    if (is_3d)
      fid << nd.x << " " << nd.y << " " << nd.z << "\n";
    else
      fid << nd.x << " " << nd.y << " 0\n";
  }

  // Cells
  int cell_list_size = 0;
  for (auto &cd : cells) cell_list_size += 1 + (int)cd.conn_idx.size();

  fid << "CELLS " << num_cells << " " << cell_list_size << "\n";
  for (auto &cd : cells) {
    fid << cd.conn_idx.size();
    for (int i : cd.conn_idx) fid << " " << i;
    fid << "\n";
  }

  // Cell types
  fid << "CELL_TYPES " << num_cells << "\n";
  for (auto &cd : cells) fid << cd.vtk_type << "\n";

  // Cell data: GroupID and GroupType
  fid << "CELL_DATA " << num_cells << "\n";
  fid << "SCALARS GroupID int 1\n";
  fid << "LOOKUP_TABLE default\n";
  for (auto &cd : cells) fid << cd.group_id << "\n";
  fid << "SCALARS GroupType int 1\n";
  fid << "LOOKUP_TABLE default\n";
  for (auto &cd : cells) fid << cd.group_type << "\n";

  fid.close();
  PRINT_INFO("Exported VTK (order %d): %s (%d nodes, %d cells)\n",
             order, filename.c_str(), num_nodes, num_cells);
  return true;
}
