#include "ExportNastranCommand.hpp"
#include "MeshData.hpp"
#include "RadiaMessageFilter.hpp"
#include "CubitMessage.hpp"
#include "utf8_path.hpp"

#include <algorithm>
#include <ctime>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

// ========================================================================
// Reorder HO mid-edge nodes from internal (VTK) order to Nastran BDF order.
// Nastran CHEXA: bottom + vertical + top (internal: bottom + top + vertical)
// Nastran CPENTA: bottom + vertical + top (internal: bottom + top + vertical)
// ========================================================================
static std::vector<int> reorder_for_bdf(const std::vector<int> &conn,
                                         int nv, ElementType type, int order)
{
  if (order < 2 || (int)conn.size() <= nv)
    return conn;

  const int *reorder = nullptr;
  int n_edge_nodes = 0;

  if (nv == 4 && (type == TETRA4 || type == TETRA)) {
    reorder = EdgeTables::tet_bdf_reorder; n_edge_nodes = 6;
  } else if (nv == 8 && (type == HEX8 || type == HEX)) {
    reorder = EdgeTables::hex_bdf_reorder; n_edge_nodes = 12;
  } else if (nv == 6 && (type == WEDGE6 || type == WEDGE)) {
    reorder = EdgeTables::wedge_bdf_reorder; n_edge_nodes = 9;
  } else if (nv == 5 && (type == PYRAMID5 || type == PYRAMID)) {
    reorder = EdgeTables::pyramid_bdf_reorder; n_edge_nodes = 8;
  } else if (nv == 3 && (type == TRI3 || type == CUBIT_TRI
                       || type == TRISHELL || type == TRISHELL3)) {
    reorder = EdgeTables::tri_bdf_reorder; n_edge_nodes = 3;
  } else if (nv == 4 && (type == QUAD4 || type == QUAD
                       || type == SHEL || type == SHELL4)) {
    reorder = EdgeTables::quad_bdf_reorder; n_edge_nodes = 4;
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

static int bdf_element_dimension(const MeshElement &elem)
{
  if (elem.nv == 4 && (elem.type == TETRA4 || elem.type == TETRA)) return 3;
  if (elem.nv == 8 && (elem.type == HEX8 || elem.type == HEX)) return 3;
  if (elem.nv == 6 && (elem.type == WEDGE6 || elem.type == WEDGE)) return 3;
  if (elem.nv == 5 && (elem.type == PYRAMID5 || elem.type == PYRAMID)) return 3;
  if (elem.nv == 3 && (elem.type == TRI3 || elem.type == CUBIT_TRI
                    || elem.type == TRISHELL || elem.type == TRISHELL3)) return 2;
  if (elem.nv == 4 && (elem.type == QUAD4 || elem.type == QUAD
                    || elem.type == SHEL || elem.type == SHELL4)) return 2;
  if (elem.type == BAR || elem.type == BAR2 || elem.type == BAR3) return 1;
  return 0;
}

static bool include_bdf_element(const MeshElement &elem, bool is_3d)
{
  const int dim = bdf_element_dimension(elem);
  return dim >= 2 && (is_3d || dim == 2);
}

static std::map<int, int> assign_sideset_property_ids(const MeshData &mesh)
{
  std::set<int> used(mesh.block_ids.begin(), mesh.block_ids.end());
  int next_pid = used.empty() ? 1 : (*used.rbegin() + 1);
  std::map<int, int> result;
  for (const auto &sg : mesh.sidesets) {
    if (sg.faces.empty()) continue;
    int pid = sg.id;
    if (used.count(pid)) {
      while (used.count(next_pid)) ++next_pid;
      pid = next_pid++;
    }
    used.insert(pid);
    result[sg.id] = pid;
  }
  return result;
}

ExportNastranCommand::ExportNastranCommand() {}
ExportNastranCommand::~ExportNastranCommand() {}

std::vector<std::string> ExportNastranCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "export nastran_bdf <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1 or 2>'>] "
    "[dimension <value:label='dimension',help='<2 or 3>'>] "
    "[nopyramid] "
    "[overwrite]"
  );
  return syntax_list;
}

std::vector<std::string> ExportNastranCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "export nastran_bdf \"filename\" [order {1|2}] [dimension {2|3}] [nopyramid] [overwrite]"
  );
  return help;
}

std::vector<std::string> ExportNastranCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Export mesh to NX Nastran BDF format.\n"
    "Block assignment is NOT required - all meshed elements are exported.\n"
    "Sidesets are exported as surface elements. Nodesets as comments.\n\n"
    "Options:\n"
    "  order 1       1st-order elements (default)\n"
    "  order 2       2nd-order elements with geometry projection\n"
    "  dimension 3   3D solid mesh output (default)\n"
    "  dimension 2   2D shell mesh output\n"
    "  nopyramid     Convert pyramid to degenerate hex (solver-compatible)\n"
    "  overwrite     Overwrite existing file without warning\n"
  );
  return help;
}

std::vector<std::string> ExportJmagNastranCommand::get_syntax()
{
  return {
    "export jmag_nastran <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1 or 2>'>] "
    "[dimension <value:label='dimension',help='<2 or 3>'>] "
    "[nopyramid] "
    "[overwrite]"
  };
}

std::vector<std::string> ExportJmagNastranCommand::get_syntax_help()
{
  return {
    "export jmag_nastran \"filename\" [order {1|2}] "
    "[dimension {2|3}] [nopyramid] [overwrite]"
  };
}

std::vector<std::string> ExportJmagNastranCommand::get_help()
{
  auto help = ExportNastranCommand::get_help();
  help.push_back(
    "Deprecated compatibility alias. Prefer export nastran_bdf."
  );
  return help;
}

bool ExportNastranCommand::execute(CubitCommandData &data)
{
  // Suppress Cubit Learn Edition's harmless 50k-cap ERROR.
  radia::ScopedLearnEditionFilter _lef_guard;

  std::string filename;
  data.get_string("filename", filename);

  int order = 1;
  data.get_value("order", order);
  if (order < 1) order = 1;
  if (order > 2) {
    PRINT_WARNING("Nastran BDF supports order <= 2. Using order 2.\n");
    order = 2;
  }

  int dim_val = 3;
  data.get_value("dimension", dim_val);
  std::string dim = (dim_val == 2) ? "2d" : "3d";

  bool pyram = !data.find_keyword("nopyramid");

  return write_nastran(filename, dim, pyram, order);
}

// ========================================================================
// Main writer
// ========================================================================
bool ExportNastranCommand::write_nastran(const std::string &filename,
                                          const std::string &dim,
                                          bool pyram, int order)
{
  MeshData mesh;
  if (!mesh.extract(order))
    return false;

  const bool is_3d = (dim == "3d");
  int selected_elements = 0;
  std::map<int, int> block_dimensions;
  for (const auto &elem : mesh.elements) {
    if (!include_bdf_element(elem, is_3d)) continue;
    ++selected_elements;
    const int elem_dim = bdf_element_dimension(elem);
    auto inserted = block_dimensions.emplace(elem.group_id, elem_dim);
    if (!inserted.second && inserted.first->second != elem_dim) {
      PRINT_ERROR("Nastran block %d mixes element dimensions %d and %d. "
                  "Split it into separate Cubit blocks before export.\n",
                  elem.group_id, inserted.first->second, elem_dim);
      return false;
    }
  }
  for (const auto &sg : mesh.sidesets)
    selected_elements += (int)sg.faces.size();
  if (selected_elements == 0) {
    PRINT_ERROR("Nastran dimension %d selected no compatible elements.\n",
                is_3d ? 3 : 2);
    return false;
  }

  const auto sideset_pids = assign_sideset_property_ids(mesh);

  std::ofstream fid(u8_string_to_path(filename));
  if (!fid.is_open()) {
    PRINT_ERROR("Cannot open file: %s\n", filename.c_str());
    return false;
  }

  write_header(fid, filename);
  write_nodes(fid, mesh, is_3d);
  int eid = write_elements(fid, mesh, is_3d, pyram, 1);
  eid = write_sidesets(fid, mesh, is_3d, eid, sideset_pids);
  write_nodesets(fid, mesh);
  write_properties(fid, mesh, is_3d, sideset_pids);

  fid << "ENDDATA\n";
  fid.close();

  PRINT_INFO("Exported Nastran BDF (order %d%s): %s (%d nodes, %d selected elements)\n",
             mesh.order, mesh.has_netgen ? ", NetgenCurver" : "",
             filename.c_str(), mesh.total_node_count(),
             selected_elements);
  return true;
}

// ========================================================================
// Header
// ========================================================================
void ExportNastranCommand::write_header(std::ofstream &fid, const std::string &filename)
{
  auto t = std::time(nullptr);
  struct tm tm_buf;
#ifdef _MSC_VER
  localtime_s(&tm_buf, &t);
#else
  tm_buf = *std::localtime(&t);
#endif
  std::ostringstream ts;
  ts << std::put_time(&tm_buf, "%d-%b-%y at %H:%M:%S");

  fid << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n";
  fid << "$\n";
  fid << "$                    Radia Cubit Plugin - Nastran Exporter\n";
  fid << "$\n";
  fid << "$   File: " << filename << "\n";
  fid << "$   Time: " << ts.str() << "\n";
  fid << "$\n";
  fid << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n";
  fid << "$\n";
  fid << "SOL 101\n";
  fid << "CEND\n";
  fid << "$\n";
  fid << "ECHO = SORT\n";
  fid << "SUBCASE = 1\n";
  fid << "LABEL = Default Set\n";
  fid << "$\n";
  fid << "BEGIN BULK\n";
  fid << "$\n";
}

// ========================================================================
// Nodes — original + high-order
// ========================================================================
void ExportNastranCommand::write_nodes(std::ofstream &fid, const MeshData &mesh,
                                        bool is_3d)
{
  (void)is_3d;
  fid << "$ Node cards";
  if (mesh.order >= 2)
    fid << " (order " << mesh.order << (mesh.has_netgen ? ", NetgenCurver" : "") << ")";
  fid << "\n$\n";

  for (auto &nd : mesh.nodes) {
    char line1[128], line2[128];
    std::snprintf(line1, sizeof(line1), "GRID*   %16d%16d%16.9e%16.9e",
                  nd.id, 0, nd.x, nd.y);
    std::snprintf(line2, sizeof(line2), "*       %16.9e", nd.z);
    fid << line1 << "\n" << line2 << "\n";
  }
}

// ========================================================================
// Single element card writer (handles both linear and HO)
// ========================================================================
int ExportNastranCommand::write_element_card(std::ofstream &fid,
                                              const MeshElement &elem,
                                              int eid, int pid,
                                              bool pyram, int order)
{
  char line[256];
  const auto &raw = (order >= 2 && !elem.ho_conn.empty()) ? elem.ho_conn : elem.conn;
  const auto c = reorder_for_bdf(raw, elem.nv, elem.type, order);
  int nv = elem.nv;
  int nn = (int)c.size();

  // --- 3D volume elements ---
  if (nv == 4 && (elem.type == TETRA4 || elem.type == TETRA)) {
    if (nn >= 10) {  // TET10
      std::snprintf(line, sizeof(line),
        "CTETRA  %8d%8d%8d%8d%8d%8d%8d%8d+",
        eid, pid, c[0], c[1], c[2], c[3], c[4], c[5]);
      fid << line << "\n";
      std::snprintf(line, sizeof(line),
        "+       %8d%8d%8d%8d", c[6], c[7], c[8], c[9]);
      fid << line << "\n";
    } else {  // TET4
      std::snprintf(line, sizeof(line),
        "CTETRA  %8d%8d%8d%8d%8d%8d",
        eid, pid, c[0], c[1], c[2], c[3]);
      fid << line << "\n";
    }
  }
  else if (nv == 8 && (elem.type == HEX8 || elem.type == HEX)) {
    if (nn >= 20) {  // HEX20
      std::snprintf(line, sizeof(line),
        "CHEXA   %8d%8d%8d%8d%8d%8d%8d%8d+",
        eid, pid, c[0], c[1], c[2], c[3], c[4], c[5]);
      fid << line << "\n";
      std::snprintf(line, sizeof(line),
        "+       %8d%8d%8d%8d%8d%8d%8d%8d+",
        c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13]);
      fid << line << "\n";
      std::snprintf(line, sizeof(line),
        "+       %8d%8d%8d%8d%8d%8d",
        c[14], c[15], c[16], c[17], c[18], c[19]);
      fid << line << "\n";
    } else {  // HEX8
      std::snprintf(line, sizeof(line),
        "CHEXA   %8d%8d%8d%8d%8d%8d%8d%8d+",
        eid, pid, c[0], c[1], c[2], c[3], c[4], c[5]);
      fid << line << "\n";
      std::snprintf(line, sizeof(line),
        "+       %8d%8d", c[6], c[7]);
      fid << line << "\n";
    }
  }
  else if (nv == 6 && (elem.type == WEDGE6 || elem.type == WEDGE)) {
    if (nn >= 15) {  // WEDGE15
      std::snprintf(line, sizeof(line),
        "CPENTA  %8d%8d%8d%8d%8d%8d%8d%8d+",
        eid, pid, c[0], c[1], c[2], c[3], c[4], c[5]);
      fid << line << "\n";
      std::snprintf(line, sizeof(line),
        "+       %8d%8d%8d%8d%8d%8d%8d%8d+",
        c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13]);
      fid << line << "\n";
      std::snprintf(line, sizeof(line), "+       %8d", c[14]);
      fid << line << "\n";
    } else {  // WEDGE6
      std::snprintf(line, sizeof(line),
        "CPENTA  %8d%8d%8d%8d%8d%8d%8d%8d",
        eid, pid, c[0], c[1], c[2], c[3], c[4], c[5]);
      fid << line << "\n";
    }
  }
  else if (nv == 5 && (elem.type == PYRAMID5 || elem.type == PYRAMID)) {
    if (pyram) {
      if (nn >= 13) {  // PYRAMID13
        std::snprintf(line, sizeof(line),
          "CPYRAM  %8d%8d%8d%8d%8d%8d%8d%8d+",
          eid, pid, c[0], c[1], c[2], c[3], c[4], c[5]);
        fid << line << "\n";
        std::snprintf(line, sizeof(line),
          "+       %8d%8d%8d%8d%8d%8d%8d",
          c[6], c[7], c[8], c[9], c[10], c[11], c[12]);
        fid << line << "\n";
      } else {  // PYRAMID5
        std::snprintf(line, sizeof(line),
          "CPYRAM  %8d%8d%8d%8d%8d%8d%8d",
          eid, pid, c[0], c[1], c[2], c[3], c[4]);
        fid << line << "\n";
      }
    } else {
      // Degenerate hex for JMAG compatibility
      std::snprintf(line, sizeof(line),
        "CHEXA   %8d%8d%8d%8d%8d%8d%8d%8d+",
        eid, pid, c[0], c[1], c[2], c[3], c[4], c[4]);
      fid << line << "\n";
      std::snprintf(line, sizeof(line), "+       %8d%8d", c[4], c[4]);
      fid << line << "\n";
    }
  }
  // --- 2D surface elements ---
  else if (nv == 3 && (elem.type == TRI3 || elem.type == CUBIT_TRI
                  || elem.type == TRISHELL || elem.type == TRISHELL3)) {
    if (nn >= 6) {  // TRI6
      std::snprintf(line, sizeof(line),
        "CTRIA6  %8d%8d%8d%8d%8d%8d%8d%8d",
        eid, pid, c[0], c[1], c[2], c[3], c[4], c[5]);
      fid << line << "\n";
    } else {  // TRI3
      std::snprintf(line, sizeof(line),
        "CTRIA3  %8d%8d%8d%8d%8d",
        eid, pid, c[0], c[1], c[2]);
      fid << line << "\n";
    }
  }
  else if (nv == 4 && (elem.type == QUAD4 || elem.type == QUAD
                  || elem.type == SHEL || elem.type == SHELL4)) {
    if (nn >= 8) {  // QUAD8
      std::snprintf(line, sizeof(line),
        "CQUAD8  %8d%8d%8d%8d%8d%8d%8d%8d+",
        eid, pid, c[0], c[1], c[2], c[3], c[4], c[5]);
      fid << line << "\n";
      std::snprintf(line, sizeof(line), "+       %8d%8d", c[6], c[7]);
      fid << line << "\n";
    } else {  // QUAD4
      std::snprintf(line, sizeof(line),
        "CQUAD4  %8d%8d%8d%8d%8d%8d",
        eid, pid, c[0], c[1], c[2], c[3]);
      fid << line << "\n";
    }
  }
  // --- 1D elements ---
  else if (elem.type == BAR || elem.type == BAR2 || elem.type == BAR3) {
    std::snprintf(line, sizeof(line),
      "CROD    %8d%8d%8d%8d", eid, pid, c[0], c[1]);
    fid << line << "\n";
  }

  return eid + 1;
}

// ========================================================================
// Block elements
// ========================================================================
int ExportNastranCommand::write_elements(std::ofstream &fid,
                                          const MeshData &mesh,
                                          bool is_3d, bool pyram,
                                          int start_eid)
{
  fid << "$\n$ Element cards\n$\n";
  int eid = start_eid;

  for (auto &elem : mesh.elements) {
    if (!include_bdf_element(elem, is_3d)) continue;
    eid = write_element_card(fid, elem, eid, elem.group_id, pyram, mesh.order);
  }

  return eid;
}

// ========================================================================
// Sideset face elements
// ========================================================================
int ExportNastranCommand::write_sidesets(std::ofstream &fid,
                                          const MeshData &mesh,
                                          bool is_3d, int start_eid,
                                          const std::map<int, int> &sideset_pids)
{
  (void)is_3d;
  int eid = start_eid;
  bool has_any = false;

  for (auto &sg : mesh.sidesets) {
    if (sg.faces.empty()) continue;
    if (!has_any) {
      fid << "$\n$ Sideset face elements\n$\n";
      has_any = true;
    }
    fid << "$ Sideset " << sg.id;
    if (!sg.name.empty()) fid << " (" << sg.name << ")";
    auto pid_it = sideset_pids.find(sg.id);
    const int pid = pid_it == sideset_pids.end() ? sg.id : pid_it->second;
    fid << ", property " << pid;
    fid << "\n";

    for (auto &face : sg.faces)
      eid = write_element_card(fid, face, eid, pid, true, mesh.order);
  }
  return eid;
}

// ========================================================================
// Nodesets as comments plus machine-readable SET1 cards
// ========================================================================
void ExportNastranCommand::write_nodesets(std::ofstream &fid, const MeshData &mesh)
{
  for (auto &ng : mesh.nodesets) {
    if (ng.node_ids.empty()) continue;
    fid << "$\n$ Nodeset " << ng.id;
    if (!ng.name.empty()) fid << " (" << ng.name << ")";
    fid << "\n$  ";
    for (int i = 0; i < (int)ng.node_ids.size(); i++) {
      if (i > 0 && i % 10 == 0) fid << "\n$  ";
      fid << " " << ng.node_ids[i];
    }
    fid << "\n";

    // Small fixed-field form: one set ID plus seven node IDs on the first
    // line, then eight node IDs on each continuation line.  A previous
    // free-field "+,..." continuation inserted blank fields in strict
    // readers such as pyNastran.
    fid << std::left << std::setw(8) << "SET1" << std::right
        << std::setw(8) << ng.id;
    int fields_used = 1;
    for (int node_id : ng.node_ids) {
      if (fields_used == 8) {
        fid << "\n" << std::setw(8) << "";
        fields_used = 0;
      }
      fid << std::setw(8) << node_id;
      ++fields_used;
    }
    fid << "\n";
  }
}

// ========================================================================
// Property cards — PSOLID for blocks, PSHELL for sidesets
// ========================================================================
void ExportNastranCommand::write_properties(
    std::ofstream &fid, const MeshData &mesh, bool is_3d,
    const std::map<int, int> &sideset_pids)
{
  fid << "$\n$ Property cards\n$\n";

  std::map<int, int> block_dimensions;
  for (const auto &elem : mesh.elements) {
    if (!include_bdf_element(elem, is_3d)) continue;
    block_dimensions[elem.group_id] = bdf_element_dimension(elem);
  }

  for (const auto &entry : block_dimensions) {
    const int bid = entry.first;
    char line[128];
    if (entry.second == 3)
      std::snprintf(line, sizeof(line), "PSOLID  %8d%8d", bid, bid);
    else
      std::snprintf(line, sizeof(line), "PSHELL  %8d%8d", bid, bid);
    fid << line << "\n";
  }

  for (auto &sg : mesh.sidesets) {
    if (sg.faces.empty()) continue;
    auto pid_it = sideset_pids.find(sg.id);
    const int pid = pid_it == sideset_pids.end() ? sg.id : pid_it->second;
    char line[128];
    std::snprintf(line, sizeof(line), "PSHELL  %8d%8d", pid, sg.id);
    fid << line << "\n";
  }
}
