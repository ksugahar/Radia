#include "ExportNastranCommand.hpp"
#include "MeshData.hpp"
#include "CubitMessage.hpp"

#include <ctime>
#include <iomanip>
#include <sstream>

ExportNastranCommand::ExportNastranCommand() {}
ExportNastranCommand::~ExportNastranCommand() {}

std::vector<std::string> ExportNastranCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "export nastran <string:label='filename',help='<filename>'> "
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
    "export nastran \"filename\" [order {1|2}] [dimension {2|3}] [nopyramid] [overwrite]"
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
    "  nopyramid     Convert pyramid to degenerate hex (JMAG compatible)\n"
    "  overwrite     Overwrite existing file without warning\n"
  );
  return help;
}

bool ExportNastranCommand::execute(CubitCommandData &data)
{
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

  std::ofstream fid(filename);
  if (!fid.is_open()) {
    PRINT_ERROR("Cannot open file: %s\n", filename.c_str());
    return false;
  }

  bool is_3d = (dim == "3d");

  write_header(fid, filename);
  write_nodes(fid, mesh, is_3d);
  int eid = write_elements(fid, mesh, is_3d, pyram, 1);
  eid = write_sidesets(fid, mesh, is_3d, eid);
  write_nodesets(fid, mesh);
  write_properties(fid, mesh);

  fid << "ENDDATA\n";
  fid.close();

  PRINT_INFO("Exported Nastran BDF (order %d%s): %s (%d nodes, %d elements)\n",
             mesh.order, mesh.has_netgen ? ", NetgenCurver" : "",
             filename.c_str(), mesh.total_node_count(),
             mesh.total_element_count());
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
  fid << "$ Node cards";
  if (mesh.order >= 2)
    fid << " (order " << mesh.order << (mesh.has_netgen ? ", NetgenCurver" : "") << ")";
  fid << "\n$\n";

  for (auto &nd : mesh.nodes) {
    char line1[128], line2[128];
    if (is_3d) {
      std::snprintf(line1, sizeof(line1), "GRID*   %16d%16d%16.5f%16.5f",
                    nd.id, 0, nd.x, nd.y);
      std::snprintf(line2, sizeof(line2), "*       %16.5f", nd.z);
    } else {
      std::snprintf(line1, sizeof(line1), "GRID*   %16d%16d%16.5f%16.5f",
                    nd.id, 0, nd.x, nd.y);
      std::snprintf(line2, sizeof(line2), "*       %16d", 0);
    }
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
  const auto &c = (order >= 2 && !elem.ho_conn.empty()) ? elem.ho_conn : elem.conn;
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

  for (auto &elem : mesh.elements)
    eid = write_element_card(fid, elem, eid, elem.group_id, pyram, mesh.order);

  return eid;
}

// ========================================================================
// Sideset face elements
// ========================================================================
int ExportNastranCommand::write_sidesets(std::ofstream &fid,
                                          const MeshData &mesh,
                                          bool is_3d, int start_eid)
{
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
    fid << "\n";

    for (auto &face : sg.faces)
      eid = write_element_card(fid, face, eid, sg.id, true, mesh.order);
  }
  return eid;
}

// ========================================================================
// Nodesets as comments
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
  }
}

// ========================================================================
// Property cards — PSOLID for blocks, PSHELL for sidesets
// ========================================================================
void ExportNastranCommand::write_properties(std::ofstream &fid, const MeshData &mesh)
{
  fid << "$\n$ Property cards\n$\n";

  for (int bid : mesh.block_ids) {
    char line[128];
    std::snprintf(line, sizeof(line), "PSOLID  %8d%8d", bid, bid);
    fid << line << "\n";
  }

  for (auto &sg : mesh.sidesets) {
    if (sg.faces.empty()) continue;
    char line[128];
    std::snprintf(line, sizeof(line), "PSHELL  %8d%8d", sg.id, sg.id);
    fid << line << "\n";
  }
}
