#include "ExportFemeemCommand.hpp"
#include "MeshData.hpp"
#include "RadiaMessageFilter.hpp"
#include "CubitMessage.hpp"
#include "utf8_path.hpp"

#include <fstream>
#include <iomanip>
#include <set>
#include <map>
#include <algorithm>
#include <cstdio>
#include <sys/stat.h>
#include <direct.h>   // _mkdir on Windows

ExportFemeemCommand::ExportFemeemCommand() {}
ExportFemeemCommand::~ExportFemeemCommand() {}

std::vector<std::string> ExportFemeemCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "export femeem <string:label='filename',help='<directory>'> "
    "[scale <value:label='scale',help='<coordinate scale>'>] "
    "[overwrite]"
  );
  return syntax_list;
}

std::vector<std::string> ExportFemeemCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "export femeem \"dirname\" [scale <value>] [overwrite]"
  );
  return help;
}

std::vector<std::string> ExportFemeemCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Export mesh to FEMEEM format (Gifu Univ. 3D FEM solver).\n"
    "Creates a directory with in.dat, sin.dat.B, sina.dat, and d3.\n"
    "Only 1st-order tetrahedral elements are exported.\n"
    "Non-tet elements (hex, wedge, pyramid) are skipped with a warning.\n\n"
    "Options:\n"
    "  scale <value>  Coordinate scale factor (default 1.0)\n"
    "                 FEMEEM stores coord * scale in in.dat\n"
    "  overwrite      Overwrite existing directory\n\n"
    "Boundary conditions:\n"
    "  Sidesets are exported as fixed boundary edges (type 8)\n"
    "  in sin.dat.B. Each sideset triangle edge becomes a BC entry.\n\n"
    "Material assignment:\n"
    "  Cubit block IDs are mapped to FEMEEM material numbers.\n"
    "  Material constants are set to placeholder values (1.0).\n"
  );
  return help;
}

bool ExportFemeemCommand::execute(CubitCommandData &data)
{
  // Suppress Cubit Learn Edition's harmless 50k-cap ERROR.
  radia::ScopedLearnEditionFilter _lef_guard;

  std::string dirname;
  data.get_string("filename", dirname);

  double scale = 1.0;
  data.get_value("scale", scale);
  if (scale <= 0.0) {
    PRINT_WARNING("Scale must be positive. Using 1.0.\n");
    scale = 1.0;
  }

  return write_femeem(dirname, scale);
}

// ========================================================================
// Main entry: create directory and write all FEMEEM files
// ========================================================================
bool ExportFemeemCommand::write_femeem(const std::string &dirname, double scale)
{
  // Extract mesh (order 1 only)
  MeshData mesh;
  if (!mesh.extract(1))
    return false;

  // Filter: count tets and warn about non-tets
  int n_tet = 0, n_other = 0;
  for (auto &elem : mesh.elements) {
    if (elem.nv == 4 && (elem.type == TETRA4 || elem.type == TETRA))
      n_tet++;
    else
      n_other++;
  }
  if (n_tet == 0) {
    PRINT_ERROR("No tetrahedral elements found. FEMEEM requires tet meshes.\n");
    return false;
  }
  if (n_other > 0) {
    PRINT_WARNING("Skipping %d non-tet elements (hex/wedge/pyramid). "
                  "FEMEEM only supports tetrahedra.\n", n_other);
  }

  // Create output directory
  struct stat st;
  if (stat(dirname.c_str(), &st) == 0) {
    // Directory exists
  } else {
    if (_mkdir(dirname.c_str()) != 0) {
      PRINT_ERROR("Cannot create directory: %s\n", dirname.c_str());
      return false;
    }
  }

  std::string sep = "/";
  std::string path_in    = dirname + sep + "in.dat";
  std::string path_sin   = dirname + sep + "sin.dat.B";
  std::string path_sina  = dirname + sep + "sina.dat";
  std::string path_d3    = dirname + sep + "d3";

  if (!write_in_dat(path_in, mesh, scale)) return false;
  if (!write_sin_dat_B(path_sin, mesh))    return false;
  if (!write_sina_dat(path_sina))          return false;
  if (!write_d3(path_d3))                  return false;

  PRINT_INFO("Exported FEMEEM: %s (%d nodes, %d tets, scale=%.3g)\n",
             dirname.c_str(), mesh.total_node_count(), n_tet, scale);
  return true;
}

// ========================================================================
// in.dat: header + tet connectivity (2 per line) + node coordinates
// ========================================================================
bool ExportFemeemCommand::write_in_dat(const std::string &path,
                                        const MeshData &mesh, double scale)
{
  std::ofstream fid(u8_string_to_path(path));
  if (!fid.is_open()) {
    PRINT_ERROR("Cannot open file: %s\n", path.c_str());
    return false;
  }

  // Collect tet elements only, build renumbered node map
  struct TetElem {
    int n[4];  // original Cubit node IDs
  };
  std::vector<TetElem> tets;
  std::set<int> used_nodes_set;

  for (auto &elem : mesh.elements) {
    if (elem.nv != 4) continue;
    if (elem.type != TETRA4 && elem.type != TETRA) continue;
    TetElem t;
    for (int i = 0; i < 4; i++) {
      t.n[i] = elem.conn[i];
      used_nodes_set.insert(elem.conn[i]);
    }
    tets.push_back(t);
  }

  // Build old_id -> new_id map (1-based, sorted by original ID)
  std::vector<int> used_nodes(used_nodes_set.begin(), used_nodes_set.end());
  std::sort(used_nodes.begin(), used_nodes.end());
  std::map<int, int> id_map;  // old Cubit node ID -> new 1-based ID
  for (int i = 0; i < (int)used_nodes.size(); i++)
    id_map[used_nodes[i]] = i + 1;

  int npoint = (int)used_nodes.size();
  int nelem  = (int)tets.size();

  // Header: npoint nelem 0 0 scale
  char buf[128];
  snprintf(buf, sizeof(buf), "%8d%8d%8d%8d       %.3E",
           npoint, nelem, 0, 0, scale);
  fid << buf << "\n";

  // Tet connectivity: 2 tets per line (8 ints, %8d each)
  // If odd number of tets, last line has only 4 ints
  for (int i = 0; i < nelem; i += 2) {
    for (int j = 0; j < 4; j++) {
      snprintf(buf, sizeof(buf), "%8d", id_map[tets[i].n[j]]);
      fid << buf;
    }
    if (i + 1 < nelem) {
      for (int j = 0; j < 4; j++) {
        snprintf(buf, sizeof(buf), "%8d", id_map[tets[i+1].n[j]]);
        fid << buf;
      }
    }
    fid << "\n";
  }

  // Node coordinates: x y z (scaled), %24.12E each
  for (int old_id : used_nodes) {
    auto it = mesh.node_id_to_index.find(old_id);
    if (it == mesh.node_id_to_index.end()) continue;
    const auto &nd = mesh.nodes[it->second];
    snprintf(buf, sizeof(buf), "%24.12E%24.12E%24.12E",
             nd.x * scale, nd.y * scale, nd.z * scale);
    fid << buf << "\n";
  }

  fid.close();
  return true;
}

// ========================================================================
// sin.dat.B: boundary conditions + element blocks + materials + current density
// ========================================================================
bool ExportFemeemCommand::write_sin_dat_B(const std::string &path,
                                           const MeshData &mesh)
{
  std::ofstream fid(u8_string_to_path(path));
  if (!fid.is_open()) {
    PRINT_ERROR("Cannot open file: %s\n", path.c_str());
    return false;
  }

  // --- Build node renumbering (same as in_dat) ---
  std::set<int> used_nodes_set;
  std::vector<int> tet_block_ids;  // block ID per tet (in order)

  for (auto &elem : mesh.elements) {
    if (elem.nv != 4) continue;
    if (elem.type != TETRA4 && elem.type != TETRA) continue;
    for (int i = 0; i < 4; i++)
      used_nodes_set.insert(elem.conn[i]);
    tet_block_ids.push_back(elem.group_id);
  }

  std::vector<int> used_nodes(used_nodes_set.begin(), used_nodes_set.end());
  std::sort(used_nodes.begin(), used_nodes.end());
  std::map<int, int> id_map;
  for (int i = 0; i < (int)used_nodes.size(); i++)
    id_map[used_nodes[i]] = i + 1;

  int nelem = (int)tet_block_ids.size();

  // --- Boundary conditions from sidesets ---
  // Extract unique edges from sideset faces, renumber, use type 8
  struct Edge {
    int n1, n2;  // renumbered 1-based, n1 < n2
    bool operator<(const Edge &o) const {
      return (n1 != o.n1) ? (n1 < o.n1) : (n2 < o.n2);
    }
  };
  std::set<Edge> bc_edges;

  for (auto &sg : mesh.sidesets) {
    for (auto &face : sg.faces) {
      // Extract edges from face (TRI3 -> 3 edges, QUAD4 -> 4 edges)
      int nv = face.nv;
      for (int i = 0; i < nv; i++) {
        int a = face.conn[i];
        int b = face.conn[(i + 1) % nv];
        // Only include if both nodes are in tet mesh
        auto it_a = id_map.find(a);
        auto it_b = id_map.find(b);
        if (it_a == id_map.end() || it_b == id_map.end()) continue;
        int na = it_a->second, nb = it_b->second;
        if (na > nb) std::swap(na, nb);
        bc_edges.insert({na, nb});
      }
    }
  }

  int nbc = (int)bc_edges.size();

  // --- Element blocks: group consecutive tets by material ID ---
  // Map Cubit block IDs to FEMEEM material numbers (1-based)
  std::set<int> unique_block_ids(tet_block_ids.begin(), tet_block_ids.end());
  std::map<int, int> block_to_mat;
  int mat_num = 1;
  for (int bid : unique_block_ids)
    block_to_mat[bid] = mat_num++;
  int nmat = (int)unique_block_ids.size();

  // Build block ranges: start_elem, end_elem, mat_id (1-based elements)
  struct BlockRange {
    int start, end, mat;
  };
  std::vector<BlockRange> block_ranges;
  if (!tet_block_ids.empty()) {
    int cur_mat = block_to_mat[tet_block_ids[0]];
    int start = 1;
    for (int i = 1; i < nelem; i++) {
      int m = block_to_mat[tet_block_ids[i]];
      if (m != cur_mat) {
        block_ranges.push_back({start, i, cur_mat});
        start = i + 1;
        cur_mat = m;
      }
    }
    block_ranges.push_back({start, nelem, cur_mat});
  }
  int nblock = (int)block_ranges.size();

  // --- Write header ---
  char buf[256];
  snprintf(buf, sizeof(buf), "%8d%8d%8d%8d%8d",
           nbc, nblock, nmat, 0, 0);
  fid << buf << "\n";

  // --- Write boundary conditions (type 8: fixed boundary, arbitrary edges) ---
  for (auto &e : bc_edges) {
    snprintf(buf, sizeof(buf), "%8d%8d%8d%16.6f", e.n1, e.n2, 8, 0.0);
    fid << buf << "\n";
  }

  // --- Write element block assignments ---
  for (auto &br : block_ranges) {
    snprintf(buf, sizeof(buf), "%8d%8d%8d", br.start, br.end, br.mat);
    fid << buf << "\n";
  }

  // --- Write material constants (placeholder: 1.0 1.0 1.0 0.0 per material) ---
  for (int i = 0; i < nmat; i++) {
    snprintf(buf, sizeof(buf), "%16.6f%16.6f%16.6f%16.6f",
             1.0, 1.0, 1.0, 0.0);
    fid << buf << "\n";
  }

  // --- Write current density section (empty) ---
  snprintf(buf, sizeof(buf), "%8d", 0);
  fid << buf << "\n";

  fid.close();
  return true;
}

// ========================================================================
// sina.dat: analysis control parameters template (AMGTYP=1 for Cubit mesh)
// ========================================================================
bool ExportFemeemCommand::write_sina_dat(const std::string &path)
{
  std::ofstream fid(u8_string_to_path(path));
  if (!fid.is_open()) {
    PRINT_ERROR("Cannot open file: %s\n", path.c_str());
    return false;
  }

  fid << "   0   0   0   0   0   0   0  = IFVOLT,IFADAT,IFAEG,IFB0,IFERR,IFMOVE,IFWOFF\n";
  fid << "   1   0   0   0   0   1   1  = NSTEP,IFWAVE,IFCONV,IFTIME,IFCONT,IFBDAT,IFJE\n";
  fid << " 0.00000E+00 0.00000E+00   0  = DT,TIME0(sec.),IFZERO\n";
  fid << " 0.00000E+00 0.00000E+00   0  = FREQ,FREQE(Hz),IFSTEP\n";
  fid << "   0  30   0   0   0   0   0  = IFMAT,MITE,IFMAG,IFLOSS,IFFAI,IFMOV2,INITNU\n";
  fid << "   0   0                      = {MatNo.,FileNo.}, ...\n";
  fid << " 0.00000E+00 0.00000E+00      = BERR,BERR2(T)\n";
  fid << " 1.00000E+00 0.00000E+00   0  = ALFA,SIG(=1.0[++],=0.0[+-]),IFPCF\n";
  fid << " 1.00000E-07 0.00000E+00   0  = DDMIN,DDCUT(%)\n";
  fid << "   4   0   0 0.00000E+00   0  = NLFLUX,IFCAP,IFBZIP,DDAMIN,IFDDAM\n";
  fid << "   1 1.00000E+00   2 1.00000E+00   3 1.00000E-05   4 1.00000E-05              = {MatNo.,I0Flux}, ...\n";
  fid << "   0   0   0 0.00000E+00   0  = IFINED,IFGAPE,IFALFF,UUAMIN,IFEEMJ\n";
  fid << "   0   0   0 0.00000E+00   0  = NMOVE,NBTYPE,MVZONE,A0PLUS,IFBDBL\n";
  fid << "   1   0   0 0.00000E+00   0  = AMGTYP,IFNLCG,DUMMY,DITE(%),IFNODF\n";
  fid << "   0   d3t   sina.dat.t    0  = IFTMAT,FND3T[A8],FNSINA[A16],IFADDT\n";
  fid << "   0   0                      = {MatNo.,FileNo.}, ...\n";
  fid << "   0 0.00000E+00 0.00000E+00  = NLCMAT,J0ERR,J0BEKI\n";
  fid << "   0                          = {IMatNo, ... }\n";
  fid << "   0   0 Results              = IFSIGMA,OUTDIR,DIRNAME[A20]\n";
  fid << "   0 0.00000E+00              = NLTYPE,NLVAL\n";

  fid.close();
  return true;
}

// ========================================================================
// d3: file reference index
// ========================================================================
bool ExportFemeemCommand::write_d3(const std::string &path)
{
  std::ofstream fid(u8_string_to_path(path));
  if (!fid.is_open()) {
    PRINT_ERROR("Cannot open file: %s\n", path.c_str());
    return false;
  }

  fid << "07out\n";
  fid << "15in.dat\n";
  fid << "25sin.dat.B\n";
  fid << "24sina.dat\n";
  fid << "45in2d.dat\n";
  fid << "48sin2d.dat\n";

  fid.close();
  return true;
}
