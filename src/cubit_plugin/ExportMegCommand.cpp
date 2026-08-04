#include "ExportMegCommand.hpp"
#include "MeshData.hpp"
#include "RadiaMessageFilter.hpp"
#include "CubitMessage.hpp"
#include "CubitInterface.hpp"
#include "utf8_path.hpp"

#include <fstream>
#include <set>
#include <cmath>
#include <algorithm>
#include <unordered_map>

ExportMegCommand::ExportMegCommand() {}
ExportMegCommand::~ExportMegCommand() {}

std::vector<std::string> ExportMegCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "export meg <string:label='filename',help='<filename>'> "
    "[threed] [twod] [axisymmetric] "
    "[labels <string:label='labels',help='<block:prefix,...>'>] "
    "[overwrite]"
  );
  return syntax_list;
}

std::vector<std::string> ExportMegCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "export meg \"filename\" [threed|twod|axisymmetric] [labels \"1:MMB,2:MWL,...\"] [overwrite]"
  );
  return help;
}

std::vector<std::string> ExportMegCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Export mesh to ELF/MAGIC MEG format.\n"
    "Block names define ELF element types (3-char physics prefix).\n"
    "Label = prefix(3) + nodecount(1) + DIM(1), e.g. MMB8T.\n\n"
    "Block naming convention (case-insensitive):\n"
    "  MMB  Magnetic body (iron core, ferromagnetic)\n"
    "  MMS  Magnetic shell\n"
    "  MMT  Magnetic thin (line/surface elements)\n"
    "  MWL  Nonlinear magnet with fixed local axis\n"
    "  MWV  Nonlinear magnet with direction vectors\n"
    "  MCO  Current conductor\n"
    "  AIR  Air region (non-magnetic, uses MMB)\n\n"
    "Node count is auto-detected: HEX8->8, TET4->4, WEDGE6->6,\n"
    "  QUAD4->4, TRI3->3, PYRAMID5->8 (degenerate hex).\n\n"
    "Options:\n"
    "  threed         3D analysis (default, DIM=T)\n"
    "  twod           2D planar analysis (DIM=K, z forced to 0)\n"
    "  axisymmetric   Axisymmetric analysis (DIM=R, y forced to 0)\n"
    "  overwrite      Overwrite existing file\n"
  );
  return help;
}

// Parse "1:MMB,2:MWL,3:MCO" -> map<int,string>
static std::unordered_map<int, std::string> parse_labels(const std::string &s)
{
  std::unordered_map<int, std::string> result;
  if (s.empty()) return result;

  size_t pos = 0;
  while (pos < s.size()) {
    size_t colon = s.find(':', pos);
    if (colon == std::string::npos) break;
    size_t comma = s.find(',', colon);
    if (comma == std::string::npos) comma = s.size();

    int bid = std::atoi(s.substr(pos, colon - pos).c_str());
    std::string prefix = s.substr(colon + 1, comma - colon - 1);
    if (bid > 0 && prefix.size() == 3)
      result[bid] = prefix;

    pos = comma + 1;
  }
  return result;
}

bool ExportMegCommand::execute(CubitCommandData &data)
{
  // Suppress Cubit Learn Edition's harmless 50k-cap ERROR.
  radia::ScopedLearnEditionFilter _lef_guard;

  std::string filename;
  data.get_string("filename", filename);

  char dim = 'T';  // default: 3D
  if (data.find_keyword("twod"))
    dim = 'K';
  else if (data.find_keyword("axisymmetric"))
    dim = 'R';

  // Parse per-block labels
  std::string labels_str;
  data.get_string("labels", labels_str);
  auto label_map = parse_labels(labels_str);

  return write_meg(filename, dim, label_map);
}

// Valid ELF 3-char physics prefixes
static const char* VALID_ELF_PREFIXES[] = {
  "MMB", "MMS", "MMT", "MMP", "MWL", "MWV", "MCO", nullptr
};

// Extract 3-char ELF prefix from block name (case-insensitive match)
// Returns "MMB" as default if no valid prefix found.
static std::string extract_elf_prefix(const std::string &block_name)
{
  std::string upper = block_name;
  for (auto &c : upper) c = (char)toupper(c);

  // Check if name starts with a valid ELF prefix
  for (int i = 0; VALID_ELF_PREFIXES[i]; i++) {
    if (upper.size() >= 3 && upper.substr(0, 3) == VALID_ELF_PREFIXES[i])
      return VALID_ELF_PREFIXES[i];
  }

  // "AIR" or unknown -> default to MMB (air is still a magnetic region in ELF)
  return "MMB";
}

// Determine effective node count for ELF label
// Pyramids are exported as degenerate 8-node hex in ELF
static int elf_node_count(int nv, ElementType type)
{
  bool is_pyramid = (nv == 5 && (type == PYRAMID5 || type == PYRAMID));
  return is_pyramid ? 8 : nv;
}

// Build ELF element type string: prefix(3) + nodecount(1) + DIM(1)
// e.g. "MMB8T" for 8-node hex magnetic body in 3D
std::string ExportMegCommand::elf_type(const std::string &prefix, int nv,
                                        ElementType type, char dim)
{
  int nc = elf_node_count(nv, type);
  return prefix + std::to_string(nc) + dim;
}

bool ExportMegCommand::write_meg(const std::string &filename, char dim,
                                  const std::unordered_map<int, std::string> &label_map)
{
  MeshData mesh;
  if (!mesh.extract(1))  // MEG is always 1st order
    return false;

  std::ofstream fid(u8_string_to_path(filename));
  if (!fid.is_open()) {
    PRINT_ERROR("Cannot open file: %s\n", filename.c_str());
    return false;
  }

  // --- ELF/MESH header ---
  fid << "BOOK  MEP  3.50\n";
  fid << "* ELF/MESH VERSION 7.3.0\n";
  fid << "* SOLVER = ELF/MAGIC\n";
  fid << "MGSC 0.001\n";

  // --- Nodes (MGR1) ---
  fid << "* NODE\n";
  int num_nodes = mesh.num_original_nodes;
  for (int i = 0; i < num_nodes; i++) {
    auto &nd = mesh.nodes[i];
    double x = nd.x;
    double y = nd.y;
    double z = nd.z;

    if (dim == 'K') {
      z = 0.0;  // 2D: force z=0
    } else if (dim == 'R') {
      y = 0.0;  // Axisymmetric: R-Z plane (y=0)
    }

    fid << "MGR1 " << nd.id << " 0 " << x << " " << y << " " << z << "\n";
  }

  // --- Elements ---
  // Build block_id -> ELF prefix map
  // Priority: label_map (from GUI) > block name (from Cubit) > default "MMB"
  // Use get_block_name (NOT get_entity_name("block", ...)) — the
  // generic accessor calls get_mref_entity which raises
  // "Invalid list type for the get_mref_entity function: block".
  std::unordered_map<int, std::string> block_prefix;
  for (int bid : mesh.block_ids) {
    auto it = label_map.find(bid);
    if (it != label_map.end()) {
      block_prefix[bid] = it->second;
    } else {
      std::string name = CubitInterface::get_block_name(bid);
      block_prefix[bid] = extract_elf_prefix(name.empty() ? "MMB" : name);
    }
  }

  // Print block->label mapping for user verification
  PRINT_INFO("MEG block -> ELF label mapping:\n");
  for (int bid : mesh.block_ids) {
    std::string name = CubitInterface::get_block_name(bid);
    PRINT_INFO("  Block %d \"%s\" -> %s*%c\n",
               bid, name.c_str(), block_prefix[bid].c_str(), dim);
  }

  fid << "* ELEMENT K\n";
  int eid = 0;

  for (auto &elem : mesh.elements) {
    eid++;
    std::string prefix = block_prefix.count(elem.group_id)
                         ? block_prefix[elem.group_id] : "MMB";
    int nc = elf_node_count(elem.nv, elem.type);
    std::string etype = prefix + std::to_string(nc) + dim;
    fid << etype << " " << eid << " 0 " << elem.group_id;

    int nv = elem.nv;
    ElementType et = elem.type;

    // Pyramid special case: ELF uses 8-node degenerate hex for pyramids
    bool is_pyramid = (nv == 5 && (et == PYRAMID5 || et == PYRAMID));
    if (is_pyramid && dim == 'T') {
      // Repeat apex node 4 times: n0 n1 n2 n3 n4 n4 n4 n4
      for (int j = 0; j < 4; j++) fid << " " << elem.conn[j];
      for (int j = 0; j < 4; j++) fid << " " << elem.conn[4];
    } else {
      for (int nid : elem.conn) fid << " " << nid;
    }
    fid << "\n";
  }

  // Sideset face elements
  for (auto &sg : mesh.sidesets) {
    std::string ss_name = sg.name.empty()
      ? "MMT" : sg.name;  // default: magnetic thin for surface elements
    for (auto &face : sg.faces) {
      eid++;
      std::string ss_prefix = extract_elf_prefix(ss_name);
      std::string etype = elf_type(ss_prefix, face.nv, face.type, dim);
      fid << etype << " " << eid << " 0 " << sg.id;
      for (int nid : face.conn) fid << " " << nid;
      fid << "\n";
    }
  }

  // --- MGR2 (spatial nodes from "SPACE" nodesets/blocks) ---
  fid << "* NODE\n";
  {
    int mgr2_id = 0;
    // Check nodesets named "SPACE" (case-insensitive)
    for (auto &ng : mesh.nodesets) {
      std::string upper_name = ng.name;
      for (auto &c : upper_name) c = (char)toupper(c);
      if (upper_name.find("SPACE") == std::string::npos) continue;
      for (int nid : ng.node_ids) {
        auto it = mesh.node_id_to_index.find(nid);
        if (it == mesh.node_id_to_index.end()) continue;
        auto &nd = mesh.nodes[it->second];
        mgr2_id++;
        double x = nd.x, y = nd.y, z = nd.z;
        if (dim == 'K') z = 0.0;
        else if (dim == 'R') y = 0.0;
        fid << "MGR2 " << mgr2_id << " 0 " << x << " " << y << " " << z << "\n";
      }
    }
    // Also check sidesets named "SPACE"
    for (auto &sg : mesh.sidesets) {
      std::string upper_name = sg.name;
      for (auto &c : upper_name) c = (char)toupper(c);
      if (upper_name.find("SPACE") == std::string::npos) continue;
      // Collect unique nodes from faces
      std::set<int> space_nodes;
      for (auto &face : sg.faces)
        for (int nid : face.conn) space_nodes.insert(nid);
      for (int nid : space_nodes) {
        auto it = mesh.node_id_to_index.find(nid);
        if (it == mesh.node_id_to_index.end()) continue;
        auto &nd = mesh.nodes[it->second];
        mgr2_id++;
        double x = nd.x, y = nd.y, z = nd.z;
        if (dim == 'K') z = 0.0;
        else if (dim == 'R') y = 0.0;
        fid << "MGR2 " << mgr2_id << " 0 " << x << " " << y << " " << z << "\n";
      }
    }
    if (mgr2_id > 0)
      PRINT_INFO("  MGR2: %d spatial nodes from SPACE nodesets/sidesets\n", mgr2_id);
  }

  // --- Footer ---
  fid << "BOOK  END\n";
  fid.close();

  const char *dim_name = (dim == 'T') ? "3D" : (dim == 'K') ? "2D" : "Axisymmetric";
  PRINT_INFO("Exported MEG (%s): %s (%d nodes, %d elements)\n",
             dim_name, filename.c_str(), num_nodes, eid);
  return true;
}
