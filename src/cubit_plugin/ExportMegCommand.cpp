#include "ExportMegCommand.hpp"
#include "MeshData.hpp"
#include "CubitMessage.hpp"

#include <fstream>
#include <set>
#include <cmath>
#include <algorithm>

ExportMegCommand::ExportMegCommand() {}
ExportMegCommand::~ExportMegCommand() {}

std::vector<std::string> ExportMegCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "export meg <string:label='filename',help='<filename>'> "
    "[threed] [twod] [axisymmetric] "
    "[overwrite]"
  );
  return syntax_list;
}

std::vector<std::string> ExportMegCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "export meg \"filename\" [threed|twod|axisymmetric] [overwrite]"
  );
  return help;
}

std::vector<std::string> ExportMegCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Export mesh to ELF/MAGIC MEG format.\n"
    "Block names are used as ELF element identifiers (first 4 chars + DIM).\n\n"
    "Options:\n"
    "  threed         3D analysis (default, DIM=T)\n"
    "  twod           2D planar analysis (DIM=K, z forced to 0)\n"
    "  axisymmetric   Axisymmetric analysis (DIM=R, y forced to 0)\n"
    "  overwrite      Overwrite existing file\n"
  );
  return help;
}

bool ExportMegCommand::execute(CubitCommandData &data)
{
  std::string filename;
  data.get_string("filename", filename);

  char dim = 'T';  // default: 3D
  if (data.find_keyword("twod"))
    dim = 'K';
  else if (data.find_keyword("axisymmetric"))
    dim = 'R';

  return write_meg(filename, dim);
}

// Build ELF element type string: first 4 chars of block name + DIM
std::string ExportMegCommand::elf_type(const std::string &block_name, char dim)
{
  std::string prefix = block_name.substr(0, std::min((size_t)4, block_name.size()));
  // Pad to 4 characters if shorter
  while (prefix.size() < 4) prefix += ' ';
  return prefix + dim;
}

bool ExportMegCommand::write_meg(const std::string &filename, char dim)
{
  MeshData mesh;
  if (!mesh.extract(1))  // MEG is always 1st order
    return false;

  std::ofstream fid(filename);
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
  // Build block_id -> block_name map
  // Use "block_{id}" as default name if Cubit name is not available
  // Block names come from MeshData; we use block_ids ordering
  std::unordered_map<int, std::string> block_names;
  for (int bid : mesh.block_ids) {
    block_names[bid] = "block_" + std::to_string(bid);
  }
  // Try to get actual Cubit block names via CubitInterface
  // (available when running inside Cubit)
  {
    // MeshData doesn't store block names, so we query CubitInterface
    // For simplicity, use the block_id-based name. The user should name
    // blocks appropriately in Cubit (e.g., "SOLI", "IRON", "AIR_").
  }

  fid << "* ELEMENT K\n";
  int eid = 0;

  // Volume elements (3D only, DIM='T')
  // 2D/Axi: tri/quad only, but we output whatever is in blocks
  for (auto &elem : mesh.elements) {
    eid++;
    std::string etype = elf_type(block_names[elem.group_id], dim);
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
      ? "side_" + std::to_string(sg.id) : sg.name;
    for (auto &face : sg.faces) {
      eid++;
      std::string etype = elf_type(ss_name, dim);
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
  // TODO: MeshData::write_companion_json(filename);
  return true;
}
