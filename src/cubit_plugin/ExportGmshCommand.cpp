#include "ExportGmshCommand.hpp"
#include "MeshData.hpp"
#include "NetgenCurver.hpp"
#include "CubitMessage.hpp"

#include <fstream>
#include <map>
#include <unordered_map>

ExportGmshCommand::ExportGmshCommand() {}
ExportGmshCommand::~ExportGmshCommand() {}

std::vector<std::string> ExportGmshCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "export gmsh <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1-4>'>] "
    "[version <value:label='version',help='<2 or 4>'>] "
    "[dimension <value:label='dimension',help='<2 or 3>'>] "
    "[overwrite]"
  );
  return syntax_list;
}

std::vector<std::string> ExportGmshCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "export gmsh \"filename\" [order {1|2|3|4}] [version {2|4}] [dimension {2|3}] [overwrite]"
  );
  return help;
}

std::vector<std::string> ExportGmshCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Export mesh to Gmsh MSH format.\n"
    "Block assignment is NOT required.\n"
    "Sidesets exported as surface elements. Nodesets as comments.\n\n"
    "Options:\n"
    "  order 1      1st-order elements (default)\n"
    "  order 2      2nd-order (edge mid-nodes)\n"
    "  order 3-4    Higher order (requires Netgen)\n"
    "  version 2    Gmsh v2.2 format (default)\n"
    "  version 4    Gmsh v4.1 format\n"
    "  dimension 3  3D mode (default)\n"
    "  dimension 2  2D mode\n"
    "  overwrite    Overwrite existing file\n"
  );
  return help;
}

bool ExportGmshCommand::execute(CubitCommandData &data)
{
  std::string filename;
  data.get_string("filename", filename);

  int order = 1;
  data.get_value("order", order);
  if (order < 1) order = 1;
  if (order > 4) {
    PRINT_WARNING("order %d not supported. Using order 4.\n", order);
    order = 4;
  }
#ifndef HAVE_NETGEN
  if (order > 2) {
    PRINT_WARNING("order %d requires Netgen. Falling back to order 2.\n", order);
    order = 2;
  }
#endif

  int ver = 2;
  data.get_value("version", ver);
  std::string version = (ver >= 4) ? "4.1" : "2.2";

  return write_gmsh(filename, version, order);
}

bool ExportGmshCommand::write_gmsh(const std::string &filename,
                                    const std::string &version,
                                    int order)
{
  MeshData mesh;
  if (!mesh.extract(order))
    return false;

  if (version == "4.1")
    return write_gmsh_v41(filename, mesh);

  return write_gmsh_v22(filename, mesh);
}

// ========================================================================
// Gmsh element type mapping
// ========================================================================
int ExportGmshCommand::gmsh_type(const MeshElement &elem, int order)
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

  if (order == 1) {
    if (is_tet)  return 4;
    if (is_hex)  return 5;
    if (is_wed)  return 6;
    if (is_pyr)  return 7;
    if (is_tri)  return 2;
    if (is_quad) return 3;
    if (is_line) return 1;
    if (is_pt)   return 15;
  }
  else if (order == 2) {
    if (is_tet)  return 11;  // TET10
    if (is_hex)  return 17;  // HEX20
    if (is_wed)  return 18;  // WEDGE15
    if (is_pyr)  return 19;  // PYRAMID13
    if (is_tri)  return 9;   // TRI6
    if (is_quad) return 16;  // QUAD8
    if (is_line) return 8;   // LINE3
  }
  else if (order == 3) {
    if (is_tet) return 29;   // TET20
    if (is_tri) return 21;   // TRI10
    // HEX, WEDGE, PYRAMID order 3 not yet mapped
  }
  else if (order == 4) {
    if (is_tet) return 30;   // TET35
    if (is_tri) return 23;   // TRI15
  }

  // Fallback: linear type
  if (is_tet)  return 4;
  if (is_hex)  return 5;
  if (is_wed)  return 6;
  if (is_pyr)  return 7;
  if (is_tri)  return 2;
  if (is_quad) return 3;
  if (is_line) return 1;
  return 0;
}

// ========================================================================
// Gmsh v2.2 writer — unified for all orders, blocks + sidesets
// ========================================================================
bool ExportGmshCommand::write_gmsh_v22(const std::string &filename,
                                        const MeshData &mesh)
{
  std::ofstream fid(filename);
  if (!fid.is_open()) {
    PRINT_ERROR("Cannot open file: %s\n", filename.c_str());
    return false;
  }

  fid << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

  // --- PhysicalNames: blocks (auto dim) + sidesets (BND) + nodesets (BBBND) ---
  //
  // Detect dimension per block from actual element types:
  //   tet/hex/wedge/pyramid -> 3, tri/quad -> 2, line -> 1, point -> 0
  std::unordered_map<int, int> block_dim;
  for (auto &elem : mesh.elements) {
    int dim = 3;  // default
    int nv = elem.nv;
    if (nv == 3 && (elem.type == TRI3 || elem.type == CUBIT_TRI
                 || elem.type == TRISHELL || elem.type == TRISHELL3))
      dim = 2;
    else if (nv == 4 && (elem.type == QUAD4 || elem.type == QUAD
                      || elem.type == SHEL || elem.type == SHELL4))
      dim = 2;
    else if (elem.type == BAR || elem.type == BAR2 || elem.type == BAR3)
      dim = 1;
    else if (nv == 1)
      dim = 0;

    auto it = block_dim.find(elem.group_id);
    if (it == block_dim.end() || dim > it->second)
      block_dim[elem.group_id] = dim;
  }

  {
    int nphys = (int)mesh.block_ids.size() + (int)mesh.sidesets.size()
                + (int)mesh.nodesets.size();
    if (nphys > 0) {
      fid << "$PhysicalNames\n" << nphys << "\n";
      // Blocks: auto-detected dimension
      for (int bid : mesh.block_ids) {
        int dim = 3;
        auto it = block_dim.find(bid);
        if (it != block_dim.end()) dim = it->second;
        fid << dim << " " << bid << " \"block_" << bid << "\"\n";
      }
      // Sidesets: dimension 2 (BND)
      for (auto &sg : mesh.sidesets) {
        fid << "2 " << sg.id << " \"";
        fid << (sg.name.empty() ? "sideset_" + std::to_string(sg.id) : sg.name);
        fid << "\"\n";
      }
      // Nodesets: dimension 0 (BBBND)
      for (auto &ng : mesh.nodesets) {
        fid << "0 " << ng.id << " \"";
        fid << (ng.name.empty() ? "nodeset_" + std::to_string(ng.id) : ng.name);
        fid << "\"\n";
      }
      fid << "$EndPhysicalNames\n";
    }
  }

  // --- Nodes ---
  fid << "$Nodes\n" << mesh.total_node_count() << "\n";
  for (auto &nd : mesh.nodes)
    fid << nd.id << " " << nd.x << " " << nd.y << " " << nd.z << "\n";
  fid << "$EndNodes\n";

  // --- Elements: block elements + sideset faces + nodeset points ---
  // Count: block elems + sideset faces + nodeset nodes (as POINT elements)
  int nodeset_point_count = 0;
  for (auto &ng : mesh.nodesets) nodeset_point_count += (int)ng.node_ids.size();
  int total = mesh.total_element_count() + nodeset_point_count;
  fid << "$Elements\n" << total << "\n";

  int eid = 0;
  int order = mesh.order;

  // Block elements
  int skipped = 0;
  for (auto &elem : mesh.elements) {
    int gtype = gmsh_type(elem, order);
    if (gtype == 0) {
      skipped++;
      PRINT_WARNING("Skipped element: type=%d nv=%d block=%d\n",
                    (int)elem.type, elem.nv, elem.group_id);
      continue;
    }
    eid++;

    const auto &c = (order >= 2 && !elem.ho_conn.empty()) ? elem.ho_conn : elem.conn;
    fid << eid << " " << gtype << " 2 " << elem.group_id << " " << elem.group_id;
    for (int nid : c) fid << " " << nid;
    fid << "\n";
  }

  // Sideset face elements
  for (auto &sg : mesh.sidesets) {
    for (auto &face : sg.faces) {
      eid++;
      int gtype = gmsh_type(face, order);
      if (gtype == 0) continue;

      const auto &c = (order >= 2 && !face.ho_conn.empty()) ? face.ho_conn : face.conn;
      fid << eid << " " << gtype << " 2 " << sg.id << " " << sg.id;
      for (int nid : c) fid << " " << nid;
      fid << "\n";
    }
  }

  // Nodeset nodes as POINT elements (BBBND for NGSolve)
  for (auto &ng : mesh.nodesets) {
    for (int nid : ng.node_ids) {
      eid++;
      fid << eid << " 15 2 " << ng.id << " " << ng.id << " " << nid << "\n";
    }
  }

  fid << "$EndElements\n";

  if (skipped > 0)
    PRINT_WARNING("Skipped %d elements with unknown type.\n", skipped);

  // Rewrite element count at the correct position
  // (we used 'total' but some may have been skipped)
  fid.close();

  PRINT_INFO("Exported Gmsh v2.2 (order %d%s): %s (%d nodes, %d elements)\n",
             mesh.order, mesh.has_netgen ? ", NetgenCurver" : "",
             filename.c_str(), mesh.total_node_count(), eid);
  // TODO: MeshData::write_companion_json(filename);
  return true;
}

// ========================================================================
// Element topological dimension
// ========================================================================
int ExportGmshCommand::elem_dim(const MeshElement &elem)
{
  int nv = elem.nv;
  if (nv == 4 && (elem.type == TETRA4 || elem.type == TETRA))   return 3;
  if (nv == 8 && (elem.type == HEX8   || elem.type == HEX))     return 3;
  if (nv == 6 && (elem.type == WEDGE6 || elem.type == WEDGE))   return 3;
  if (nv == 5 && (elem.type == PYRAMID5 || elem.type == PYRAMID)) return 3;
  if (nv == 3 && (elem.type == TRI3   || elem.type == CUBIT_TRI
               || elem.type == TRISHELL || elem.type == TRISHELL3)) return 2;
  if (nv == 4 && (elem.type == QUAD4  || elem.type == QUAD
               || elem.type == SHEL  || elem.type == SHELL4))   return 2;
  if (elem.type == BAR || elem.type == BAR2 || elem.type == BAR3) return 1;
  if (nv == 1) return 0;
  return 3; // default
}

// ========================================================================
// Gmsh v4.1 writer
//
// Format reference: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
// Sections: $MeshFormat, $PhysicalNames, $Entities, $Nodes, $Elements
// ========================================================================
bool ExportGmshCommand::write_gmsh_v41(const std::string &filename,
                                        const MeshData &mesh)
{
  std::ofstream fid(filename);
  if (!fid.is_open()) {
    PRINT_ERROR("Cannot open file: %s\n", filename.c_str());
    return false;
  }

  int order = mesh.order;

  // --- Detect block dimensions and build bounding boxes ---
  struct BlockInfo {
    int dim;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    BlockInfo() : dim(0),
      xmin(1e30), ymin(1e30), zmin(1e30),
      xmax(-1e30), ymax(-1e30), zmax(-1e30) {}
    void expand(double x, double y, double z) {
      if (x < xmin) xmin = x; if (x > xmax) xmax = x;
      if (y < ymin) ymin = y; if (y > ymax) ymax = y;
      if (z < zmin) zmin = z; if (z > zmax) zmax = z;
    }
  };

  // bid -> BlockInfo
  std::map<int, BlockInfo> block_info;
  for (int bid : mesh.block_ids)
    block_info[bid] = BlockInfo();

  for (auto &el : mesh.elements) {
    auto &bi = block_info[el.group_id];
    int d = elem_dim(el);
    if (d > bi.dim) bi.dim = d;
    for (int nid : el.conn) {
      auto it = mesh.node_id_to_index.find(nid);
      if (it != mesh.node_id_to_index.end()) {
        auto &nd = mesh.nodes[it->second];
        bi.expand(nd.x, nd.y, nd.z);
      }
    }
  }

  // Sideset info
  std::map<int, BlockInfo> sideset_info;
  for (auto &sg : mesh.sidesets) {
    auto &si = sideset_info[sg.id];
    si.dim = 2;
    for (auto &face : sg.faces) {
      for (int nid : face.conn) {
        auto it = mesh.node_id_to_index.find(nid);
        if (it != mesh.node_id_to_index.end()) {
          auto &nd = mesh.nodes[it->second];
          si.expand(nd.x, nd.y, nd.z);
        }
      }
    }
  }

  // Nodeset info
  std::map<int, BlockInfo> nodeset_info;
  for (auto &ng : mesh.nodesets) {
    auto &ni = nodeset_info[ng.id];
    ni.dim = 0;
    for (int nid : ng.node_ids) {
      auto it = mesh.node_id_to_index.find(nid);
      if (it != mesh.node_id_to_index.end()) {
        auto &nd = mesh.nodes[it->second];
        ni.expand(nd.x, nd.y, nd.z);
      }
    }
  }

  // ---- $MeshFormat ----
  fid << "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n";

  // ---- $PhysicalNames ----
  int nphys = (int)block_info.size() + (int)sideset_info.size()
              + (int)nodeset_info.size();
  if (nphys > 0) {
    fid << "$PhysicalNames\n" << nphys << "\n";
    for (auto &[bid, bi] : block_info)
      fid << bi.dim << " " << bid << " \"block_" << bid << "\"\n";
    for (auto &sg : mesh.sidesets) {
      fid << "2 " << sg.id << " \"";
      fid << (sg.name.empty() ? "sideset_" + std::to_string(sg.id) : sg.name);
      fid << "\"\n";
    }
    for (auto &ng : mesh.nodesets) {
      fid << "0 " << ng.id << " \"";
      fid << (ng.name.empty() ? "nodeset_" + std::to_string(ng.id) : ng.name);
      fid << "\"\n";
    }
    fid << "$EndPhysicalNames\n";
  }

  // ---- $Entities ----
  // Count entities per dimension
  int n_ent[4] = {0, 0, 0, 0}; // dim 0..3
  // Each block/sideset/nodeset is one entity
  // Assign unique entity tags: blocks first, then sidesets, then nodesets
  struct EntityEntry {
    int dim, tag, phys_tag;
    double xmin, ymin, zmin, xmax, ymax, zmax;
  };
  std::vector<EntityEntry> entities;

  for (auto &[bid, bi] : block_info) {
    n_ent[bi.dim]++;
    entities.push_back({bi.dim, bid, bid,
                        bi.xmin, bi.ymin, bi.zmin,
                        bi.xmax, bi.ymax, bi.zmax});
  }
  for (auto &[sid, si] : sideset_info) {
    n_ent[si.dim]++;
    entities.push_back({si.dim, sid, sid,
                        si.xmin, si.ymin, si.zmin,
                        si.xmax, si.ymax, si.zmax});
  }
  for (auto &[nid, ni] : nodeset_info) {
    n_ent[ni.dim]++;
    entities.push_back({ni.dim, nid, nid,
                        ni.xmin, ni.ymin, ni.zmin,
                        ni.xmax, ni.ymax, ni.zmax});
  }

  fid << "$Entities\n";
  fid << n_ent[0] << " " << n_ent[1] << " " << n_ent[2] << " " << n_ent[3] << "\n";

  // Write entities grouped by dimension (0, 1, 2, 3)
  for (int d = 0; d <= 3; d++) {
    for (auto &e : entities) {
      if (e.dim != d) continue;
      if (d == 0) {
        // point entity: tag x y z numPhysicalTags [phys...]
        fid << e.tag << " "
            << e.xmin << " " << e.ymin << " " << e.zmin << " "
            << "1 " << e.phys_tag << "\n";
      } else {
        // curve/surface/volume entity:
        // tag xMin yMin zMin xMax yMax zMax numPhysicalTags [phys...] numBounding [bound...]
        fid << e.tag << " "
            << e.xmin << " " << e.ymin << " " << e.zmin << " "
            << e.xmax << " " << e.ymax << " " << e.zmax << " "
            << "1 " << e.phys_tag << " 0\n";
      }
    }
  }
  fid << "$EndEntities\n";

  // ---- $Nodes ----
  // v4.1: numEntityBlocks numNodes minNodeTag maxNodeTag
  // We put all nodes in one entity block per block/entity.
  // Simpler: one big entity block with all nodes.
  int n_total_nodes = mesh.total_node_count();
  int min_nid = mesh.nodes.empty() ? 0 : mesh.nodes.front().id;
  int max_nid = mesh.nodes.empty() ? 0 : mesh.nodes.back().id;
  for (auto &nd : mesh.nodes) {
    if (nd.id < min_nid) min_nid = nd.id;
    if (nd.id > max_nid) max_nid = nd.id;
  }

  fid << "$Nodes\n";
  // One entity block containing all nodes (dim=3, entityTag=1, parametric=0)
  // Pick the first volume entity tag, or 1 as fallback
  int node_entity_tag = 1;
  if (!block_info.empty())
    node_entity_tag = block_info.begin()->first;
  int node_entity_dim = 3;
  if (!block_info.empty())
    node_entity_dim = block_info.begin()->second.dim;

  fid << "1 " << n_total_nodes << " " << min_nid << " " << max_nid << "\n";
  fid << node_entity_dim << " " << node_entity_tag << " 0 " << n_total_nodes << "\n";
  // Node tags
  for (auto &nd : mesh.nodes)
    fid << nd.id << "\n";
  // Node coordinates
  for (auto &nd : mesh.nodes)
    fid << nd.x << " " << nd.y << " " << nd.z << "\n";
  fid << "$EndNodes\n";

  // ---- $Elements ----
  // Build entity blocks: group elements by (entity_tag, gmsh_type)
  struct ElemBlock {
    int dim;
    int entity_tag;
    int gmsh_type;
    std::vector<std::vector<int>> conns; // each element's connectivity
  };

  // Key: (entity_tag, gmsh_type)
  using BlockKey = std::pair<int, int>;
  std::map<BlockKey, ElemBlock> elem_blocks;

  int total_elems = 0;
  int skipped = 0;

  // Block elements
  for (auto &el : mesh.elements) {
    int gtype = gmsh_type(el, order);
    if (gtype == 0) { skipped++; continue; }
    BlockKey key(el.group_id, gtype);
    auto &blk = elem_blocks[key];
    blk.dim = elem_dim(el);
    blk.entity_tag = el.group_id;
    blk.gmsh_type = gtype;
    const auto &c = (order >= 2 && !el.ho_conn.empty()) ? el.ho_conn : el.conn;
    blk.conns.push_back(c);
    total_elems++;
  }

  // Sideset face elements
  for (auto &sg : mesh.sidesets) {
    for (auto &face : sg.faces) {
      int gtype = gmsh_type(face, order);
      if (gtype == 0) continue;
      BlockKey key(sg.id, gtype);
      auto &blk = elem_blocks[key];
      blk.dim = elem_dim(face);
      blk.entity_tag = sg.id;
      blk.gmsh_type = gtype;
      const auto &c = (order >= 2 && !face.ho_conn.empty()) ? face.ho_conn : face.conn;
      blk.conns.push_back(c);
      total_elems++;
    }
  }

  // Nodeset nodes as POINT elements
  for (auto &ng : mesh.nodesets) {
    if (ng.node_ids.empty()) continue;
    BlockKey key(ng.id, 15);
    auto &blk = elem_blocks[key];
    blk.dim = 0;
    blk.entity_tag = ng.id;
    blk.gmsh_type = 15;
    for (int nid : ng.node_ids) {
      blk.conns.push_back({nid});
      total_elems++;
    }
  }

  int min_eid = 1;
  int max_eid = total_elems;

  fid << "$Elements\n";
  fid << elem_blocks.size() << " " << total_elems
      << " " << min_eid << " " << max_eid << "\n";

  int eid = 0;
  for (auto &[key, blk] : elem_blocks) {
    fid << blk.dim << " " << blk.entity_tag << " " << blk.gmsh_type
        << " " << blk.conns.size() << "\n";
    for (auto &c : blk.conns) {
      eid++;
      fid << eid;
      for (int nid : c) fid << " " << nid;
      fid << "\n";
    }
  }

  fid << "$EndElements\n";
  fid.close();

  if (skipped > 0)
    PRINT_WARNING("Skipped %d elements with unknown type.\n", skipped);

  PRINT_INFO("Exported Gmsh v4.1 (order %d%s): %s (%d nodes, %d elements, %d blocks)\n",
             mesh.order, mesh.has_netgen ? ", NetgenCurver" : "",
             filename.c_str(), n_total_nodes, total_elems, (int)elem_blocks.size());
  // TODO: MeshData::write_companion_json(filename);
  return true;
}
