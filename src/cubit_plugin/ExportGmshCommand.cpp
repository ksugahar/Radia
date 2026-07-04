#include "ExportGmshCommand.hpp"
#include "MeshData.hpp"
#include "NetgenCurver.hpp"
#include "RadiaMessageFilter.hpp"
#include "CubitMessage.hpp"
#include "CubitInterface.hpp"
#include "utf8_path.hpp"

#include <fstream>
#include <map>
#include <unordered_map>

namespace {

std::string replace_extension(const std::string &path, const std::string &ext)
{
  auto slash = path.find_last_of("/\\");
  auto dot = path.find_last_of('.');
  if (dot == std::string::npos || (slash != std::string::npos && dot < slash))
    return path + ext;
  return path.substr(0, dot) + ext;
}

std::string basename_of(const std::string &path)
{
  auto sep = path.find_last_of("/\\");
  if (sep == std::string::npos)
    return path;
  return path.substr(sep + 1);
}

void write_gmsh_display_options(std::ofstream &out)
{
  out << "General.Orthographic = 1;\n";
  out << "General.Trackball = 0;\n";
  out << "General.RotationX = -68;\n";
  out << "General.RotationY = 0;\n";
  out << "General.RotationZ = 0;\n";
  out << "General.RotationCenterGravity = 1;\n";
  out << "Mesh.NumSubEdges = 4;\n";
  out << "Mesh.SurfaceFaces = 1;\n";
  out << "Mesh.SurfaceEdges = 1;\n";
  out << "Mesh.VolumeEdges = 0;\n";
  out << "Mesh.VolumeFaces = 0;\n";
  out << "Mesh.ColorCarousel = 2;\n";
}

bool write_gmsh_launch_companions(const std::string &msh_filename)
{
  const std::string msh_base = basename_of(msh_filename);
  const std::string geo = replace_extension(msh_filename, ".geo");
  const std::string geo_opt = geo + ".opt";
  const std::string msh_opt = msh_filename + ".opt";

  std::ofstream gf(u8_string_to_path(geo));
  if (!gf.is_open()) {
    PRINT_WARNING("GMSH companion: cannot write %s\n", geo.c_str());
    return false;
  }
  gf << "// Auto-generated Radia GMSH launch companion for " << msh_base << "\n";
  gf << "// Open this .geo for normal review; open .msh only for raw data inspection.\n";
  gf << "Merge \"" << msh_base << "\";\n\n";
  write_gmsh_display_options(gf);
  gf.close();

  std::ofstream go(u8_string_to_path(geo_opt));
  if (go.is_open()) {
    go << "// Auto-generated display options for " << basename_of(geo) << "\n";
    write_gmsh_display_options(go);
    go.close();
  } else {
    PRINT_WARNING("GMSH companion: cannot write %s\n", geo_opt.c_str());
  }

  std::ofstream mo(u8_string_to_path(msh_opt));
  if (mo.is_open()) {
    mo << "// Auto-generated raw mesh/data inspection options for " << msh_base << "\n";
    write_gmsh_display_options(mo);
    mo.close();
  } else {
    PRINT_WARNING("GMSH companion: cannot write %s\n", msh_opt.c_str());
  }

  PRINT_INFO("Companions: %s, %s, %s\n",
             geo.c_str(), geo_opt.c_str(), msh_opt.c_str());
  return true;
}

} // namespace

ExportGmshCommand::ExportGmshCommand() {}
ExportGmshCommand::~ExportGmshCommand() {}

std::vector<std::string> ExportGmshCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "export gmsh <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1-5>'>] "
    "[dimension <value:label='dimension',help='<2 or 3>'>] "
    "[overwrite]"
  );
  return syntax_list;
}

std::vector<std::string> ExportGmshCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "export gmsh \"filename\" [order {1|2|3|4|5}] [dimension {2|3}] [overwrite]"
  );
  return help;
}

std::vector<std::string> ExportGmshCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Export mesh to Gmsh MSH v4.1 format (lab-wide standard).\n"
    "Block assignment is NOT required.\n"
    "Sidesets exported as surface elements. Nodesets as comments.\n\n"
    "Writes filename.geo, filename.geo.opt, and filename.msh.opt.\n"
    "Open the .geo file for normal Radia review; the .msh is raw data.\n\n"
    "Options:\n"
    "  order 1      1st-order elements (default)\n"
    "  order 2      2nd-order (edge mid-nodes)\n"
    "  order 3-5    Higher order (requires Netgen, wedge limited to order 2)\n"
    "  dimension 3  3D mode (default)\n"
    "  dimension 2  2D mode\n"
    "  overwrite    Overwrite existing file\n"
  );
  return help;
}

bool ExportGmshCommand::execute(CubitCommandData &data)
{
  // Suppress Cubit Learn Edition's harmless 50k-cap ERROR.  Radia
  // completes the export regardless; the ERROR is misleading noise.
  radia::ScopedLearnEditionFilter _lef_guard;

  std::string filename;
  data.get_string("filename", filename);

  int order = 1;
  data.get_value("order", order);
  if (order < 1) order = 1;
  // GMSH export supports order 1-3. Order 4-5 face/volume interior nodes
  // are not reliably extracted by NetgenCurver (missing faces filled with
  // linear interpolation -> negative Jacobians in GMSH). Use .vol for p>=4.
  if (order > 3) {
    PRINT_ERROR("GMSH export: order %d not supported. "
                "Use 'export netgen' for order 4-5.\n", order);
    return false;
  }
#ifndef HAVE_NETGEN
  if (order > 1) {
    PRINT_WARNING("order %d requires Netgen. Falling back to order 1.\n", order);
    order = 1;
  }
#endif

  return write_gmsh(filename, order);
}

bool ExportGmshCommand::write_gmsh(const std::string &filename,
                                    int order)
{
  MeshData mesh;
  if (!mesh.extract(order))
    return false;

  bool ok = write_gmsh_v41(filename, mesh);

  if (ok)
    write_gmsh_launch_companions(filename);

  return ok;
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

  // Serendipity (edge-only HO nodes) types from GMSH API.
  // NetgenCurver produces serendipity elements.
  // WEDGE order 3-5: not supported by GMSH (FaceClosureFull not implemented).
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
    if (is_tet)  return 11;   // TET10
    if (is_hex)  return 17;   // HEX20
    if (is_wed)  return 18;   // WEDGE15 (Prism 15)
    if (is_pyr)  return 19;   // PYRAMID13
    if (is_tri)  return 9;    // TRI6
    if (is_quad) return 16;   // QUAD8
    if (is_line) return 8;    // LINE3
  }
  else if (order == 3) {
    if (is_tet)  return 29;   // TET20
    if (is_hex)  return 99;   // HEX32 (serendipity)
    if (is_pyr)  return 125;  // PYRAMID21 (serendipity)
    if (is_tri)  return 21;   // TRI10
    if (is_quad) return 39;   // QUAD12 (serendipity)
    if (is_line) return 26;   // LINE4
  }
  // Order 4-5: not supported for GMSH (face/volume node extraction unreliable).
  // Use export netgen for order 4-5.

  // Fallback: WEDGE order 3-5 not supported by GMSH, warn and use linear
  if (is_wed && order >= 3) {
    PRINT_WARNING("GMSH: Wedge/Prism order %d not supported (GMSH limitation). "
                  "Exporting as linear.\n", order);
  }
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
// Reorder HO edge nodes from Nastran canonical to GMSH file format order.
// GMSH .msh uses a different edge ordering than Nastran/VTK (see gmsh.info
// Section 10.2). Vertices (first nv nodes) are unchanged.
// ========================================================================
static std::vector<int> reorder_for_gmsh(const std::vector<int> &conn,
                                          int nv, ElementType type, int order,
                                          const MeshData *mesh = nullptr)
{
  if (order < 2 || (int)conn.size() <= nv)
    return conn;

  const int *reorder = nullptr;
  const bool *flip = nullptr;
  int n_edges = 0;  // number of edges (NOT number of HO edge nodes)

  if (nv == 4 && (type == TETRA4 || type == TETRA)) {
    reorder = EdgeTables::tet_gmsh_reorder; flip = EdgeTables::tet_gmsh_flip; n_edges = 6;
  } else if (nv == 8 && (type == HEX8 || type == HEX)) {
    reorder = EdgeTables::hex_gmsh_reorder; flip = EdgeTables::hex_gmsh_flip; n_edges = 12;
  } else if (nv == 6 && (type == WEDGE6 || type == WEDGE)) {
    reorder = EdgeTables::wedge_gmsh_reorder; flip = EdgeTables::wedge_gmsh_flip; n_edges = 9;
  } else if (nv == 5 && (type == PYRAMID5 || type == PYRAMID)) {
    reorder = EdgeTables::pyramid_gmsh_reorder; flip = EdgeTables::pyramid_gmsh_flip; n_edges = 8;
  } else if (nv == 3 && (type == TRI3 || type == CUBIT_TRI
                       || type == TRISHELL || type == TRISHELL3)) {
    reorder = EdgeTables::tri_gmsh_reorder; flip = EdgeTables::tri_gmsh_flip; n_edges = 3;
  } else if (nv == 4 && (type == QUAD4 || type == QUAD
                       || type == SHEL || type == SHELL4)) {
    reorder = EdgeTables::quad_gmsh_reorder; flip = EdgeTables::quad_gmsh_flip; n_edges = 4;
  }

  if (!reorder) return conn;

  int npe = order - 1;  // nodes per edge
  int n_total_edge_nodes = n_edges * npe;
  int n_mid = (int)conn.size() - nv;
  if (n_mid < n_total_edge_nodes) return conn;

  std::vector<int> out(conn.size());
  // Vertices
  for (int i = 0; i < nv; i++)
    out[i] = conn[i];
  // Edge nodes: reorder blocks of npe nodes per edge, with optional direction flip
  for (int e = 0; e < n_edges; e++) {
    int src_start = nv + e * npe;
    int dst_start = nv + reorder[e] * npe;
    if (flip && flip[e]) {
      // Reverse within-edge node order (GMSH edge direction opposite to internal)
      for (int k = 0; k < npe; k++)
        out[dst_start + k] = conn[src_start + (npe - 1 - k)];
    } else {
      for (int k = 0; k < npe; k++)
        out[dst_start + k] = conn[src_start + k];
    }
  }
  // Face + volume interior nodes: pass through for order <= 3.
  // Order 3 TET has 1 face node per face (centroid, no ordering issue).
  // Order 4-5 are not supported in GMSH export (clamped to 3 above).
  int face_start = nv + n_total_edge_nodes;
  for (int i = face_start; i < (int)conn.size(); i++)
    out[i] = conn[i];
  return out;
}

// ========================================================================
// gmsh_conn — get GMSH-ordered connectivity for a single element
// ========================================================================
static std::vector<int> gmsh_conn(const MeshElement &el, const MeshData &mesh)
{
  int order = mesh.order;
#ifdef HAVE_NETGEN
  auto nc = mesh.get_netgen_curver();
  if (order >= 2 && nc) {
    auto c = MeshData::build_ho_conn_gmsh(el.conn, el.nv, el.type, *nc, mesh);
    if (!c.empty()) return c;
  }
#endif
  // Fallback: use pre-built ho_conn (if any) with legacy reorder
  const auto &raw = (order >= 2 && !el.ho_conn.empty()) ? el.ho_conn : el.conn;
  return (order >= 2) ? reorder_for_gmsh(raw, el.nv, el.type, order, &mesh)
                      : std::vector<int>(raw.begin(), raw.end());
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
  std::ofstream fid(u8_string_to_path(filename));
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
    for (auto &[bid, bi] : block_info) {
      // Use the Cubit block name (material label) when present, matching
      // the netgen exporter; fall back to "block_<id>" only if unnamed.
      std::string bname = CubitInterface::get_block_name(bid);
      if (bname.empty())
        bname = "block_" + std::to_string(bid);
      fid << bi.dim << " " << bid << " \"" << bname << "\"\n";
    }
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
    blk.conns.push_back(gmsh_conn(el, mesh));
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
      blk.conns.push_back(gmsh_conn(face, mesh));
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

  return true;
}
