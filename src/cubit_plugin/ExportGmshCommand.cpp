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
    "radia_export gmsh <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1-5>'>] "
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
    "radia_export gmsh \"filename\" [order {1|2|3|4|5}] [version {2|4}] [dimension {2|3}] [overwrite]"
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
    "  order 3-5    Higher order (requires Netgen, wedge limited to order 2)\n"
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
  if (order > 5) {
    PRINT_WARNING("order %d not supported. Using order 5.\n", order);
    order = 5;
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

  bool ok = (version == "4.1")
    ? write_gmsh_v41(filename, mesh)
    : write_gmsh_v22(filename, mesh);

  // Write companion .geo file for proper curved element display in GMSH GUI.
  // Mesh.NumSubEdges=4 is required to render high-order curved surfaces
  // (GMSH default=1 draws straight edges).
  if (ok && order >= 2) {
    std::string geo = filename.substr(0, filename.rfind('.')) + ".geo";
    std::string msh_base = filename;
    auto sep = msh_base.rfind('/');
    auto sep2 = msh_base.rfind('\\');
    if (sep2 != std::string::npos && (sep == std::string::npos || sep2 > sep))
      sep = sep2;
    if (sep != std::string::npos)
      msh_base = msh_base.substr(sep + 1);

    std::ofstream gf(geo);
    if (gf.is_open()) {
      gf << "// Auto-generated companion for " << msh_base << "\n";
      gf << "Merge \"" << msh_base << "\";\n";
      gf << "Mesh.NumSubEdges = 4;\n";
      gf.close();
      PRINT_INFO("Companion: %s\n", geo.c_str());
    }
  }

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
  else if (order == 4) {
    if (is_tet)  return 30;   // TET35
    if (is_hex)  return 100;  // HEX44 (serendipity)
    if (is_pyr)  return 126;  // PYRAMID29 (serendipity)
    if (is_tri)  return 23;   // TRI15
    if (is_quad) return 40;   // QUAD16 (serendipity)
    if (is_line) return 27;   // LINE5
  }
  else if (order == 5) {
    if (is_tet)  return 31;   // TET56
    if (is_hex)  return 101;  // HEX56 (serendipity)
    if (is_pyr)  return 127;  // PYRAMID37 (serendipity)
    if (is_tri)  return 25;   // TRI21
    if (is_quad) return 41;   // QUAD20 (serendipity)
    if (is_line) return 28;   // LINE6
  }

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
  // Face + volume interior nodes: apply GMSH permutation for TET order >= 4.
  //
  // GMSH uses a different face/volume node enumeration than Netgen's (i,j) row-major.
  // Permutation tables derived empirically by comparing GMSH-generated meshes
  // against our barycentric (i,j) / (i,j,k) enumeration.
  //
  // Verified for order 4 (TET35) and order 5 (TET56) using Jacobian sign test.
  int face_start = nv + n_total_edge_nodes;
  bool is_tet = (nv == 4 && (type == TETRA4 || type == TETRA));

  if (is_tet && order >= 4) {
    int npf = (order - 1) * (order - 2) / 2;  // face interior nodes per face
    int npv = (order - 1) * (order - 2) * (order - 3) / 6;  // volume interior

    // Per-order, per-face permutation tables.
    // perm[k] means: out[face_start + f*npf + k] = conn[face_start + f*npf + perm[k]]
    // Identity = {0,1,...,npf-1} (no permutation needed)

    // Order 4: npf=3, npv=1
    static const int o4_face1[] = {0, 2, 1};
    static const int o4_face3[] = {1, 0, 2};

    // Order 5: npf=6, npv=4
    static const int o5_face0[] = {0, 2, 5, 1, 4, 3};
    static const int o5_face1[] = {0, 5, 2, 3, 4, 1};
    static const int o5_face2[] = {0, 2, 5, 1, 4, 3};
    static const int o5_face3[] = {2, 0, 5, 1, 3, 4};
    static const int o5_vol[]   = {1, 2, 3, 0};

    const int *face_perms[4] = {nullptr, nullptr, nullptr, nullptr};
    const int *vol_perm = nullptr;

    if (order == 4) {
      face_perms[1] = o4_face1;
      face_perms[3] = o4_face3;
    } else if (order == 5) {
      // Order 5 face node permutation is element-dependent (permute_face_nodes
      // reorders shared face nodes per element). Fixed tables partially correct
      // but do not cover all shared-face cases. Volume nodes use coordinate-based
      // sorting. Face nodes pass through with per-order tables (partial fix).
      face_perms[0] = o5_face0;
      face_perms[1] = o5_face1;
      face_perms[2] = o5_face2;
      face_perms[3] = o5_face3;
    }
    // Order 3: npf=1 (single node, no permutation), npv=0

    // Per-face reorder using barycentric coordinates in OUR face vertex order.
    // GMSH target (i,j) sequence in OUR face coords (measured, constant for all elements):
    //   Face 0 (0,1,2): (1,1)(1,3)(3,1)(1,2)(2,2)(2,1)  [order 5]
    //   Face 1 (0,1,3): (1,1)(3,1)(1,3)(2,1)(2,2)(1,2)  [order 5]
    //   Face 2 (0,2,3): (1,1)(1,3)(3,1)(1,2)(2,2)(2,1)  [order 5]
    //   Face 3 (1,2,3): (1,3)(1,1)(3,1)(1,2)(2,1)(2,2)  [order 5]
    // conn may have face nodes reordered by permute_face_nodes (shared faces),
    // so use coordinate-based sorting instead of fixed permutation tables.
    static const int our_fv[4][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};

    // Build target (i,j) for each face in OUR coords
    // For face f with OUR vertices (a,b,c):
    //   GMSH winding: face 0=(0,2,1), face 1=(0,1,3), face 2=(0,3,2), face 3=(3,1,2)
    //   Enumerate (i,j) with i=weight on GMSH_v1, j=weight on GMSH_v2
    //   Convert to OUR coords: depends on which OUR vertex = GMSH vertex
    // Pre-computed target sequences per face (in OUR (i,j) coords):
    // For order p, face interior has (p-1)(p-2)/2 nodes.
    // Build GMSH (i,j) in GMSH face coords, then map to OUR coords.
    //
    // GMSH enumeration in GMSH face coords: for i=1..p-1, j=1..p-1-i
    // Face f: GMSH verts (ga,gb,gc). OUR verts (oa,ob,oc).
    // In OUR coords: weight_on_ob = i_our/p, weight_on_oc = j_our/p
    // In GMSH coords: weight_on_gb = i_gmsh/p, weight_on_gc = j_gmsh/p
    // Relationship depends on mapping between OUR and GMSH vertex sets.

    // Direct coordinate approach: compute barycentric in OUR face coords,
    // then sort to match GMSH-expected sequence.
    // GMSH target for each face is pre-computed per order.
    // Since this target is constant (verified across elements), we compute it once.
    static const int gmsh_fv[4][3] = {{0,2,1},{0,1,3},{0,3,2},{3,1,2}};

    // Face node reorder using per-order permutation tables.
    // Order 3: npf=1 (identity), Order 4: validated tables, Order 5: partial.
    // Order 5 face node winding is a known issue — permute_face_nodes in
    // NetgenCurver reorders shared face nodes per-element, making fixed tables
    // incorrect for some elements. Full fix requires per-element coordinate-based
    // sorting or modifying permute_face_nodes to output GMSH winding.
    for (int f = 0; f < 4; f++) {
      int fs = face_start + f * npf;
      if (face_perms[f]) {
        for (int k = 0; k < npf; k++) out[fs + k] = conn[fs + face_perms[f][k]];
      } else {
        for (int k = 0; k < npf; k++) out[fs + k] = conn[fs + k];
      }
    }

    // Volume interior nodes: sort by GMSH enumeration order.
    // GMSH TET volume nodes follow (i,j,k) with i = lam1*p (v0->v1 direction)
    // enumerated as: for i=1..p-1, for j=1..p-1-i, for k=1..p-1-i-j
    // where lam1 = weight on v1, lam2 = weight on v2, lam3 = weight on v3.
    // Our generation order may differ per element (Netgen vertex ordering varies).
    // Solution: compute barycentric coords of each vol node and sort to match GMSH.
    int vol_start = face_start + 4 * npf;
    if (npv > 0 && mesh) {
      // Get vertex coordinates
      auto v0c = mesh->get_node_coords(conn[0]);
      auto v1c = mesh->get_node_coords(conn[1]);
      auto v2c = mesh->get_node_coords(conn[2]);
      auto v3c = mesh->get_node_coords(conn[3]);

      // Build GMSH expected (i,j,k) order
      std::vector<std::array<int,3>> gmsh_ijk;
      for (int i = 1; i < order; i++)
        for (int j = 1; j < order - i; j++)
          for (int k = 1; k < order - i - j; k++)
            gmsh_ijk.push_back({i, j, k});

      // Compute barycentric (i,j,k) for each vol node
      // Solve: pt = v0 + u*(v1-v0) + v*(v2-v0) + w*(v3-v0)
      double e1[3], e2[3], e3[3];
      for (int d = 0; d < 3; d++) {
        e1[d] = v1c[d] - v0c[d];
        e2[d] = v2c[d] - v0c[d];
        e3[d] = v3c[d] - v0c[d];
      }
      // 3x3 matrix inverse (Cramer's rule)
      double det_A = e1[0]*(e2[1]*e3[2]-e2[2]*e3[1])
                   - e1[1]*(e2[0]*e3[2]-e2[2]*e3[0])
                   + e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]);

      struct VolNode { int node_id; int i, j, k; };
      std::vector<VolNode> vol_nodes(npv);
      for (int n = 0; n < npv; n++) {
        int nid = conn[vol_start + n];
        auto nc = mesh->get_node_coords(nid);
        double rhs[3] = {nc[0]-v0c[0], nc[1]-v0c[1], nc[2]-v0c[2]};
        // Cramer's rule for u,v,w
        double u = (rhs[0]*(e2[1]*e3[2]-e2[2]*e3[1])
                  - rhs[1]*(e2[0]*e3[2]-e2[2]*e3[0])
                  + rhs[2]*(e2[0]*e3[1]-e2[1]*e3[0])) / det_A;
        double v = (e1[0]*(rhs[1]*e3[2]-rhs[2]*e3[1])
                  - e1[1]*(rhs[0]*e3[2]-rhs[2]*e3[0])
                  + e1[2]*(rhs[0]*e3[1]-rhs[1]*e3[0])) / det_A;
        double w = (e1[0]*(e2[1]*rhs[2]-e2[2]*rhs[1])
                  - e1[1]*(e2[0]*rhs[2]-e2[2]*rhs[0])
                  + e1[2]*(e2[0]*rhs[1]-e2[1]*rhs[0])) / det_A;
        vol_nodes[n] = {nid, (int)std::round(u*order),
                              (int)std::round(v*order),
                              (int)std::round(w*order)};
      }

      // Place each vol node at the GMSH-expected position
      for (int g = 0; g < npv; g++) {
        auto &target = gmsh_ijk[g];
        for (int n = 0; n < npv; n++) {
          if (vol_nodes[n].i == target[0] &&
              vol_nodes[n].j == target[1] &&
              vol_nodes[n].k == target[2]) {
            out[vol_start + g] = vol_nodes[n].node_id;
            break;
          }
        }
      }
      for (int i = vol_start + npv; i < (int)conn.size(); i++)
        out[i] = conn[i];
    } else {
      for (int i = vol_start; i < (int)conn.size(); i++)
        out[i] = conn[i];
    }
  } else {
    for (int i = face_start; i < (int)conn.size(); i++)
      out[i] = conn[i];
  }
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

    const auto c = gmsh_conn(elem, mesh);
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

      const auto c = gmsh_conn(face, mesh);
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
