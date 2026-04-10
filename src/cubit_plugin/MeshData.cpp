#include "MeshData.hpp"
#ifdef HAVE_NETGEN
#include "NetgenCurver.hpp"
#endif
#include "MeshExportInterface.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"

#include <algorithm>
#include <set>
#include <sstream>

// ================================================================
// extract() — main entry point
// ================================================================
bool MeshData::extract(int req_order)
{
  MeshExportInterface *iface =
    dynamic_cast<MeshExportInterface*>(CubitInterface::get_interface("MeshExport"));
  if (!iface) {
    PRINT_ERROR("Unable to get MeshExport interface.\n");
    return false;
  }

  iface->set_use_sequential_ids(false);
  iface->initialize_export();

  if (iface->get_num_nodes() == 0) {
    PRINT_WARNING("No mesh nodes found. Nothing to export.\n");
    CubitInterface::release_interface(iface);
    return false;
  }

  order = req_order;

  // Extract linear mesh data first (single pass through iface)
  extract_nodes(iface);
  extract_elements(iface);
  extract_sidesets(iface);
  extract_nodesets(iface);

  // Build high-order from extracted MeshData — no more iface iteration
  if (order >= 2) {
#ifdef HAVE_NETGEN
    auto nc = std::make_shared<NetgenCurver>();
    PRINT_INFO("build_high_order: trying NetgenCurver (order %d)...\n", order);
    if (nc->build(*this, order)) {
      PRINT_INFO("build_high_order: NetgenCurver succeeded (%d nodes).\n",
                 nc->get_num_nodes());
      has_netgen = true;
      netgen_curver_ = nc;  // keep alive for pybind11 export
      apply_ho(*nc);
    } else {
      PRINT_ERROR("NetgenCurver failed for order %d. "
                  "No fallback — check mesh and geometry.\n", order);
      CubitInterface::release_interface(iface);
      return false;
#else
    {
      PRINT_ERROR("Order %d requires Netgen support (not built).\n", order);
      CubitInterface::release_interface(iface);
      return false;
#endif
    }
  }

  CubitInterface::release_interface(iface);
  return true;
}

// ================================================================
// extract_nodes
// ================================================================
void MeshData::extract_nodes(MeshExportInterface *iface)
{
  const int buf = 1000;
  std::vector<double> xc(buf), yc(buf), zc(buf);
  std::vector<int> ids(buf);

  int start = 0, found = 0;
  int idx = 0;
  while ((found = iface->get_coords(start, buf, xc, yc, zc, ids)) > 0) {
    for (int i = 0; i < found; i++) {
      nodes.push_back({ids[i], xc[i], yc[i], zc[i]});
      node_id_to_index[ids[i]] = idx++;
      if (ids[i] > max_original_node_id)
        max_original_node_id = ids[i];
    }
    start += found;
  }
  num_original_nodes = (int)nodes.size();
}

// ================================================================
// extract_elements — from blocks
// ================================================================
void MeshData::extract_elements(MeshExportInterface *iface)
{
  std::vector<BlockHandle> blocks;
  iface->get_block_list(blocks);

  std::set<int> seen_bids;
  // Deduplicate by sorted connectivity (ElementHandle may differ across blocks)
  std::set<std::vector<int>> seen_conn;

  const int ebuf = 100;
  std::vector<ElementType> etypes(ebuf);
  std::vector<ElementHandle> handles(ebuf);

  int dup_count = 0;
  for (size_t bi = 0; bi < blocks.size(); bi++) {
    BlockHandle block = blocks[bi];
    int bid = iface->id_from_handle(block);
    if (seen_bids.insert(bid).second)
      block_ids.push_back(bid);

    int estart = 0, count = 0;
    while ((count = iface->get_block_elements(estart, ebuf, block, etypes, handles)) > 0) {
      for (int i = 0; i < count; i++) {
        std::vector<int> conn(27);
        int nn = iface->get_connectivity(handles[i], conn);
        conn.resize(nn);

        // Skip duplicate elements (same connectivity in multiple blocks)
        std::vector<int> key(conn.begin(), conn.begin() + nn);
        std::sort(key.begin(), key.end());
        if (!seen_conn.insert(key).second) {
          dup_count++;
          continue;
        }

        MeshElement elem;
        elem.group_id = bid;
        elem.type = etypes[i];
        elem.nv = nn;
        elem.conn = std::move(conn);
        elements.push_back(std::move(elem));
      }
      estart += count;
    }
  }
  if (dup_count > 0)
    PRINT_WARNING("Skipped %d duplicate elements across blocks.\n", dup_count);
}

// ================================================================
// extract_sidesets — face elements from sideset surfaces
//
// ExodusII side table (get_sideset_sides + get_nodes_per_side) fails
// for wedge elements (returns nv=0). Instead, get the surfaces in
// each sideset and collect tri/face elements via parse_cubit_list,
// which works for all element types (tet, hex, wedge, pyramid).
// ================================================================
void MeshData::extract_sidesets(MeshExportInterface *iface)
{
  std::vector<SidesetHandle> ss_handles;
  iface->get_sideset_list(ss_handles);

  for (auto &ssh : ss_handles) {
    SidesetGroup sg;
    sg.id = iface->id_from_handle(ssh);
    sg.name = iface->get_sideset_name(ssh);

    // Get surfaces in this sideset
    std::vector<int> ss_surfs = CubitInterface::get_sideset_surfaces(sg.id);

    for (int sid : ss_surfs) {
      // Triangles on this surface (tet/wedge boundary faces)
      std::vector<int> tri_ids = CubitInterface::parse_cubit_list(
          "tri", "in surface " + std::to_string(sid));
      for (int tid : tri_ids) {
        std::vector<int> conn = CubitInterface::get_connectivity("tri", tid);
        if ((int)conn.size() < 3) continue;

        MeshElement face;
        face.group_id = sg.id;
        face.type = TRI3;
        face.nv = 3;
        face.conn = {conn[0], conn[1], conn[2]};
        sg.faces.push_back(std::move(face));
      }

      // Quads on this surface (hex/wedge boundary faces)
      std::vector<int> face_ids = CubitInterface::parse_cubit_list(
          "face", "in surface " + std::to_string(sid));
      for (int fid : face_ids) {
        std::vector<int> conn = CubitInterface::get_connectivity("face", fid);
        if ((int)conn.size() < 4) continue;

        MeshElement face;
        face.group_id = sg.id;
        face.type = QUAD4;
        face.nv = 4;
        face.conn = {conn[0], conn[1], conn[2], conn[3]};
        sg.faces.push_back(std::move(face));
      }
    }
    sidesets.push_back(std::move(sg));
  }
}

// ================================================================
// extract_nodesets
// ================================================================
void MeshData::extract_nodesets(MeshExportInterface *iface)
{
  std::vector<NodesetHandle> ns_handles;
  iface->get_nodeset_list(ns_handles);

  const int nbuf = 1000;
  std::vector<int> nids(nbuf);

  for (auto &nsh : ns_handles) {
    NodesetGroup ng;
    ng.id = iface->id_from_handle(nsh);
    ng.name = iface->get_nodeset_name(nsh);

    int nstart = 0, ncount = 0;
    while ((ncount = iface->get_nodeset_node_ids(nsh, nstart, nbuf, nids)) > 0) {
      for (int i = 0; i < ncount; i++)
        ng.node_ids.push_back(nids[i]);
      nstart += ncount;
    }
    nodesets.push_back(std::move(ng));
  }
}

// ================================================================
// apply_ho(NetgenCurver) — fill ho_conn + append generated nodes
// ================================================================
#ifdef HAVE_NETGEN
void MeshData::apply_ho(const NetgenCurver &nc)
{
  // Append generated HO nodes
  int n_ho = nc.get_num_nodes() - num_original_nodes;
  int idx = (int)nodes.size();
  for (int i = 0; i < n_ho; i++) {
    int nid = max_original_node_id + 1 + i;
    auto c = nc.get_node_coords(nid);
    nodes.push_back({nid, c[0], c[1], c[2]});
    node_id_to_index[nid] = idx++;
  }

  // Fill ho_conn for block elements
  for (auto &elem : elements)
    elem.ho_conn = build_ho_conn_nc(elem.conn, elem.nv, elem.type, nc);

  // Fill ho_conn for sideset faces
  for (auto &sg : sidesets)
    for (auto &face : sg.faces)
      face.ho_conn = build_ho_conn_nc(face.conn, face.nv, face.type, nc);
}

// ================================================================
// build_ho_conn_nc — assemble HO connectivity from NetgenCurver edges
// ================================================================
std::vector<int> MeshData::build_ho_conn_nc(
    const std::vector<int> &conn, int nv, ElementType type,
    const NetgenCurver &nc)
{
  // Select edge table
  const int (*edges)[2] = nullptr;
  int num_edges = 0;

  if (nv == 4 && (type == TETRA4 || type == TETRA)) {
    edges = EdgeTables::tet; num_edges = 6;
  } else if (nv == 8 && (type == HEX8 || type == HEX)) {
    edges = EdgeTables::hex; num_edges = 12;
  } else if (nv == 6 && (type == WEDGE6 || type == WEDGE)) {
    edges = EdgeTables::wedge; num_edges = 9;
  } else if (nv == 5 && (type == PYRAMID5 || type == PYRAMID)) {
    edges = EdgeTables::pyramid; num_edges = 8;
  } else if (nv == 3 && (type == TRI3 || type == CUBIT_TRI
                      || type == TRISHELL || type == TRISHELL3)) {
    edges = EdgeTables::tri; num_edges = 3;
  } else if (nv == 4 && (type == QUAD4 || type == QUAD
                      || type == SHEL || type == SHELL4)) {
    edges = EdgeTables::quad; num_edges = 4;
  } else {
    PRINT_WARNING("build_ho_conn_nc: unsupported type=%d nv=%d (element skipped)\n",
                  (int)type, nv);
    return {};  // unsupported element type
  }

  // Build: vertices + edge HO nodes (order-1 per edge)
  std::vector<int> ho;
  ho.reserve(nv + num_edges * (nc.get_order() - 1));

  // Vertices
  for (int j = 0; j < nv; j++)
    ho.push_back(conn[j]);

  // Edge nodes
#ifdef DEBUG_HO_NODES
  static int dbg_elem_count_ = 0;
  bool dbg_this = (dbg_elem_count_ < 3);
  dbg_elem_count_++;
#endif
  for (int e = 0; e < num_edges; e++) {
    auto en = nc.get_edge_nodes(conn[edges[e][0]], conn[edges[e][1]]);
    if (en.empty()) {
      // Edge not found in NetgenCurver — element not in curved mesh
      // (e.g., standalone tri block not part of any volume).
      // Fall back to linear connectivity.
      return {};
    }
#ifdef DEBUG_HO_NODES
    if (dbg_this && !en.empty()) {
      PRINT_INFO("  ho_conn edge %d: verts=(%d,%d) mid_id=%d\n",
        e, conn[edges[e][0]], conn[edges[e][1]], en[0]);
    }
#endif
    for (int nid : en)
      ho.push_back(nid);
    int expected = nc.get_order() - 1;
    for (int k = (int)en.size(); k < expected; k++)
      ho.push_back(en.back());  // repeat last node (safe fallback)
  }

  // Face interior nodes (order >= 3, TET and TRI only — HEX/QUAD/PYR use serendipity)
  if (nc.get_order() >= 3) {
    bool is_tet = (nv == 4 && (type == TETRA4 || type == TETRA));
    bool is_tri = (nv == 3 && (type == TRI3 || type == CUBIT_TRI
                             || type == TRISHELL || type == TRISHELL3));
    if (is_tet) {
      // GMSH TET face order: same vertex sets as ours, but winding may differ.
      // For order 3 (1 face node/face): position is unique, winding irrelevant.
      // For order 4+ (3+ face nodes/face): winding affects node ordering.
      // Our face order: (0,1,2), (0,1,3), (0,2,3), (1,2,3)
      // GMSH face order: (0,2,1), (0,1,3), (0,3,2), (3,1,2) [MTetrahedron.h]
      // get_face_nodes + permute_face_nodes handles vertex remapping,
      // but the output order depends on the generation winding.
      // We pass OUR face vertices — permute_face_nodes will reorder internally.
      static const int tet_faces[][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
      for (auto &f : tet_faces) {
        auto fn = nc.get_face_nodes(conn[f[0]], conn[f[1]], conn[f[2]]);
        for (int nid : fn)
          ho.push_back(nid);
      }
    }
    if (is_tri) {
      auto fn = nc.get_face_nodes(conn[0], conn[1], conn[2]);
      for (int nid : fn)
        ho.push_back(nid);
    }
  }

  // Volume interior nodes (order >= 4, TET only)
  if (nc.get_order() >= 4) {
    bool is_tet = (nv == 4 && (type == TETRA4 || type == TETRA));
    if (is_tet) {
      auto vn = nc.get_vol_nodes(conn[0], conn[1], conn[2], conn[3]);
      for (int nid : vn)
        ho.push_back(nid);
    }
  }

  return ho;
}
// ================================================================
// build_ho_conn_gmsh — GMSH-ordered connectivity in one pass
// ================================================================
std::vector<int> MeshData::build_ho_conn_gmsh(
    const std::vector<int> &conn, int nv, ElementType type,
    const NetgenCurver &nc, const MeshData &mesh)
{
  int order = nc.get_order();
  if (order < 2) return conn;

  bool is_tet = (nv == 4 && (type == TETRA4 || type == TETRA));
  bool is_tri = (nv == 3 && (type == TRI3 || type == CUBIT_TRI
                           || type == TRISHELL || type == TRISHELL3));

  if (is_tet) {
    // --- TET: build entirely in GMSH order ---
    // GMSH TET edges (gmsh.info Section 10.2, verified):
    //   pos 0: (0,1), pos 1: (1,2), pos 2: (2,0),
    //   pos 3: (3,0), pos 4: (3,2), pos 5: (3,1)
    static const int gmsh_edges[][2] = {
      {0,1},{1,2},{2,0},{3,0},{3,2},{3,1}
    };
    // GMSH TET faces (MTetrahedron.h, verified via Jacobian test):
    //   face 0: (0,2,1), face 1: (0,1,3), face 2: (0,3,2), face 3: (3,1,2)
    // get_face_nodes with GMSH winding returns face nodes directly in GMSH (i,j) order.
    static const int gmsh_faces[][3] = {
      {0,2,1},{0,1,3},{0,3,2},{3,1,2}
    };

    int npe = order - 1;  // nodes per edge
    int npf = (order - 1) * (order - 2) / 2;  // face interior per face
    int npv = (order - 1) * (order - 2) * (order - 3) / 6;  // volume interior
    int total = 4 + 6 * npe + 4 * npf + npv;

    std::vector<int> ho;
    ho.reserve(total);

    // Vertices (same order as GMSH)
    for (int j = 0; j < 4; j++) ho.push_back(conn[j]);

    // Edge nodes in GMSH order and direction
    for (int e = 0; e < 6; e++) {
      auto en = nc.get_edge_nodes(conn[gmsh_edges[e][0]],
                                   conn[gmsh_edges[e][1]]);
      if (en.empty()) return {};  // not in curved mesh
      for (int nid : en) ho.push_back(nid);
      for (int k = (int)en.size(); k < npe; k++)
        ho.push_back(en.back());  // safe fallback
    }

    // Face interior nodes in GMSH winding (order >= 3)
    if (order >= 3) {
      for (int f = 0; f < 4; f++) {
        auto fn = nc.get_face_nodes(conn[gmsh_faces[f][0]],
                                     conn[gmsh_faces[f][1]],
                                     conn[gmsh_faces[f][2]]);
        for (int nid : fn) ho.push_back(nid);
      }
    }

    // Volume interior nodes (order >= 4) — coordinate-based GMSH sorting
    if (order >= 4 && npv > 0) {
      auto vn = nc.get_vol_nodes(conn[0], conn[1], conn[2], conn[3]);
      if ((int)vn.size() >= npv) {
        // GMSH enumerates vol nodes as (i,j,k) for i=1..p-1, j=1..p-1-i, k=1..p-1-i-j
        // where i = lam1 * p (weight on v1), j = lam2 * p, k = lam3 * p.
        std::vector<std::array<int,3>> gmsh_ijk;
        for (int i = 1; i < order; i++)
          for (int j = 1; j < order - i; j++)
            for (int k = 1; k < order - i - j; k++)
              gmsh_ijk.push_back({i, j, k});

        // Compute barycentric (i,j,k) for each vol node via Cramer's rule
        auto v0c = mesh.get_node_coords(conn[0]);
        auto v1c = mesh.get_node_coords(conn[1]);
        auto v2c = mesh.get_node_coords(conn[2]);
        auto v3c = mesh.get_node_coords(conn[3]);

        double e1[3], e2[3], e3[3];
        for (int d = 0; d < 3; d++) {
          e1[d] = v1c[d] - v0c[d];
          e2[d] = v2c[d] - v0c[d];
          e3[d] = v3c[d] - v0c[d];
        }
        double det_A = e1[0]*(e2[1]*e3[2]-e2[2]*e3[1])
                     - e1[1]*(e2[0]*e3[2]-e2[2]*e3[0])
                     + e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]);

        struct VolNode { int node_id; int i, j, k; };
        std::vector<VolNode> vol_nodes(npv);
        for (int n = 0; n < npv; n++) {
          auto nc_coord = mesh.get_node_coords(vn[n]);
          double rhs[3] = {nc_coord[0]-v0c[0], nc_coord[1]-v0c[1], nc_coord[2]-v0c[2]};
          double u = (rhs[0]*(e2[1]*e3[2]-e2[2]*e3[1])
                    - rhs[1]*(e2[0]*e3[2]-e2[2]*e3[0])
                    + rhs[2]*(e2[0]*e3[1]-e2[1]*e3[0])) / det_A;
          double v = (e1[0]*(rhs[1]*e3[2]-rhs[2]*e3[1])
                    - e1[1]*(rhs[0]*e3[2]-rhs[2]*e3[0])
                    + e1[2]*(rhs[0]*e3[1]-rhs[1]*e3[0])) / det_A;
          double w = (e1[0]*(e2[1]*rhs[2]-e2[2]*rhs[1])
                    - e1[1]*(e2[0]*rhs[2]-e2[2]*rhs[0])
                    + e1[2]*(e2[0]*rhs[1]-e2[1]*rhs[0])) / det_A;
          vol_nodes[n] = {vn[n], (int)std::round(u * order),
                                 (int)std::round(v * order),
                                 (int)std::round(w * order)};
        }

        // Place each vol node at the GMSH target position
        for (int g = 0; g < npv; g++) {
          auto &target = gmsh_ijk[g];
          bool found = false;
          for (int n = 0; n < npv; n++) {
            if (vol_nodes[n].i == target[0] &&
                vol_nodes[n].j == target[1] &&
                vol_nodes[n].k == target[2]) {
              ho.push_back(vol_nodes[n].node_id);
              found = true;
              break;
            }
          }
          if (!found) ho.push_back(vn[g]);  // fallback
        }
      } else {
        for (int nid : vn) ho.push_back(nid);
      }
    }

    return ho;
  }

  // --- Non-TET: build with internal order, then reorder edges ---
  // Non-TET elements use serendipity (edge-only HO), so no face/vol reordering.
  // TRI/QUAD: GMSH edge order matches internal (identity reorder).
  auto ho = build_ho_conn_nc(conn, nv, type, nc);
  if (ho.empty()) return {};

  // Select reorder and flip tables
  const int *reorder = nullptr;
  const bool *flip = nullptr;
  int n_edges = 0;

  if (nv == 8 && (type == HEX8 || type == HEX)) {
    reorder = EdgeTables::hex_gmsh_reorder; flip = EdgeTables::hex_gmsh_flip; n_edges = 12;
  } else if (nv == 6 && (type == WEDGE6 || type == WEDGE)) {
    reorder = EdgeTables::wedge_gmsh_reorder; flip = EdgeTables::wedge_gmsh_flip; n_edges = 9;
  } else if (nv == 5 && (type == PYRAMID5 || type == PYRAMID)) {
    reorder = EdgeTables::pyramid_gmsh_reorder; flip = EdgeTables::pyramid_gmsh_flip; n_edges = 8;
  } else if (is_tri) {
    // TRI: identity reorder, no flip — return as-is
    return ho;
  } else if (nv == 4 && (type == QUAD4 || type == QUAD
                       || type == SHEL || type == SHELL4)) {
    // QUAD: identity reorder — return as-is
    return ho;
  }

  if (!reorder) return ho;

  int npe = order - 1;
  int n_total_edge_nodes = n_edges * npe;
  int n_mid = (int)ho.size() - nv;
  if (n_mid < n_total_edge_nodes) return ho;

  std::vector<int> out(ho.size());
  // Vertices unchanged
  for (int i = 0; i < nv; i++) out[i] = ho[i];

  // Reorder edge nodes
  for (int e = 0; e < n_edges; e++) {
    int src_start = nv + e * npe;
    int dst_start = nv + reorder[e] * npe;
    if (flip && flip[e]) {
      for (int k = 0; k < npe; k++)
        out[dst_start + k] = ho[src_start + (npe - 1 - k)];
    } else {
      for (int k = 0; k < npe; k++)
        out[dst_start + k] = ho[src_start + k];
    }
  }

  // Copy remaining nodes (face interior for TRI — but TRI already returned above)
  int edge_end = nv + n_total_edge_nodes;
  for (int i = edge_end; i < (int)ho.size(); i++)
    out[i] = ho[i];

  return out;
}
#endif // HAVE_NETGEN

// ================================================================
// total_element_count
// ================================================================
int MeshData::total_element_count() const
{
  int total = (int)elements.size();
  for (auto &sg : sidesets)
    total += (int)sg.faces.size();
  return total;
}

