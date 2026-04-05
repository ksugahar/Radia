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
// extract_sidesets — face elements from sideset side definitions
// ================================================================
void MeshData::extract_sidesets(MeshExportInterface *iface)
{
  std::vector<SidesetHandle> ss_handles;
  iface->get_sideset_list(ss_handles);

  const int sbuf = 100;
  std::vector<ElementHandle> ehandles(sbuf);
  std::vector<ElementType> etypes(sbuf);
  std::vector<int> side_ids(sbuf);

  for (auto &ssh : ss_handles) {
    SidesetGroup sg;
    sg.id = iface->id_from_handle(ssh);
    sg.name = iface->get_sideset_name(ssh);

    int sstart = 0, scount = 0;
    while ((scount = iface->get_sideset_sides(ssh, sstart, sbuf,
                                               ehandles, etypes, side_ids)) > 0) {
      for (int i = 0; i < scount; i++) {
        // Get parent element connectivity
        std::vector<int> parent_conn(27);
        int parent_nn = iface->get_connectivity(ehandles[i], parent_conn);
        (void)parent_nn;

        // Get face node indices from ExodusII side table
        int side_nn = 0;
        iface->get_nodes_per_side(etypes[i], side_ids[i], side_nn);
        std::vector<int> idx(side_nn);
        iface->get_node_indices_per_side(etypes[i], side_ids[i], idx);

        // Build face element
        MeshElement face;
        face.group_id = sg.id;
        face.nv = side_nn;
        face.conn.resize(side_nn);
        for (int j = 0; j < side_nn; j++)
          face.conn[j] = parent_conn[idx[j]];

        // Determine face element type
        if (side_nn == 3)
          face.type = TRI3;
        else if (side_nn == 4)
          face.type = QUAD4;
        else
          face.type = INVALID_ELEMENT_TYPE;

        sg.faces.push_back(std::move(face));
      }
      sstart += scount;
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
  for (int e = 0; e < num_edges; e++) {
    auto en = nc.get_edge_nodes(conn[edges[e][0]], conn[edges[e][1]]);
    if (en.empty()) {
      // Edge not found in NetgenCurver — element not in curved mesh
      // (e.g., standalone tri block not part of any volume).
      // Fall back to linear connectivity.
      return {};
    }
    for (int nid : en)
      ho.push_back(nid);
    int expected = nc.get_order() - 1;
    for (int k = (int)en.size(); k < expected; k++)
      ho.push_back(en.back());  // repeat last node (safe fallback)
  }

  // TODO: face interior nodes and volume interior nodes for order >= 3
  // These are needed for Gmsh TET20, TET35, etc.
  // Currently NetgenCurver generates them but doesn't expose per-element.

  return ho;
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

