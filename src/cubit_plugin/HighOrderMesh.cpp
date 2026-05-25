#include "HighOrderMesh.hpp"
#include "MeshExportInterface.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"

// ACIS geometry access for direct projection
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "CubitVector.hpp"
// 2026-05-25: Cubit 2025.12 removed GeometryQueryTool::instance() static
// accessor; CubitRefAccessor wraps
//   CubitCoreModel::instance()->cgmapp()->factory()->get_ref_{face,edge}()
// so the 2 lookup sites below stay one-liners.
#include "CubitRefAccessor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>

HighOrderMesh::HighOrderMesh() {}

bool HighOrderMesh::build(MeshExportInterface *iface)
{
  // -------------------------------------------------------
  // Step 1: Cache all original node coordinates
  // -------------------------------------------------------
  int num_nodes = iface->get_num_nodes();
  if (num_nodes == 0) return false;

  const int buf = 1000;
  std::vector<double> xc(buf), yc(buf), zc(buf);
  std::vector<int> ids(buf);
  int start = 0, found = 0;
  int max_id = 0;

  while ((found = iface->get_coords(start, buf, xc, yc, zc, ids)) > 0) {
    for (int i = 0; i < found; i++) {
      // Note: we don't store orig_coords separately; we use CubitInterface
      // to query as needed. But we do track max_id for new node numbering.
      if (ids[i] > max_id) max_id = ids[i];
    }
    start += found;
  }

  next_node_id_ = max_id + 1;
  total_nodes_ = num_nodes;  // will grow as mid-nodes are created

  // -------------------------------------------------------
  // Step 2: Classify nodes by geometric owner
  // Suppress "Only bodies currently accepted" errors from
  // CubitInterface geometry queries on mesh nodes.
  // -------------------------------------------------------
  // Suppress "Only bodies currently accepted" errors from geometry queries
  CubitMessage::instance()->set_error_flag(false);
  classify_nodes();

  // -------------------------------------------------------
  // Step 3: Pre-scan all elements to create mid-nodes
  // -------------------------------------------------------
  std::vector<BlockHandle> blocks;
  iface->get_block_list(blocks);

  const int ebuf = 100;
  std::vector<ElementType> etypes(ebuf);
  std::vector<ElementHandle> handles(ebuf);

  for (size_t bi = 0; bi < blocks.size(); bi++) {
    BlockHandle block = blocks[bi];
    int estart = 0, count = 0;
    while ((count = iface->get_block_elements(estart, ebuf, block, etypes, handles)) > 0) {
      for (int i = 0; i < count; i++) {
        std::vector<int> conn(27);
        int nn = iface->get_connectivity(handles[i], conn);
        ElementType et = etypes[i];

        // Create mid-nodes for all edges of this element
        if (nn == 4 && (et == TETRA4 || et == TETRA)) {
          // Tet4: 6 edges
          static const int edges[][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
          for (auto &e : edges) get_or_create_midnode(conn[e[0]], conn[e[1]]);
        }
        else if (nn == 8 && (et == HEX8 || et == HEX)) {
          // Hex8: 12 edges
          static const int edges[][2] = {
            {0,1},{1,2},{2,3},{3,0}, {4,5},{5,6},{6,7},{7,4},
            {0,4},{1,5},{2,6},{3,7}
          };
          for (auto &e : edges) get_or_create_midnode(conn[e[0]], conn[e[1]]);
        }
        else if (nn == 6 && (et == WEDGE6 || et == WEDGE)) {
          // Wedge6: 9 edges
          static const int edges[][2] = {
            {0,1},{1,2},{2,0}, {3,4},{4,5},{5,3}, {0,3},{1,4},{2,5}
          };
          for (auto &e : edges) get_or_create_midnode(conn[e[0]], conn[e[1]]);
        }
        else if (nn == 5 && (et == PYRAMID5 || et == PYRAMID)) {
          // Pyramid5: 8 edges
          static const int edges[][2] = {
            {0,1},{1,2},{2,3},{3,0}, {0,4},{1,4},{2,4},{3,4}
          };
          for (auto &e : edges) get_or_create_midnode(conn[e[0]], conn[e[1]]);
        }
        else if (nn == 3 && (et == TRI3 || et == CUBIT_TRI
                         || et == TRISHELL || et == TRISHELL3)) {
          static const int edges[][2] = {{0,1},{1,2},{2,0}};
          for (auto &e : edges) get_or_create_midnode(conn[e[0]], conn[e[1]]);
        }
        else if (nn == 4 && (et == QUAD4 || et == QUAD
                         || et == SHEL || et == SHELL4)) {
          static const int edges[][2] = {{0,1},{1,2},{2,3},{3,0}};
          for (auto &e : edges) get_or_create_midnode(conn[e[0]], conn[e[1]]);
        }
      }
      estart += count;
    }
  }

  CubitMessage::instance()->set_error_flag(true);

  PRINT_INFO("High-order mesh: %d original nodes + %d mid-nodes = %d total\n",
             num_nodes, total_nodes_ - num_nodes, total_nodes_);
  return true;
}

// ============================================================
// Node classification using CubitInterface geometry queries
// ============================================================
void HighOrderMesh::classify_nodes()
{
  // Get all geometry entity IDs
  std::vector<int> vertex_ids = CubitInterface::get_entities("vertex");
  std::vector<int> curve_ids  = CubitInterface::get_entities("curve");
  std::vector<int> surface_ids = CubitInterface::get_entities("surface");

  CubitMessage::instance()->set_error_flag(true);  // temporarily restore for debug
  PRINT_INFO("classify_nodes: %d vertices, %d curves, %d surfaces\n",
             (int)vertex_ids.size(), (int)curve_ids.size(),
             (int)surface_ids.size());
  CubitMessage::instance()->set_error_flag(false);  // re-suppress

  // Vertex nodes
  for (int vid : vertex_ids) {
    int nid = CubitInterface::get_vertex_node(vid);
    if (nid > 0) {
      node_owner_[nid] = {GeomOwner::VERTEX, vid};
      // A vertex belongs to its parent curves
      std::vector<int> parent_curves =
        CubitInterface::get_relatives("vertex", vid, "curve");
      node_curves_[nid] = parent_curves;
      // And parent surfaces
      std::vector<int> parent_surfaces =
        CubitInterface::get_relatives("vertex", vid, "surface");
      node_surfaces_[nid] = parent_surfaces;
    }
  }

  // Curve nodes (excludes vertex nodes)
  for (int cid : curve_ids) {
    std::vector<int> cnodes = CubitInterface::get_curve_nodes(cid);
    // Parent surfaces of this curve
    std::vector<int> parent_surfaces =
      CubitInterface::get_relatives("curve", cid, "surface");
    for (int nid : cnodes) {
      node_owner_[nid] = {GeomOwner::CURVE, cid};
      node_curves_[nid].push_back(cid);
      for (int sid : parent_surfaces) {
        node_surfaces_[nid].push_back(sid);
      }
    }
  }

  // Surface nodes (excludes curve and vertex nodes)
  for (int sid : surface_ids) {
    std::vector<int> snodes = CubitInterface::get_surface_nodes(sid);
    for (int nid : snodes) {
      node_owner_[nid] = {GeomOwner::SURFACE, sid};
      node_surfaces_[nid].push_back(sid);
    }
  }

  // Any node not in the map is a volume-interior node
  CubitMessage::instance()->set_error_flag(true);
  PRINT_INFO("classify_nodes: %d nodes classified (curve: %d, surface: %d)\n",
             (int)node_owner_.size(), (int)node_curves_.size(),
             (int)node_surfaces_.size());
  CubitMessage::instance()->set_error_flag(false);
}

// ============================================================
// Get or create a mid-node for edge (n0, n1)
// ============================================================
int HighOrderMesh::get_or_create_midnode(int n0, int n1)
{
  EdgeKey key(n0, n1);
  auto it = edge_midnodes_.find(key);
  if (it != edge_midnodes_.end())
    return it->second;

  // Get endpoint coordinates
  auto c0 = CubitInterface::get_nodal_coordinates(n0);
  auto c1 = CubitInterface::get_nodal_coordinates(n1);

  // Linear midpoint
  double mx = 0.5 * (c0[0] + c1[0]);
  double my = 0.5 * (c0[1] + c1[1]);
  double mz = 0.5 * (c0[2] + c1[2]);

  // Determine projection target using ACIS geometry directly
  // Priority: curve > surface > none (interior stays linear)

  // Check if both nodes share a common curve
  auto ci0 = node_curves_.find(n0);
  auto ci1 = node_curves_.find(n1);
  if (ci0 != node_curves_.end() && ci1 != node_curves_.end()) {
    // Find intersection of curve sets
    for (int c : ci0->second) {
      for (int d : ci1->second) {
        if (c == d) {
          // Both on same curve → project to curve via RefEdge
          RefEdge *re = radia::cubit_get_ref_edge(c);
          if (re) {
            CubitVector input_pt(mx, my, mz);
            CubitVector closest_pt;
            if (re->closest_point_trimmed(input_pt, closest_pt) == CUBIT_SUCCESS) {
              mx = closest_pt.x();
              my = closest_pt.y();
              mz = closest_pt.z();
            }
          }
          goto done;
        }
      }
    }
  }

  // Check if both nodes share a common surface
  {
    auto si0 = node_surfaces_.find(n0);
    auto si1 = node_surfaces_.find(n1);
    if (si0 != node_surfaces_.end() && si1 != node_surfaces_.end()) {
      for (int s : si0->second) {
        for (int t : si1->second) {
          if (s == t) {
            // Both on same surface → project to surface via RefFace
            RefFace *rf = radia::cubit_get_ref_face(s);
            if (rf) {
              CubitVector input_pt(mx, my, mz);
              CubitVector closest_pt;
              rf->find_closest_point_trimmed(input_pt, closest_pt);
              double ox = mx, oy = my, oz = mz;
              mx = closest_pt.x();
              my = closest_pt.y();
              mz = closest_pt.z();
              {
                static int s_count = 0;
                s_count++;
                if (s_count <= 3) {
                  PRINT_INFO("snap_surface[%d]: surf=%d "
                             "before=(%.6f,%.6f,%.6f) after=(%.6f,%.6f,%.6f)\n",
                             s_count, s,
                             ox, oy, oz, mx, my, mz);
                }
              }
            }
            goto done;
          }
        }
      }
    }
  }

done:
  int mid_id = next_node_id_++;
  total_nodes_++;
  edge_midnodes_[key] = mid_id;
  midnode_coords_[mid_id] = {mx, my, mz};
  return mid_id;
}

// ============================================================
// Get coordinates for any node (original or generated)
// ============================================================
std::array<double,3> HighOrderMesh::get_node_coords(int node_id) const
{
  // Check generated mid-nodes first
  auto it = midnode_coords_.find(node_id);
  if (it != midnode_coords_.end())
    return it->second;

  // Original node: query CubitInterface
  return CubitInterface::get_nodal_coordinates(node_id);
}

// ============================================================
// Element conversion functions
//
// Nastran 2nd-order node ordering:
//   CTETRA(10): v0 v1 v2 v3  m01 m12 m02 m03 m13 m23
//   CHEXA(20):  v0..v7  m01 m12 m23 m30  m45 m56 m67 m74  m04 m15 m26 m37
//   CPENTA(15): v0..v5  m01 m12 m20  m34 m45 m53  m03 m14 m25
//   CTRIA6:     v0 v1 v2  m01 m12 m20
//   CQUAD8:     v0 v1 v2 v3  m01 m12 m23 m30
//   CPYRAM(13): v0..v4  m01 m12 m23 m30  m04 m14 m24 m34
// ============================================================

// Helper to look up an existing mid-node
int HighOrderMesh_lookup(const std::unordered_map<EdgeKey, int, EdgeKeyHash> &table,
                          int a, int b)
{
  auto it = table.find(EdgeKey(a, b));
  return (it != table.end()) ? it->second : -1;
}

std::vector<int> HighOrderMesh::tet4_to_tet10(const int *c) const
{
  return {
    c[0], c[1], c[2], c[3],
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[1]),  // m01
    HighOrderMesh_lookup(edge_midnodes_, c[1], c[2]),  // m12
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[2]),  // m02
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[3]),  // m03
    HighOrderMesh_lookup(edge_midnodes_, c[1], c[3]),  // m13
    HighOrderMesh_lookup(edge_midnodes_, c[2], c[3])   // m23
  };
}

std::vector<int> HighOrderMesh::hex8_to_hex20(const int *c) const
{
  return {
    c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7],
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[1]),  // m01
    HighOrderMesh_lookup(edge_midnodes_, c[1], c[2]),  // m12
    HighOrderMesh_lookup(edge_midnodes_, c[2], c[3]),  // m23
    HighOrderMesh_lookup(edge_midnodes_, c[3], c[0]),  // m30
    HighOrderMesh_lookup(edge_midnodes_, c[4], c[5]),  // m45
    HighOrderMesh_lookup(edge_midnodes_, c[5], c[6]),  // m56
    HighOrderMesh_lookup(edge_midnodes_, c[6], c[7]),  // m67
    HighOrderMesh_lookup(edge_midnodes_, c[7], c[4]),  // m74
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[4]),  // m04
    HighOrderMesh_lookup(edge_midnodes_, c[1], c[5]),  // m15
    HighOrderMesh_lookup(edge_midnodes_, c[2], c[6]),  // m26
    HighOrderMesh_lookup(edge_midnodes_, c[3], c[7])   // m37
  };
}

std::vector<int> HighOrderMesh::wedge6_to_wedge15(const int *c) const
{
  return {
    c[0], c[1], c[2], c[3], c[4], c[5],
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[1]),  // m01
    HighOrderMesh_lookup(edge_midnodes_, c[1], c[2]),  // m12
    HighOrderMesh_lookup(edge_midnodes_, c[2], c[0]),  // m20
    HighOrderMesh_lookup(edge_midnodes_, c[3], c[4]),  // m34
    HighOrderMesh_lookup(edge_midnodes_, c[4], c[5]),  // m45
    HighOrderMesh_lookup(edge_midnodes_, c[5], c[3]),  // m53
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[3]),  // m03
    HighOrderMesh_lookup(edge_midnodes_, c[1], c[4]),  // m14
    HighOrderMesh_lookup(edge_midnodes_, c[2], c[5])   // m25
  };
}

std::vector<int> HighOrderMesh::tri3_to_tri6(const int *c) const
{
  return {
    c[0], c[1], c[2],
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[1]),
    HighOrderMesh_lookup(edge_midnodes_, c[1], c[2]),
    HighOrderMesh_lookup(edge_midnodes_, c[2], c[0])
  };
}

std::vector<int> HighOrderMesh::quad4_to_quad8(const int *c) const
{
  return {
    c[0], c[1], c[2], c[3],
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[1]),
    HighOrderMesh_lookup(edge_midnodes_, c[1], c[2]),
    HighOrderMesh_lookup(edge_midnodes_, c[2], c[3]),
    HighOrderMesh_lookup(edge_midnodes_, c[3], c[0])
  };
}

std::vector<int> HighOrderMesh::pyramid5_to_pyramid13(const int *c) const
{
  return {
    c[0], c[1], c[2], c[3], c[4],
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[1]),  // m01
    HighOrderMesh_lookup(edge_midnodes_, c[1], c[2]),  // m12
    HighOrderMesh_lookup(edge_midnodes_, c[2], c[3]),  // m23
    HighOrderMesh_lookup(edge_midnodes_, c[3], c[0]),  // m30
    HighOrderMesh_lookup(edge_midnodes_, c[0], c[4]),  // m04
    HighOrderMesh_lookup(edge_midnodes_, c[1], c[4]),  // m14
    HighOrderMesh_lookup(edge_midnodes_, c[2], c[4]),  // m24
    HighOrderMesh_lookup(edge_midnodes_, c[3], c[4])   // m34
  };
}

// ============================================================
// Write generated mid-nodes to Nastran BDF format
// ============================================================
void HighOrderMesh::write_midnodes(std::ofstream &fid, bool is_3d) const
{
  fid << "$\n$ Generated mid-edge nodes\n$\n";
  for (auto &pair : midnode_coords_) {
    int nid = pair.first;
    auto &c = pair.second;
    char line1[128], line2[128];
    if (is_3d) {
      std::snprintf(line1, sizeof(line1), "GRID*   %16d%16d%16.5f%16.5f",
                    nid, 0, c[0], c[1]);
      std::snprintf(line2, sizeof(line2), "*       %16.5f", c[2]);
    } else {
      std::snprintf(line1, sizeof(line1), "GRID*   %16d%16d%16.5f%16.5f",
                    nid, 0, c[0], c[1]);
      std::snprintf(line2, sizeof(line2), "*       %16d", 0);
    }
    fid << line1 << "\n" << line2 << "\n";
  }
}

void HighOrderMesh::write_midnodes_gmsh(std::ofstream &fid) const
{
  for (auto &pair : midnode_coords_) {
    fid << pair.first << " "
        << pair.second[0] << " "
        << pair.second[1] << " "
        << pair.second[2] << "\n";
  }
}
