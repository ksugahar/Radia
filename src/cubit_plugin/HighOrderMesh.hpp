#ifndef HIGH_ORDER_MESH_HPP
#define HIGH_ORDER_MESH_HPP

#include <vector>
#include <array>
#include <iosfwd>
#include <unordered_map>
#include <cstdint>

class MeshExportInterface;

// Geometric entity classification for mesh nodes
enum class GeomOwner { VERTEX, CURVE, SURFACE, VOLUME };

struct NodeGeomInfo {
  GeomOwner type;
  int entity_id;
};

// Edge key: ordered pair of node IDs (for deduplication)
struct EdgeKey {
  int n0, n1;
  EdgeKey(int a, int b) : n0(a < b ? a : b), n1(a < b ? b : a) {}
  bool operator==(const EdgeKey &o) const { return n0 == o.n0 && n1 == o.n1; }
};

struct EdgeKeyHash {
  size_t operator()(const EdgeKey &k) const {
    return std::hash<int64_t>()(((int64_t)k.n0 << 32) | (int64_t)(unsigned)k.n1);
  }
};

// ============================================================
// HighOrderMesh
//
// Converts a 1st-order Cubit mesh to 2nd-order by inserting
// mid-edge nodes and projecting them onto ACIS geometry via
// CubitInterface API.
//
// All elements get uniform order (no mixed-order support).
// ============================================================
class HighOrderMesh
{
public:
  HighOrderMesh();

  // Build high-order mesh data. Call after iface->initialize_export().
  bool build(MeshExportInterface *iface);

  // Total node count (original + generated mid-nodes)
  int get_num_nodes() const { return total_nodes_; }

  // Get coordinates for any node ID (original or generated)
  std::array<double,3> get_node_coords(int node_id) const;

  // Element conversion: 1st-order → 2nd-order
  // Nastran node ordering for each element type.
  // Input: vertex node IDs from MeshExportInterface connectivity.
  std::vector<int> tet4_to_tet10(const int *conn) const;
  std::vector<int> hex8_to_hex20(const int *conn) const;
  std::vector<int> wedge6_to_wedge15(const int *conn) const;
  std::vector<int> tri3_to_tri6(const int *conn) const;
  std::vector<int> quad4_to_quad8(const int *conn) const;
  std::vector<int> pyramid5_to_pyramid13(const int *conn) const;

  // Write generated mid-nodes to a Nastran-format stream
  void write_midnodes(std::ofstream &fid, bool is_3d) const;

  // Write generated mid-nodes to Gmsh format (id x y z)
  void write_midnodes_gmsh(std::ofstream &fid) const;

private:
  void classify_nodes();
  int get_or_create_midnode(int n0, int n1);

  // Node classification
  std::unordered_map<int, NodeGeomInfo> node_owner_;
  // Node → list of curves it touches (for vertex nodes, includes parent curves)
  std::unordered_map<int, std::vector<int>> node_curves_;
  // Node → list of surfaces it touches
  std::unordered_map<int, std::vector<int>> node_surfaces_;

  // Edge mid-node table (mutable because get_or_create modifies it)
  mutable std::unordered_map<EdgeKey, int, EdgeKeyHash> edge_midnodes_;

  // Generated mid-node coordinates
  mutable std::unordered_map<int, std::array<double,3>> midnode_coords_;

  int total_nodes_ = 0;
  mutable int next_node_id_ = 0;
};

#endif
