#ifndef NETGEN_CURVER_HPP
#define NETGEN_CURVER_HPP

// ============================================================
// NetgenCurver
//
// High-order mesh curving via Netgen's BuildCurvedElements.
// Bridges Cubit's ACIS geometry (via cubit_geom) with Netgen's
// CallbackGeometry to produce properly curved elements at any
// polynomial order >= 2.
//
// For order 2: edge mid-nodes curved onto geometry
// For order 3+: edge + face + interior nodes via polynomial fit
// ============================================================

#include <vector>
#include <array>
#include <unordered_map>
#include <iosfwd>
#include <memory>
#include <cstdint>

class MeshData;

// Forward declare Netgen types (avoid including heavy headers here)
namespace netgen {
  class Mesh;
  class CallbackGeometry;
  class NetgenGeometry;
}

// Edge key for deduplication (same as HighOrderMesh)
struct HoEdgeKey {
  int n0, n1;
  HoEdgeKey(int a, int b) : n0(a < b ? a : b), n1(a < b ? b : a) {}
  bool operator==(const HoEdgeKey &o) const { return n0 == o.n0 && n1 == o.n1; }
};
struct HoEdgeKeyHash {
  size_t operator()(const HoEdgeKey &k) const {
    return std::hash<int64_t>()(((int64_t)k.n0 << 32) | (int64_t)(unsigned)k.n1);
  }
};

class NetgenCurver
{
public:
  NetgenCurver();
  ~NetgenCurver();

  // Build curved mesh data for the given polynomial order (>= 2).
  // Reads linear mesh from MeshData; geometry queries via CubitInterface.
  // Returns false on failure.
  bool build(const MeshData &md, int order);

  // Total node count (original + high-order nodes)
  int get_num_nodes() const { return total_nodes_; }

  int get_order() const { return order_; }

  // Access the internal netgen::Mesh (for pybind11 export)
  std::shared_ptr<netgen::Mesh> get_ng_mesh() const { return ng_mesh_; }

  // Get coordinates for any node ID (original or generated)
  std::array<double,3> get_node_coords(int node_id) const;

  // Get edge high-order nodes between two vertex node IDs.
  // Returns (order-1) node IDs along the edge from n0 to n1.
  std::vector<int> get_edge_nodes(int n0, int n1) const;

  // Write all high-order nodes (beyond vertex nodes) to Gmsh format
  void write_ho_nodes_gmsh(std::ofstream &fid) const;

  // Write all high-order nodes to Nastran format (order 2 only)
  void write_ho_nodes_nastran(std::ofstream &fid, bool is_3d) const;

  // Get the Gmsh element type code for a given base type and order
  static int gmsh_element_type(int base_nn, int elem_type, int order);

private:
  // Phase 1: Convert linear mesh (from MeshData) to netgen::Mesh
  bool build_netgen_mesh(const MeshData &md);

  // Phase 2: Build CallbackGeometry with ACIS projection lambdas
  bool attach_callback_geometry();

  // Phase 3: Curve and extract
  bool curve_and_extract(int order);

  // Surface mapping: Cubit surface ID <-> Netgen FaceDescriptor index
  std::unordered_map<int, int> cubit_sid_to_ng_fd_;
  std::unordered_map<int, int> ng_fd_to_cubit_sid_;

  // Node mapping: Cubit node ID <-> Netgen PointIndex
  std::unordered_map<int, int> cubit_nid_to_ng_pi_;
  std::unordered_map<int, int> ng_pi_to_cubit_nid_;

  // The netgen mesh
  std::shared_ptr<netgen::Mesh> ng_mesh_;

  // Generated high-order node coordinates (keyed by new node ID)
  std::unordered_map<int, std::array<double,3>> ho_node_coords_;

  // Edge -> ordered list of high-order node IDs (order-1 nodes per edge)
  // Key is ordered (min,max) pair; value is nodes from min→max direction
  std::unordered_map<HoEdgeKey, std::vector<int>, HoEdgeKeyHash> edge_ho_nodes_;

  int total_nodes_ = 0;
  int next_node_id_ = 0;
  int order_ = 2;
};

#endif
