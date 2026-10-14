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
#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

// Face key for deduplication (triangular face — 3 sorted vertex IDs)
struct HoFaceKey3 {
  int n0, n1, n2;  // sorted ascending
  HoFaceKey3(int a, int b, int c) {
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    n0 = a; n1 = b; n2 = c;
  }
  bool operator==(const HoFaceKey3 &o) const {
    return n0 == o.n0 && n1 == o.n1 && n2 == o.n2;
  }
};
struct HoFaceKey3Hash {
  size_t operator()(const HoFaceKey3 &k) const {
    size_t h = std::hash<int>()(k.n0);
    h ^= std::hash<int>()(k.n1) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= std::hash<int>()(k.n2) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

// Volume key for deduplication (tet — 4 sorted vertex IDs)
struct HoElemKey4 {
  int n0, n1, n2, n3;  // sorted ascending
  HoElemKey4(int a, int b, int c, int d) {
    int v[4] = {a, b, c, d};
    std::sort(v, v + 4);
    n0 = v[0]; n1 = v[1]; n2 = v[2]; n3 = v[3];
  }
  bool operator==(const HoElemKey4 &o) const {
    return n0 == o.n0 && n1 == o.n1 && n2 == o.n2 && n3 == o.n3;
  }
};
struct HoElemKey4Hash {
  size_t operator()(const HoElemKey4 &k) const {
    size_t h = std::hash<int>()(k.n0);
    h ^= std::hash<int>()(k.n1) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= std::hash<int>()(k.n2) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= std::hash<int>()(k.n3) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
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

  // Cubit node ID -> Netgen PointIndex (1-based).  Needed by the
  // exporter to turn Cubit nodesets (GND, etc.) into Netgen CD3
  // (BBBND vertex) names that round-trip through .vol load.
  const std::unordered_map<int, int>& get_cubit_nid_to_ng_pi() const {
    return cubit_nid_to_ng_pi_;
  }

  // Get edge high-order nodes between two vertex node IDs.
  // Returns (order-1) node IDs along the edge from n0 to n1.
  std::vector<int> get_edge_nodes(int n0, int n1) const;

  // Get face interior nodes for a triangular face (order >= 3).
  // Vertices n0,n1,n2 in element-local order. Returns nodes in
  // GMSH convention relative to the given vertex ordering.
  std::vector<int> get_face_nodes(int n0, int n1, int n2) const;

  // Get volume interior nodes for a tet (order >= 4).
  // Vertices n0..n3 in element-local order.
  std::vector<int> get_vol_nodes(int n0, int n1, int n2, int n3) const;

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

  // Face -> ordered list of face interior node IDs (order >= 3)
  // Key: sorted 3-vertex triple. Value: nodes generated with the
  // canonical vertex ordering (v_sorted[0] -> v_sorted[1] -> v_sorted[2]).
  // Stored vertices for orientation recovery:
  struct FaceData {
    std::vector<int> node_ids;   // interior node IDs in generation order
    int gen_v0, gen_v1, gen_v2;  // vertex IDs in the order used during generation
  };
  std::unordered_map<HoFaceKey3, FaceData, HoFaceKey3Hash> face_ho_nodes_;

  // Volume interior nodes (order >= 4, TET only)
  // Key: sorted 4-vertex tuple. Value: node IDs + generation vertex order.
  struct VolData {
    std::vector<int> node_ids;
    int gen_v0, gen_v1, gen_v2, gen_v3;
  };
  std::unordered_map<HoElemKey4, VolData, HoElemKey4Hash> vol_ho_nodes_;

  int total_nodes_ = 0;
  int next_node_id_ = 0;
  int order_ = 2;

  // ---- Diagnostic counters for ACIS projection (Phase 2 callback) ----
  // Populated during curve_and_extract -> BuildCurvedElements -> project().
  // "reject" = projection displacement exceeded threshold, fallback applied.
  mutable long diag_project_calls_ = 0;
  mutable long diag_project_rejects_ = 0;
  mutable long diag_refine_success_ = 0;  // trim-then-refine recoveries
  mutable long diag_refine_tries_ = 0;
  mutable double diag_refine_max_tdisp_ = 0.0;
  mutable double diag_project_max_disp_ = 0.0;
  mutable int diag_project_max_disp_surf_ = -1;
  mutable std::unordered_map<int, long> diag_project_reject_per_surf_;
  // Per-path call counts: index 0..3 corresponds to the `path` value in
  // project() — 0=uv_direct_eval, 1=uv_guess_fallback, 2=nearest-vertex hint,
  // 3=no_hint_table_trimmed.
  mutable long diag_path_calls_[4] = {0, 0, 0, 0};
  mutable double diag_path_max_disp_[4] = {0, 0, 0, 0};
  // Path=0 displacement histogram (UV-hinted surface projections).
  // Bins: [0, 0.01), [0.01, 0.05), [0.05, 0.1), [0.1, 0.2), [0.2, 0.5),
  //       [0.5, 1.0), [1.0, inf).  Tells us the true distribution of
  //       curvature offsets Netgen is generating.
  mutable long diag_path0_hist_[7] = {0, 0, 0, 0, 0, 0, 0};
  mutable long diag_edge_calls_ = 0;
  mutable long diag_edge_rejects_ = 0;
  mutable double diag_edge_max_disp_ = 0.0;
  // Reject thresholds (mm). If projection displaces the query point by more
  // than this, fall back to identity (linear interpolation will pick up the
  // edge in the final pass).
  //
  // Surface threshold is tight because wrong-side projections on thin wire
  // cross-sections (full diameter ~6.3 mm) must be caught; raising to 5 mm
  // blows volume up +13%.
  //
  // Edge/curve threshold is *also* tight on 3turnCoil: accepting chord
  // midpoints that project more than ~1.5 mm away routinely picks up
  // cross-curve jumps (6.3 mm = full wire cross-section). With edge
  // threshold at 10 mm: volume +25.6%. At 1.5 mm: -0.80%.
  double project_reject_threshold_ = 1.5;
  double edge_reject_threshold_ = 1.5;

  // Accumulated reject detail rows (written to file at end of build())
  // Each string: "surfnr,cubit_sid,path,in_x,in_y,in_z,out_x,out_y,out_z,u_hint,v_hint,out_u,out_v,disp"
  mutable std::vector<std::string> project_reject_log_entries_;
  int project_reject_dump_max_ = 2000;

  // Phase1e diagnostic: per-surface vertex UV spread (max|u1-u0|, max|v1-v0| within tris)
  // Key = cubit_sid
  mutable std::unordered_map<int, double> diag_tri_max_du_;
  mutable std::unordered_map<int, double> diag_tri_max_dv_;
  mutable std::unordered_map<int, int> diag_tri_count_;

  // Per-surface "UV hint table": one entry per surface vertex, storing
  // the 3D position alongside its canonical UV. The project() callback uses
  // this to derive a UV guess on its no-hint first call (u_hint==v_hint==0)
  // via the nearest stored vertex — letting `closest_point_uv_guess` start a
  // local Newton iteration instead of the global `closest_point_trimmed`
  // that can jump to the wrong surface branch on thin wire cross-sections.
  // Key = cubit_sid. Each row = (x, y, z, u, v).
  std::unordered_map<int, std::vector<std::array<double,5>>> surf_vertex_uvs_;

  // Polar (singular) surfaces: those with at least one singular pole on the
  // underlying parametric surface (RefFace::num_poles() > 0) AND a planar
  // trim — the classic wire-terminal end-cap geometry where the parameter
  // space has a radial center that maps to a single 3D point ("center UV
  // arbitrary"). On these faces `closest_point_uv_guess` degenerates: the
  // Newton step cannot distinguish the true UV from a boundary-clamped one,
  // and displacements of wire-diameter magnitude appear (observed: 61 mm on
  // 3turnCoil's disk surf 1). Planar disks need no high-order curving in
  // the first place, so we short-circuit the callback to identity on them.
  // Populated in Phase1 from rf_cache.
  std::unordered_set<int> polar_surfaces_;
  mutable long diag_polar_skips_ = 0;

  // Per-surface axis-aligned bounding box derived from that surface's mesh
  // vertices plus a pad of ~1 wire-thickness.  Used by project() to reject
  // queries for inputs obviously outside the surface — a common failure
  // mode after `unite volume all` where Netgen's curving pipeline
  // occasionally asks us to project an interior midpoint to a distant
  // terminal face, producing tens-of-mm "projection" spikes.
  struct BBox { double xmin, ymin, zmin, xmax, ymax, zmax; };
  std::unordered_map<int, BBox> surf_bbox_;
  double surf_bbox_pad_ = 0.0;  // computed once from geometry diameter
};

#endif
