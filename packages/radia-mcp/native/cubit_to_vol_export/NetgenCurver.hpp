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
namespace cubit_to_vol_export { class CurveIdentity; }


class NetgenCurver
{
public:
  NetgenCurver();
  ~NetgenCurver();

  // Build curved mesh data for the given polynomial order (>= 2).
  // Reads linear mesh from MeshData; geometry queries via CubitInterface.
  // Returns false on failure.
  bool build(const MeshData &md, int order);

  // Shared production mapping entry; also exercised without loading Cubit CAD.
  static cubit_to_vol_export::CurveIdentity curve_identities(
      const netgen::Mesh & mesh, const std::unordered_map<int, int> & curve_descriptors);
  std::string edge_diagnostics_json() const;

  int get_order() const { return order_; }

  // Access the internal netgen::Mesh (for the .vol writer)
  std::shared_ptr<netgen::Mesh> get_ng_mesh() const { return ng_mesh_; }

  // Hand the mesh over to its LiveCubitGeometry: the geometry then owns this
  // curver, so the curver must not keep owning the mesh (no reference cycle).
  void release_mesh() { ng_mesh_.reset(); }

  // Cubit node ID -> Netgen PointIndex (1-based).  Needed by the
  // exporter to turn Cubit nodesets (GND, etc.) into Netgen CD3
  // (BBBND vertex) names that round-trip through .vol load.
  const std::unordered_map<int, int>& get_cubit_nid_to_ng_pi() const {
    return cubit_nid_to_ng_pi_;
  }

  // Only curves whose complete linear segmentation was exported are present.
  // Values are 1-based Netgen EdgeDescriptor indices, also used by CD2 names.
  const std::unordered_map<int, int>& get_cubit_curve_to_ng_edge() const {
    return cubit_curve_to_ng_edge_;
  }


private:
  std::unordered_map<int, int> cubit_curve_to_ng_edge_;
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

  int order_ = 2;

  // ---- Diagnostic counters for ACIS projection (Phase 2 callback) ----
  // Populated during curve_and_extract -> BuildCurvedElements -> project().
  // "reject" = projection displacement exceeded threshold, fallback applied.
  mutable long diag_project_calls_ = 0;
  mutable long diag_project_rejects_ = 0;
  mutable double diag_project_max_disp_ = 0.0;
  mutable int diag_project_max_disp_surf_ = -1;
  mutable long diag_edge_calls_ = 0;
  mutable long diag_edge_rejects_ = 0;
  mutable double diag_edge_max_disp_ = 0.0;
  mutable long diag_edge_failures_ = 0;
  mutable std::string diag_edge_worst_;
  std::string last_error_;
  // Reject thresholds in the CAD model's length units; no unit conversion.
  // Known-edge reprojection beyond edge_reject_threshold_ fails the export.
  // Arc interpolation with displacement greater than its chord also fails.
  // Neither case substitutes identity or linear interpolation (user decision
  // 2026-10-24). Surface projection rejection has its separate policy below.
  //
  // Historical 3turnCoil provenance (millimetre CAD model, former identity
  // fallback): a surface threshold of 5 mm produced +13% volume; an edge
  // threshold of 10 mm produced +25.6%, and 1.5 mm produced -0.80%.
  // Wire cross-section diameter was about 6.3 mm. These historical fallback
  // results do not certify the current fail-loud export path.
  double project_reject_threshold_ = 1.5;
  double edge_reject_threshold_ = 1.5;

  // Accumulated reject detail rows (written to file at end of build())
  // Each string: "surfnr,cubit_sid,path,in_x,in_y,in_z,out_x,out_y,out_z,u_hint,v_hint,out_u,out_v,disp"

  // Phase1e diagnostic: per-surface vertex UV spread (max|u1-u0|, max|v1-v0| within tris)
  // Key = cubit_sid

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
