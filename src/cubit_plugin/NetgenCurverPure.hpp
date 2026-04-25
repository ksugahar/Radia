#ifndef NETGEN_CURVER_PURE_HPP
#define NETGEN_CURVER_PURE_HPP

/* ===================================================================
 * NetgenCurverPure — Cubit-free high-order mesh curving via Netgen
 *
 * Accepts mesh data and projection/normal callbacks from Python.
 * Does NOT include any Cubit headers or link to Cubit libraries.
 *
 * Used by the pybind11 module (radia_cubit_mesh.pyd) to enable
 * high-order curved meshes from external Python 3.12 via
 * cubit Python API + Netgen C++.
 * =================================================================== */

#include <memory>
#include <vector>
#include <array>
#include <unordered_map>
#include <functional>

namespace netgen { class Mesh; }

// Callback types matching CallbackGeometry (must match callbackgeom.hpp)
// has_hint distinguishes "caller supplied a UV hint" from "no hint, use
// global projection" — added 2026 to fix UV=(0,0) ambiguity at parameter
// origin. See callbackgeom.hpp comment near line 30 for rationale.
using ProjectFunc = std::function<std::tuple<double,double,double,double,double>
                                  (int surfnr, double x, double y, double z,
                                   double u_hint, double v_hint, bool has_hint)>;
using NormalFunc = std::function<std::tuple<double,double,double>
                                (int surfnr, double x, double y, double z)>;
using EdgeProjectFunc = std::function<std::tuple<double,double,double>
                                      (int surfnr1, int surfnr2,
                                       double x, double y, double z)>;

// Per-surface data passed from Python
struct SurfaceInfo {
  int id;                                    // Cubit surface ID
  std::vector<std::array<int,3>> tris;       // triangle connectivity (node IDs)
  std::vector<std::array<int,4>> quads;      // quad connectivity (node IDs)
  std::vector<std::pair<int, std::array<double,2>>> uvs;  // (node_id, {u, v})
};

// Volume element types
enum PureElemType { PURE_TET = 0, PURE_HEX = 1, PURE_PRISM = 2, PURE_PYRAMID = 3 };

struct VolumeElement {
  PureElemType type;
  std::vector<int> conn;   // node IDs
  int domain = 1;          // domain index (1-based, maps to material)
};

// Edge (segment) elements on geometry curves
struct EdgeInfo {
  int surfnr1, surfnr2;                    // 1-based FaceDescriptor indices
  std::vector<std::array<int,2>> segments;  // node ID pairs
  std::vector<std::array<double,2>> dists;  // normalized curve parameter [0,1] per segment endpoint
};

class NetgenCurverPure
{
public:
  std::shared_ptr<netgen::Mesh> build(
      const std::vector<double>& node_coords,
      const std::vector<int>&    node_ids,
      const std::vector<VolumeElement>& vol_elems,
      const std::vector<SurfaceInfo>& surfaces,
      ProjectFunc project_func,
      NormalFunc  normal_func,
      int order,
      EdgeProjectFunc edge_project_func = nullptr,
      const std::vector<EdgeInfo>& edge_elements = {});

  // Access the built mesh
  std::shared_ptr<netgen::Mesh> get_mesh() const { return ng_mesh_; }

private:
  bool build_netgen_mesh(
      const std::vector<double>& node_coords,
      const std::vector<int>&    node_ids,
      const std::vector<VolumeElement>& vol_elems,
      const std::vector<SurfaceInfo>& surfaces,
      const std::vector<EdgeInfo>& edge_elements);

  bool attach_callback_geometry(ProjectFunc proj, NormalFunc norm, int num_surfaces,
                                EdgeProjectFunc edge_proj = nullptr);
  bool curve(int order);

  std::shared_ptr<netgen::Mesh> ng_mesh_;
  std::unordered_map<int, int> cubit_nid_to_ng_pi_;  // cubit node ID -> netgen PointIndex
};

#endif // NETGEN_CURVER_PURE_HPP
