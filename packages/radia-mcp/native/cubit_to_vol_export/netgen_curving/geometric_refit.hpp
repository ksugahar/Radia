#pragma once
// Geometric refit of curved-element coefficients of order >= 3.
// Moved out of the Radia Netgen fork (libsrc/meshing/geometric_refit.*).
namespace netgen { class Mesh; }

namespace radia_netgen_curving {

struct RefitStats
{
  int edges_tried = 0, edges_accepted = 0, edges_failed = 0, edges_rejected_distortion = 0;
  int faces_tried = 0, faces_accepted = 0, faces_failed = 0, faces_rejected_distortion = 0;
  double edge_dist_before = 0, edge_dist_after = 0;   // max sampled distance
  double face_dist_before = 0, face_dist_after = 0;
};

// Refit the coefficients to minimise the distance to the mesh geometry; call
// after Curve(order).  Orders <= 2 are left unchanged.  A failing geometry
// callback leaves every entity at either its original or its accepted value.
RefitStats geometric_refit (netgen::Mesh & mesh);

} // namespace radia_netgen_curving
