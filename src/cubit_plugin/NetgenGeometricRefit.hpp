#ifndef NETGEN_GEOMETRIC_REFIT_HPP
#define NETGEN_GEOMETRIC_REFIT_HPP

// ============================================================
// Geometric refit of Netgen curved-element coefficients
//
// Netgen's BuildCurvedElements fits each edge/face to the displacement of
// the chord point onto the geometry, at fixed parameters, in L2.  On a
// symmetric edge the normal part of that displacement is even, so adding an
// odd polynomial degree cannot reduce the shape error: geometry error stalls
// from p=2 to p=3 and from p=4 to p=5.
//
// This pass keeps Netgen's coefficient layout and basis and replaces the
// coefficients by a least-squares minimiser of the distance between the
// curved edge/face and the geometry (Gauss-Newton, then a Levenberg-Marquardt
// polish).  The tangential freedom of the parametrisation is what lets an odd
// degree improve the shape; an edge's lateral slide inside its surface, which
// the distance cannot see, stays at Netgen's value.  A refit is kept only when
// the sampled distance to the geometry decreases and the parametrisation does
// not degenerate; otherwise Netgen's coefficients stay untouched.
// ============================================================

namespace netgen { class Mesh; }

struct GeometricRefitStats
{
  // failed: a projection was not trustworthy.  rejected_distortion: the
  // refit was closer to the geometry but would distort an element.
  int edges_tried = 0, edges_accepted = 0, edges_failed = 0, edges_rejected_distortion = 0;
  int faces_tried = 0, faces_accepted = 0, faces_failed = 0, faces_rejected_distortion = 0;
  double edge_dist_before = 0, edge_dist_after = 0;   // max sampled distance
  double face_dist_before = 0, face_dist_after = 0;
};

// Refit the curved-element coefficients of `mesh` in place.  The mesh must
// carry its projection geometry and already have BuildCurvedElements applied.
GeometricRefitStats geometric_refit(netgen::Mesh & mesh);

#endif
