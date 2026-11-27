#pragma once
// The only place that touches CurvedElements coefficient storage.
//
// GeometricRefit and the periodic coefficient copy live here, outside Netgen;
// they need read/write access to the per-edge and per-face coefficients, which
// upstream CurvedElements keeps private.  The Radia Netgen fork (r17) provides
// that access through the members used below.  When upstream grows an
// equivalent accessor, change only this header.
#include "../netgen_meshing.hpp"
#include <cstddef>

namespace radia_netgen_curving {

struct CurvedAccess
{
  netgen::CurvedElements & ce;

  int order () const { return ce.GetOrder(); }
  bool high_order () const { return ce.IsHighOrder(); }
  bool rational () const { return ce.PeriodicRational(); }

  // Throws unless the edge/face storage matches the given topology sizes and
  // every entity span agrees with its order.  Call before indexing entities.
  void validate (std::size_t nedges, std::size_t nfaces) const
  { ce.ValidatePeriodicStorage(nedges, nfaces); }

  // Indices are 0-based topology edge/face numbers.  Pointers are valid only
  // when the matching count is positive.
  int edge_count (int e) const { return ce.NumEdgeCoefficients(e); }
  int face_count (int f) const { return ce.NumFaceCoefficients(f); }
  netgen::Vec<3> * edge_coeffs (int e) { return ce.EdgeCoefficients(e); }
  netgen::Vec<3> * face_coeffs (int f) { return ce.FaceCoefficients(f); }
  int edge_order (int e) const { return ce.PeriodicEdgeOrder(e); }
  int face_order (int f) const { return ce.PeriodicFaceOrder(f); }
};

inline CurvedAccess curved_access (netgen::Mesh & mesh)
{ return CurvedAccess{mesh.GetCurvedElements()}; }

} // namespace radia_netgen_curving
