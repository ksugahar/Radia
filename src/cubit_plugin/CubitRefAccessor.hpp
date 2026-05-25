#pragma once
//=============================================================================
// CubitRefAccessor: thin wrappers around the Cubit 2025.12 RefEntity factory.
//
// Background (2026-05-25):
//   Cubit 2025.12 (CUBIT_VERSION_ALLINT 17040+) removed the convenience
//   static `GeometryQueryTool::instance()` accessor that older Radia
//   plugin code relied on for `get_ref_face(int)` / `get_ref_edge(int)`
//   / `get_ref_volume(int)` lookups by ID.  The replacement path is
//
//     CubitCoreModel::instance()        // process singleton (still here)
//       ->cgmapp()                       // CGMApp* owned by the model
//       ->factory()                      // RefEntityFactory* per CGMApp
//       ->get_ref_face(id) / _edge / _volume
//
//   That chain is verbose and repeated; collect it here so the call sites
//   stay one-liners and the migration to any future Cubit-N+1 accessor
//   pattern is a single-file edit.
//
// Usage:
//   #include "CubitRefAccessor.hpp"
//   RefFace*   rf = radia::cubit_get_ref_face(sid);
//   RefEdge*   re = radia::cubit_get_ref_edge(cid);
//   RefVolume* rv = radia::cubit_get_ref_volume(vid);
//
// All accessors return nullptr if the ID does not resolve (same contract as
// the old GeometryQueryTool methods).  Callers MUST null-check.
//=============================================================================

#include "CubitCoreModel.hpp"
#include "CGMApp.hpp"
#include "RefEntityFactory.hpp"

class RefFace;
class RefEdge;
class RefVolume;

namespace radia {

inline RefEntityFactory* cubit_factory() {
  CubitCoreModel* model = CubitCoreModel::instance();
  if (!model) return nullptr;
  CGMApp* app = model->cgmapp();
  if (!app) return nullptr;
  return app->factory();
}

inline RefFace* cubit_get_ref_face(int id) {
  RefEntityFactory* f = cubit_factory();
  return f ? f->get_ref_face(id) : nullptr;
}

inline RefEdge* cubit_get_ref_edge(int id) {
  RefEntityFactory* f = cubit_factory();
  return f ? f->get_ref_edge(id) : nullptr;
}

inline RefVolume* cubit_get_ref_volume(int id) {
  RefEntityFactory* f = cubit_factory();
  return f ? f->get_ref_volume(id) : nullptr;
}

}  // namespace radia
