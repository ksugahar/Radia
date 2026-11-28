#ifndef CUBIT_TO_VOL_EXPORT_LIVE_VOL_MESH_HPP
#define CUBIT_TO_VOL_EXPORT_LIVE_VOL_MESH_HPP

// One curving pipeline for `export vol` and the in-memory mesh: Netgen's
// BuildCurvedElements on the Cubit callback geometry, the geometric refit for
// order >= 3 and, for a named periodic group, the master-to-slave coefficient
// copy.  The in-memory mesh keeps its LiveCubitGeometry, so it can be curved
// again at another order while the Cubit session that built it is alive.

#include "callbackgeom.hpp"
#include "PeriodicCurvedGeometry.hpp"
#include "netgen_curving/geometric_refit.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace cubit_to_vol_export {

class LiveCubitGeometry : public netgen::CallbackGeometry
{
public:
  using netgen::CallbackGeometry::CallbackGeometry;

  // Owner of the state the callbacks read (the NetgenCurver).
  std::shared_ptr<void> owner;
  // Re-applied after every curving of a mesh with a named periodic group.
  std::vector<PeriodicTrianglePair> periodic_triangles;
  // The Cubit model the callbacks query; checked before every re-curving.
  std::function<std::string()> model_fingerprint;
  std::string built_fingerprint;
};

// BuildCurvedElements(order) and, for order >= 3, the geometric refit.
// Throws when order >= 2 produces no high-order data.
radia_netgen_curving::RefitStats curve_mesh(netgen::Mesh &mesh, int order);

// Curve an in-memory mesh again at `order` (1..5).  Throws when the mesh has
// no LiveCubitGeometry or the Cubit model changed since the mesh was built.
void recurve(netgen::Mesh &mesh, int order);

}  // namespace cubit_to_vol_export

#endif
