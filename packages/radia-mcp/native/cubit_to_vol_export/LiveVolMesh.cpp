#include "LiveVolMesh.hpp"

#include <stdexcept>

namespace cubit_to_vol_export {

radia_netgen_curving::RefitStats curve_mesh(netgen::Mesh &mesh, int order)
{
  if (order < 1 || order > 5)
    throw std::invalid_argument("curve order must be 1-5, got " + std::to_string(order));
  mesh.BuildCurvedElements(order);
  if (order == 1) return {};
  if (!mesh.GetCurvedElements().IsHighOrder())
    throw std::runtime_error("BuildCurvedElements(" + std::to_string(order) +
                             ") produced no high-order data");
  return radia_netgen_curving::geometric_refit(mesh);
}

void recurve(netgen::Mesh &mesh, int order)
{
  auto geometry = std::dynamic_pointer_cast<LiveCubitGeometry>(mesh.GetGeometry());
  if (!geometry)
    throw std::runtime_error("this mesh carries no live Cubit geometry; only a mesh from "
                             "mesh_from_cubit can be curved again");
  if (!geometry->model_fingerprint)
    throw std::runtime_error("the live Cubit geometry has no model fingerprint");
  // A bad order is refused before the mesh is touched, so it does not mark the mesh failed.
  if (order < 1 || order > 5)
    throw std::invalid_argument("curve order must be 1-5, got " + std::to_string(order));
  if (!geometry->failure.empty())
    throw std::runtime_error("an earlier curving of this mesh failed (" + geometry->failure +
                             "); build the mesh again");
  const std::string now = geometry->model_fingerprint();
  if (now != geometry->built_fingerprint)
    throw std::runtime_error("the Cubit model changed since this mesh was built (" +
                             geometry->built_fingerprint + " -> " + now +
                             "); build the mesh again");
  try {
    curve_mesh(mesh, order);
    if (!geometry->periodic_triangles.empty())
      synchronize_periodic_curvature(mesh, geometry->periodic_triangles);
  } catch (const std::exception &error) {
    geometry->failure = error.what();
    throw std::runtime_error(std::string("curving at order ") + std::to_string(order) +
                             " failed and left this mesh unusable (" + error.what() +
                             "); build the mesh again");
  }
}

}  // namespace cubit_to_vol_export
