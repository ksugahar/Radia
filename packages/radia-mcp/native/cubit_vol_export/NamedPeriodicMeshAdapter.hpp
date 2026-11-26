#ifndef CUBIT_NAMED_PERIODIC_MESH_ADAPTER_HPP
#define CUBIT_NAMED_PERIODIC_MESH_ADAPTER_HPP

#include "netgen_meshing.hpp"
#include <stdexcept>
#include <string>

namespace cubit_mesh_export {

// External point and descriptor numbers use the VOL format's one-based scheme.
inline netgen::PointIndex periodic_point_index(int number)
{
  if (number < 1) throw std::runtime_error("Invalid one-based point number");
  return netgen::PointIndex::FromNr1(number);
}

inline int periodic_point_number(netgen::PointIndex index) { return index.Nr1(); }

inline netgen::FaceRegionIndex periodic_face_index(int number)
{
  if (number < 1) throw std::runtime_error("Invalid one-based face descriptor");
  return netgen::FaceRegionIndex::FromNr1(number);
}
inline int periodic_face_number(netgen::FaceRegionIndex index) { return index.Nr1(); }

// Resolve a boundary-condition property, never reinterpret it as a face index.
// Netgen stores boundary names on the face descriptors (face regions); a BC
// property may own several descriptors, and they must agree on the name.
inline std::string periodic_boundary_label(const netgen::Mesh &mesh, int bc)
{
  if (bc < 1) throw std::runtime_error("Invalid boundary-condition property");
  bool found = false;
  std::string label;
  for (int fd = 1; fd <= mesh.GetNFD(); ++fd) {
    const auto &descriptor = mesh.GetFaceDescriptor(periodic_face_index(fd));
    if (descriptor.BCProperty() != bc) continue;
    const auto &name = descriptor.GetBCName();
    if (found && label != name)
      throw std::runtime_error("Conflicting descriptor labels for BC " + std::to_string(bc));
    label = name;
    found = true;
  }
  if (!found) throw std::runtime_error("Missing descriptor for BC " + std::to_string(bc));
  return label;
}

} // namespace cubit_mesh_export
#endif
