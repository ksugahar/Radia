#ifndef CUBIT_NAMED_PERIODIC_MESH_ADAPTER_HPP
#define CUBIT_NAMED_PERIODIC_MESH_ADAPTER_HPP

#include <meshing.hpp>
#include <stdexcept>
#include <string>

namespace cubit_mesh_export {

// These are the two declared build variants, not a runtime fallback. External
// point and descriptor numbers always use the VOL format's one-based scheme.
inline netgen::PointIndex periodic_point_index(int number)
{
  if (number < 1) throw std::runtime_error("Invalid one-based point number");
#ifdef COMPACT_NETGEN
  return netgen::PointIndex(number);
#else
  return netgen::PointIndex::FromNr1(number);
#endif
}

inline int periodic_point_number(netgen::PointIndex index)
{
#ifdef COMPACT_NETGEN
  return static_cast<int>(index);
#else
  return index.Nr1();
#endif
}

#ifdef COMPACT_NETGEN
inline int periodic_face_index(int number)
{
  if (number < 1) throw std::runtime_error("Invalid one-based face descriptor");
  return number;
}
inline int periodic_face_number(int index) { return index; }
#else
inline netgen::FaceRegionIndex periodic_face_index(int number)
{
  if (number < 1) throw std::runtime_error("Invalid one-based face descriptor");
  return netgen::FaceRegionIndex::FromNr1(number);
}
inline int periodic_face_number(netgen::FaceRegionIndex index) { return index.Nr1(); }
#endif

// Resolve a boundary-condition property, never reinterpret it as a face index.
// Compact Netgen owns a separate BC-name table; full Netgen stores names on
// face descriptors. Multiple descriptors for one BC must agree on the name.
inline std::string periodic_boundary_label(const netgen::Mesh &mesh, int bc)
{
  if (bc < 1) throw std::runtime_error("Invalid boundary-condition property");
#ifdef COMPACT_NETGEN
  return mesh.GetBCName(bc - 1);
#else
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
#endif
}

} // namespace cubit_mesh_export
#endif
