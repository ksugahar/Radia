#ifndef CUBIT_PERIODIC_CURVED_GEOMETRY_HPP
#define CUBIT_PERIODIC_CURVED_GEOMETRY_HPP

#include "NamedPeriodicMeshAdapter.hpp"
#include "netgen_curving/curved_access.hpp"
#include "netgen_curving/periodic_copy.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace cubit_to_vol_export {

// Missing topology reports its identity instead of an opaque map::at failure.
template <class Key, class Value>
inline const Value &periodic_curved_lookup(const std::map<Key, Value> &entries,
    const Key &key, const std::string &entity)
{
  const auto found = entries.find(key);
  if (found == entries.end()) throw std::runtime_error("Missing periodic " + entity);
  return found->second;
}

template <std::size_t N>
inline std::string periodic_curved_entity(const char *kind, const std::array<int, N> &key)
{
  std::string text = std::string(kind) + " [";
  for (std::size_t i = 0; i < N; ++i)
    text += (i ? "," : "") + std::to_string(key[i]);
  return text + "]";
}

struct PeriodicTrianglePair {
  std::array<int, 3> master, slave; // slave entries follow the master mapping
};

// Connectivity is checked per BC in both directions, before any coefficient
// or identification write. A node-set match alone does not prove face coverage.
inline std::vector<PeriodicTrianglePair> periodic_triangle_plan(
    const netgen::Mesh &mesh, const std::map<int, int> &bc_pairs,
    const std::map<int, int> &points)
{
  using Triangle = std::array<int, 3>;
  std::map<int, std::set<Triangle>> patches;
  std::set<int> selected;
  for (const auto &p : bc_pairs) { selected.insert(p.first); selected.insert(p.second); }
  for (const auto &surface : mesh.SurfaceElements()) {
    const int bc = mesh.GetFaceDescriptor(surface.GetIndex()).BCProperty();
    if (!selected.count(bc)) continue;
    if (surface.GetType() != netgen::TRIG || surface.GetNP() != 3)
      throw std::runtime_error("Named periodic geometry requires triangular linear connectivity");
    Triangle tri;
    for (int j = 0; j < 3; ++j) tri[j] = periodic_point_number(surface[j]);
    std::sort(tri.begin(), tri.end());
    if (tri[0] == tri[1] || tri[1] == tri[2] || !patches[bc].insert(tri).second)
      throw std::runtime_error("Duplicate or degenerate periodic triangle");
  }
  std::vector<PeriodicTrianglePair> result;
  for (const auto &patch : bc_pairs) {
    std::set<Triangle> covered;
    for (const auto &master : periodic_curved_lookup(patches, patch.first, "BC patch " + std::to_string(patch.first))) {
      Triangle mapped;
      for (int j = 0; j < 3; ++j) mapped[j] = periodic_curved_lookup(points, master[j], "point " + std::to_string(master[j]));
      auto key = mapped; std::sort(key.begin(), key.end());
      if (!periodic_curved_lookup(patches, patch.second, "BC patch " + std::to_string(patch.second)).count(key) || !covered.insert(key).second)
        throw std::runtime_error("Periodic triangle correspondence is incomplete or ambiguous");
      result.push_back({master, mapped});
    }
    if (covered != periodic_curved_lookup(patches, patch.second, "BC patch " + std::to_string(patch.second)))
      throw std::runtime_error("Incomplete reverse periodic triangle coverage");
  }
  return result;
}

struct PeriodicCurvedTransfer {
  bool face;
  int master, slave, mode, sign;
  std::array<double, 3> original, expected;
};

// q2 has one even edge mode; q3 adds an odd edge mode and a symmetric
// triangular face bubble. Translation belongs to corners, not displacements.
inline std::vector<PeriodicCurvedTransfer> synchronize_periodic_curvature(
    netgen::Mesh &mesh, const std::vector<PeriodicTrianglePair> &triangles)
{
  auto &curved = mesh.GetCurvedElements();
  const int order = curved.GetOrder();
  if (order < 1 || order > 3)
    throw std::runtime_error("Named periodic geometry supports polynomial orders 1 through 3 only");
  for (const auto &element : mesh.VolumeElements())
    if (element.GetType() != netgen::TET || element.GetNP() != 4)
      throw std::runtime_error("Periodic export requires linear tetrahedral connectivity");
  if (radia_netgen_curving::CurvedAccess{curved}.rational())
    throw std::runtime_error("Rational periodic geometry is unsupported");
  if (!curved.IsHighOrder()) return {};
  std::vector<std::array<int, 3>> master, slave;
  for (const auto &pair : triangles) { master.push_back(pair.master); slave.push_back(pair.slave); }
  std::vector<PeriodicCurvedTransfer> plan;
  for (const auto &row : radia_netgen_curving::copy_periodic_coefficients(mesh, master, slave))
    plan.push_back({row.face, row.master, row.slave, row.mode, row.sign, row.original, row.expected});
  return plan;
}

} // namespace cubit_to_vol_export
#endif
