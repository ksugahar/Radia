#ifndef CUBIT_PERIODIC_CURVED_GEOMETRY_HPP
#define CUBIT_PERIODIC_CURVED_GEOMETRY_HPP

#include "NamedPeriodicMeshAdapter.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace cubit_mesh_export {

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
#ifdef COMPACT_NETGEN
  if (curved.PeriodicRational())
    throw std::runtime_error("Rational periodic geometry is unsupported");
#endif
  if (!curved.IsHighOrder()) return {};
#ifndef COMPACT_NETGEN
  throw std::runtime_error("Curved periodic synchronization requires the compact exporter coefficient API");
#else
  const auto &top = mesh.GetTopology();
  curved.ValidatePeriodicStorage(top.GetNEdges(), top.GetNFaces());
  using Edge = std::array<int, 2>;
  using Face = std::array<int, 3>;
  std::map<Edge, int> edges;
  std::map<Face, int> faces;
  for (int e = 0; e < int(top.GetNEdges()); ++e) {
    const auto vertices = top.GetEdgeVertices(e);
    Edge key{periodic_point_number(vertices[0]), periodic_point_number(vertices[1])};
    if (key[0] >= key[1] || !edges.emplace(key, e).second)
      throw std::runtime_error("Noncanonical or duplicate topology edge");
  }
  for (int f = 0; f < int(top.GetNFaces()); ++f) {
    const auto vertices = top.GetFaceVertices(f);
    if (vertices.Size() != 3) continue;
    Face key;
    for (int j = 0; j < 3; ++j) key[j] = periodic_point_number(vertices[j]);
    std::sort(key.begin(), key.end());
    if (!faces.emplace(key, f).second)
      throw std::runtime_error("Duplicate topology face");
  }
  struct Owner { int master, sign; };
  std::map<int, Owner> edge_owners, face_owners;
  auto own = [](auto &owners, int slave, int master, int sign) {
    auto inserted = owners.emplace(slave, Owner{master, sign});
    if (!inserted.second && (inserted.first->second.master != master ||
                            inserted.first->second.sign != sign))
      throw std::runtime_error("Conflicting periodic seam coefficient ownership");
  };
  for (const auto &pair : triangles) {
    Face slave = pair.slave; std::sort(slave.begin(), slave.end());
    own(face_owners, periodic_curved_lookup(faces, slave, periodic_curved_entity("face", slave)), periodic_curved_lookup(faces, pair.master, periodic_curved_entity("face", pair.master)), 1);
    for (int a = 0; a < 3; ++a) for (int b = a + 1; b < 3; ++b) {
      Edge master{pair.master[a], pair.master[b]};
      Edge target{pair.slave[a], pair.slave[b]};
      const int sign = target[0] < target[1] ? 1 : -1;
      std::sort(target.begin(), target.end());
      own(edge_owners, periodic_curved_lookup(edges, target, periodic_curved_entity("edge", target)), periodic_curved_lookup(edges, master, periodic_curved_entity("edge", master)), sign);
    }
  }
  std::vector<PeriodicCurvedTransfer> plan;
  auto stage = [&](bool face, const auto &owners) {
    std::set<int> sources;
    for (const auto &row : owners) sources.insert(row.second.master);
    for (const auto &row : owners) {
      const int slave = row.first, master = row.second.master;
      if (sources.count(slave))
        throw std::runtime_error("Periodic coefficient source and destination overlap");
      const int count = face ? curved.NumFaceCoefficients(master) : curved.NumEdgeCoefficients(master);
      const int target_count = face ? curved.NumFaceCoefficients(slave) : curved.NumEdgeCoefficients(slave);
      const int source_order = face ? curved.PeriodicFaceOrder(master) : curved.PeriodicEdgeOrder(master);
      const int target_order = face ? curved.PeriodicFaceOrder(slave) : curved.PeriodicEdgeOrder(slave);
      if (count != target_count || source_order != target_order ||
          source_order < 1 || source_order > order ||
          count != (face ? (source_order-1)*(source_order-2)/2 : source_order-1))
        throw std::runtime_error("Periodic coefficient order or layout mismatch");
      if (!count) continue;
      auto *src = face ? curved.FaceCoefficients(master) : curved.EdgeCoefficients(master);
      auto *dst = face ? curved.FaceCoefficients(slave) : curved.EdgeCoefficients(slave);
      for (int mode = 0; mode < count; ++mode) {
        const int sign = !face && mode == 1 ? row.second.sign : 1;
        PeriodicCurvedTransfer transfer{face, master, slave, mode, sign, {}, {}};
        for (int axis = 0; axis < 3; ++axis) {
          transfer.original[axis] = dst[mode](axis);
          transfer.expected[axis] = sign * src[mode](axis);
          if (!std::isfinite(transfer.original[axis]) || !std::isfinite(transfer.expected[axis]))
            throw std::runtime_error("Nonfinite periodic coefficient");
        }
        plan.push_back(transfer);
      }
    }
  };
  stage(false, edge_owners); stage(true, face_owners);
  auto coefficient = [&](const PeriodicCurvedTransfer &row) {
    return (row.face ? curved.FaceCoefficients(row.slave) : curved.EdgeCoefficients(row.slave)) + row.mode;
  };
  // Every operation that can reject topology/layout precedes the first write.
  try {
    for (const auto &row : plan)
      for (int axis = 0; axis < 3; ++axis) (*coefficient(row))(axis) = row.expected[axis];
    for (const auto &row : plan)
      for (int axis = 0; axis < 3; ++axis)
        if ((*coefficient(row))(axis) != row.expected[axis])
          throw std::runtime_error("Periodic paired coefficient equality failed");
  } catch (...) {
    for (const auto &row : plan)
      for (int axis = 0; axis < 3; ++axis) (*coefficient(row))(axis) = row.original[axis];
    throw;
  }
  return plan;
#endif
}

} // namespace cubit_mesh_export
#endif
