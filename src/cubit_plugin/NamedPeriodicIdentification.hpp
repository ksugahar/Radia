#ifndef CUBIT_NAMED_PERIODIC_IDENTIFICATION_HPP
#define CUBIT_NAMED_PERIODIC_IDENTIFICATION_HPP

#include "NamedPeriodicMeshAdapter.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cubit_mesh_export {

// BC numbers are one-based boundary-condition properties, never FD indices.
struct NamedPeriodicContract {
  std::string name, master_label, slave_label, coordinate_unit;
  std::array<double, 3> translation{};
  double matching_distance = 0;
  std::vector<std::pair<int, int>> boundary_pairs;
};

struct NamedPeriodicResult {
  std::vector<std::pair<int, int>> point_pairs;
  std::map<int, int> descriptor_to_boundary;
  double maximum_distance = 0;
};

inline std::string periodic_geometry_serialization(const netgen::Mesh &mesh)
{
  std::ostringstream saved;
  mesh.Save(saved);
  std::istringstream input(saved.str());
  std::ostringstream geometry;
  std::string line;
  while (std::getline(input, line)) {
    if (line == "identifications" || line == "identificationtypes" ||
        line == "identificationnames") {
      std::string count_line;
      if (!std::getline(input, count_line))
        throw std::runtime_error("Truncated periodic metadata count");
      std::istringstream count_input(count_line);
      int count = -1;
      if (!(count_input >> count) || count < 0)
        throw std::runtime_error("Invalid periodic metadata count");
      const int rows = line == "identifications" ? count : 1;
      if (line != "identifications" && count != 1)
        throw std::runtime_error("Only one named periodic group is supported");
      for (int row = 0; row < rows; ++row)
        if (!std::getline(input, count_line))
          throw std::runtime_error("Truncated periodic metadata rows");
      continue;
    }
    geometry << line << '\n';
  }
  return geometry.str();
}

inline NamedPeriodicResult apply_named_periodic_identification(
    netgen::Mesh &mesh, const NamedPeriodicContract &contract)
{
  if (contract.name != "kt" || contract.coordinate_unit != "m" ||
      contract.master_label.empty() || contract.slave_label.empty() ||
      contract.master_label == contract.slave_label ||
      !std::isfinite(contract.matching_distance) ||
      contract.matching_distance <= 0 || contract.boundary_pairs.empty())
    throw std::runtime_error("Explicit named periodic contract is incomplete");
  for (double value : contract.translation)
    if (!std::isfinite(value))
      throw std::runtime_error("Nonfinite periodic translation");

  auto &ident = mesh.GetIdentifications();
  if (ident.GetMaxNr() != 0 || !ident.GetName(1).empty() ||
      ident.GetType(1) != netgen::Identifications::UNDEFINED)
    throw std::runtime_error("Pre-existing periodic metadata is unsupported");

  std::map<int, int> master_to_slave;
  std::set<int> slave_boundaries;
  for (const auto &patch : contract.boundary_pairs) {
    if (patch.first < 1 || patch.second < 1 || patch.first == patch.second ||
        !master_to_slave.emplace(patch.first, patch.second).second ||
        !slave_boundaries.insert(patch.second).second)
      throw std::runtime_error("Periodic patch map is not one-to-one");
  }
  for (const auto &patch : master_to_slave)
    if (slave_boundaries.count(patch.first))
      throw std::runtime_error("Periodic master and slave BC sets overlap");

  NamedPeriodicResult result;
  std::map<int, std::set<int>> boundary_points;
  std::map<int, std::set<int>> master_incidence;
  for (int fd = 1; fd <= mesh.GetNFD(); ++fd) {
    const int bc = mesh.GetFaceDescriptor(periodic_face_index(fd)).BCProperty();
    if (bc < 1) throw std::runtime_error("Invalid boundary-condition property");
    result.descriptor_to_boundary[fd] = bc;
    const auto &label = periodic_boundary_label(mesh, bc);
    if ((label == contract.master_label && !master_to_slave.count(bc)) ||
        (label == contract.slave_label && !slave_boundaries.count(bc)))
      throw std::runtime_error("Periodic patch map omits a labeled boundary");
  }
  // The container supplies native ElIndex values internally; no integer-to-
  // ElIndex construction or raw zero-/one-based reinterpretation occurs here.
  for (const auto &surface : mesh.SurfaceElements()) {
    const auto descriptor = result.descriptor_to_boundary.find(periodic_face_number(surface.GetIndex()));
    if (descriptor == result.descriptor_to_boundary.end())
      throw std::runtime_error("Surface element references an unknown descriptor");
    const int bc = descriptor->second;
    const bool master = master_to_slave.count(bc) != 0;
    if (!master && !slave_boundaries.count(bc)) continue;
    const auto &label = periodic_boundary_label(mesh, bc);
    if (label != (master ? contract.master_label : contract.slave_label))
      throw std::runtime_error("Periodic BC label mismatch at BC " + std::to_string(bc));
    for (int j = 0; j < surface.GetNP(); ++j) {
      const int pi = periodic_point_number(surface[j]);
      if (pi < 1 || pi > mesh.GetNP())
        throw std::runtime_error("Periodic surface point index out of range");
      const auto &point = mesh.Point(periodic_point_index(pi));
      for (int axis = 0; axis < 3; ++axis)
        if (!std::isfinite(point(axis)))
          throw std::runtime_error("Nonfinite periodic point coordinate");
      boundary_points[bc].insert(pi);
      if (master) master_incidence[pi].insert(bc);
    }
  }
  std::set<int> masters, slaves;
  for (const auto &patch : master_to_slave) {
    const auto &mp = boundary_points[patch.first];
    const auto &sp = boundary_points[patch.second];
    if (mp.empty() || sp.empty())
      throw std::runtime_error("Empty periodic patch at master BC " + std::to_string(patch.first));
    masters.insert(mp.begin(), mp.end()); slaves.insert(sp.begin(), sp.end());
  }
  if (masters.size() != slaves.size())
    throw std::runtime_error("Periodic point-set cardinality mismatch");
  for (int pi : masters)
    if (slaves.count(pi))
      throw std::runtime_error("Periodic point index belongs to both sides");

  std::set<int> assigned_slaves;
  for (int master : masters) {
    const auto &mp = mesh.Point(periodic_point_index(master));
    std::vector<std::pair<int, double>> candidates;
    for (int slave : slaves) {
      bool incident = true;
      for (int bc : master_incidence.at(master))
        if (!boundary_points.at(master_to_slave.at(bc)).count(slave)) incident = false;
      if (!incident) continue;
      const auto &sp = mesh.Point(periodic_point_index(slave));
      const double distance = std::hypot(
          sp(0) - (mp(0) + contract.translation[0]),
          sp(1) - (mp(1) + contract.translation[1]),
          sp(2) - (mp(2) + contract.translation[2]));
      if (distance <= contract.matching_distance)
        candidates.emplace_back(slave, distance);
    }
    if (candidates.size() != 1)
      throw std::runtime_error("Periodic master point " + std::to_string(master) +
                               " has " + std::to_string(candidates.size()) + " candidates");
    if (!assigned_slaves.insert(candidates[0].first).second)
      throw std::runtime_error("Periodic target point assigned more than once");
    result.point_pairs.emplace_back(master, candidates[0].first);
    result.maximum_distance = (std::max)(result.maximum_distance, candidates[0].second);
  }
  if (assigned_slaves != slaves)
    throw std::runtime_error("Incomplete periodic slave coverage");
  const std::map<int, int> mapping(result.point_pairs.begin(), result.point_pairs.end());
  for (const auto &patch : master_to_slave) {
    std::set<int> image;
    for (int master : boundary_points.at(patch.first)) image.insert(mapping.at(master));
    if (image != boundary_points.at(patch.second))
      throw std::runtime_error("Incomplete periodic patch coverage at master BC " +
                               std::to_string(patch.first));
  }

  // Validation precedes all mutations. Only identification metadata changes.
  const std::string geometry_before = periodic_geometry_serialization(mesh);
  for (const auto &pair : result.point_pairs)
    ident.Add(periodic_point_index(pair.first), periodic_point_index(pair.second), 1);
  ident.SetType(1, netgen::Identifications::PERIODIC);
  ident.SetName(1, contract.name);
  if (ident.GetMaxNr() != 1 || ident.GetName(1) != "kt" ||
      ident.GetType(1) != netgen::Identifications::PERIODIC)
    throw std::runtime_error("Named periodic metadata insertion failed");
  if (periodic_geometry_serialization(mesh) != geometry_before)
    throw std::runtime_error("Periodic insertion changed serialized geometry");
  return result;
}

} // namespace cubit_mesh_export
#endif
