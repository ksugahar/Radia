#ifndef CUBIT_NAMED_PERIODIC_IDENTIFICATION_HPP
#define CUBIT_NAMED_PERIODIC_IDENTIFICATION_HPP

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include "NamedPeriodicMeshAdapter.hpp"
#include "PeriodicCurvedGeometry.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cerrno>
#include <system_error>
#include <filesystem>
#include <fstream>
#include <iterator>
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
  std::vector<PeriodicCurvedTransfer> curved_transfers;
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

  const auto triangles = periodic_triangle_plan(mesh, master_to_slave, mapping);
  result.curved_transfers = synchronize_periodic_curvature(mesh, triangles);
  // Curvature synchronization is intentional. ID insertion itself must still
  // preserve the already synchronized geometric serialization exactly.
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

struct NamedPeriodicOutputPaths {
  std::filesystem::path final, temporary, sidecar;
};

inline NamedPeriodicOutputPaths named_periodic_output_paths(const std::filesystem::path &filename)
{
  if (filename.extension() != ".vol")
    throw std::runtime_error("Named periodic output requires an explicit .vol suffix");
  auto temporary = filename; temporary += ".partial.vol";
  auto sidecar = filename; sidecar += ".json";
  return {filename, temporary, sidecar};
}

inline void require_fresh_named_periodic_output(const NamedPeriodicOutputPaths &paths)
{
  for (const auto &path : {paths.final, paths.temporary, paths.sidecar})
    if (std::filesystem::symlink_status(path).type() != std::filesystem::file_type::not_found)
      throw std::runtime_error("Named periodic output already exists: " + path.u8string() +
          ". Inspect the earlier failure; remove only your own unqualified leftovers before retrying.");
}

// No replace: the OS checks the destination and publishes atomically. A
// pre-publication exists() check alone cannot close the Windows rename race.
inline void move_named_periodic_no_replace(const std::filesystem::path &source,
                                           const std::filesystem::path &target)
{
#ifdef _WIN32
  if (!MoveFileExW(source.c_str(), target.c_str(), 0))
    throw std::system_error(static_cast<int>(GetLastError()), std::system_category(),
                            "Exclusive periodic publication failed");
#else
  std::filesystem::create_hard_link(source, target); // atomic no-replace
  std::filesystem::remove(source);
#endif
}

inline std::string named_periodic_retry_guidance(const NamedPeriodicOutputPaths &paths)
{
  return " Failed outputs are unqualified. Inspect the original error and remove only "
         "your own leftover files before a fresh attempt: " + paths.final.u8string() +
         ", " + paths.temporary.u8string() + ", " + paths.sidecar.u8string();
}

inline void publish_named_periodic_output(const NamedPeriodicOutputPaths &paths)
{
  if (std::filesystem::symlink_status(paths.sidecar).type() != std::filesystem::file_type::not_found)
    throw std::runtime_error("Periodic sidecar target already exists");
  move_named_periodic_no_replace(paths.temporary, paths.final);
}

// Serialize the complete record before exclusively creating the sidecar. Never
// truncate a pre-existing file; flush/close failures are export failures too.
struct WrittenPeriodicSidecar {
  std::filesystem::path path;
  std::string bytes;
  std::array<std::uint64_t, 3> identity{};
};

template <typename Writer>
inline WrittenPeriodicSidecar write_named_periodic_sidecar(const std::filesystem::path &path, Writer writer)
{
  std::ostringstream record;
  writer(record);
  if (!record) throw std::runtime_error("Periodic sidecar serialization failed");
  const std::string bytes = record.str();
  if (bytes.empty()) throw std::runtime_error("Empty periodic sidecar record");
  WrittenPeriodicSidecar receipt{path, bytes, {}};
#ifdef _WIN32
  HANDLE file = CreateFileW(path.c_str(), GENERIC_WRITE, FILE_SHARE_READ, nullptr,
                            CREATE_NEW, FILE_ATTRIBUTE_NORMAL, nullptr);
  if (file == INVALID_HANDLE_VALUE)
    throw std::system_error(static_cast<int>(GetLastError()), std::system_category(),
                            "Cannot exclusively open periodic sidecar");
  DWORD error = ERROR_SUCCESS;
  std::size_t offset = 0;
  while (offset < bytes.size()) {
    const DWORD count = static_cast<DWORD>(std::min<std::size_t>(bytes.size() - offset, 0x7fffffff));
    DWORD written = 0;
    if (!WriteFile(file, bytes.data() + offset, count, &written, nullptr) || written != count) {
      error = GetLastError();
      if (error == ERROR_SUCCESS) error = ERROR_WRITE_FAULT;
      break;
    }
    offset += written;
  }
  if (error == ERROR_SUCCESS && !FlushFileBuffers(file)) {
    error = GetLastError();
    if (error == ERROR_SUCCESS) error = ERROR_WRITE_FAULT;
  }
  BY_HANDLE_FILE_INFORMATION identity{};
  if (error == ERROR_SUCCESS && !GetFileInformationByHandle(file, &identity)) {
    error = GetLastError();
    if (error == ERROR_SUCCESS) error = ERROR_READ_FAULT;
  }
  receipt.identity = {identity.dwVolumeSerialNumber, identity.nFileIndexHigh, identity.nFileIndexLow};
  if (!CloseHandle(file) && error == ERROR_SUCCESS) {
    error = GetLastError();
    if (error == ERROR_SUCCESS) error = ERROR_WRITE_FAULT;
  }
  if (error != ERROR_SUCCESS)
    throw std::system_error(static_cast<int>(error), std::system_category(),
                            "Periodic sidecar write/flush/close failed");
#else
  const int file = ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0600);
  if (file < 0) throw std::system_error(errno, std::generic_category(), "Cannot exclusively open periodic sidecar");
  int error = 0;
  std::size_t offset = 0;
  while (offset < bytes.size()) {
    const auto written = ::write(file, bytes.data() + offset, bytes.size() - offset);
    if (written < 0 && errno == EINTR) continue;
    if (written <= 0) { error = written < 0 ? errno : EIO; break; }
    offset += static_cast<std::size_t>(written);
  }
  if (!error && ::fsync(file) != 0) error = errno;
  struct stat identity{};
  if (!error && ::fstat(file, &identity) != 0) error = errno;
  receipt.identity = {std::uint64_t(identity.st_dev), std::uint64_t(identity.st_ino), 0};
  if (::close(file) != 0 && !error) error = errno;
  if (error) throw std::system_error(error, std::generic_category(), "Periodic sidecar write/flush/close failed");
#endif
  return receipt;
}

// Bundle finalization accepts only the exact exclusive record from this run.
// Keep a read handle open across publication; on Windows its sharing mode
// prevents another process from modifying/removing/replacing that record.
inline void publish_named_periodic_bundle(const NamedPeriodicOutputPaths &paths,
                                          const WrittenPeriodicSidecar &receipt)
{
  if (receipt.path != paths.sidecar || receipt.bytes.empty())
    throw std::runtime_error("Periodic bundle record ownership mismatch");
  std::string actual;
#ifdef _WIN32
  HANDLE file = CreateFileW(paths.sidecar.c_str(), GENERIC_READ, FILE_SHARE_READ,
      nullptr, OPEN_EXISTING, FILE_FLAG_OPEN_REPARSE_POINT, nullptr);
  if (file == INVALID_HANDLE_VALUE)
    throw std::system_error(static_cast<int>(GetLastError()), std::system_category(), "Cannot verify owned periodic record");
  struct Close { HANDLE file; ~Close() { CloseHandle(file); } } close{file};
  BY_HANDLE_FILE_INFORMATION info{};
  if (!GetFileInformationByHandle(file, &info) || (info.dwFileAttributes & (FILE_ATTRIBUTE_REPARSE_POINT | FILE_ATTRIBUTE_DIRECTORY)) ||
      receipt.identity != std::array<std::uint64_t,3>{info.dwVolumeSerialNumber, info.nFileIndexHigh, info.nFileIndexLow})
    throw std::runtime_error("Periodic bundle record identity changed");
  char buffer[4096]; DWORD count = 0;
  do {
    if (!ReadFile(file, buffer, sizeof(buffer), &count, nullptr))
      throw std::runtime_error("Cannot read owned periodic record");
    actual.append(buffer, count);
    if (actual.size() > receipt.bytes.size())
      throw std::runtime_error("Periodic bundle record content changed");
  } while (count);
#else
  const int file = ::open(paths.sidecar.c_str(), O_RDONLY | O_NOFOLLOW);
  if (file < 0) throw std::system_error(errno, std::generic_category(), "Cannot verify owned periodic record");
  struct Close { int file; ~Close() { ::close(file); } } close{file};
  struct stat info{};
  if (::fstat(file, &info) != 0 || !S_ISREG(info.st_mode) ||
      receipt.identity != std::array<std::uint64_t,3>{std::uint64_t(info.st_dev), std::uint64_t(info.st_ino), 0})
    throw std::runtime_error("Periodic bundle record identity changed");
  char buffer[4096];
  for (;;) {
    const auto count = ::read(file, buffer, sizeof(buffer));
    if (count < 0 && errno == EINTR) continue;
    if (count < 0) throw std::runtime_error("Cannot read owned periodic record");
    if (!count) break;
    actual.append(buffer, count);
    if (actual.size() > receipt.bytes.size())
      throw std::runtime_error("Periodic bundle record content changed");
  }
  struct stat current{};
  if (::lstat(paths.sidecar.c_str(), &current) != 0 ||
      current.st_dev != info.st_dev || current.st_ino != info.st_ino)
    throw std::runtime_error("Periodic bundle record path changed");
#endif
  if (actual != receipt.bytes)
    throw std::runtime_error("Periodic bundle record content changed");
  move_named_periodic_no_replace(paths.temporary, paths.final);
}

// Single-owner output: preserve the native .vol text writer and its newline
// convention. Unknown suffixes cause Mesh::Save(path) to append .vol.gz.
inline NamedPeriodicOutputPaths stage_named_periodic_vol(
    const netgen::Mesh &mesh, const std::filesystem::path &filename)
{
  const auto paths = named_periodic_output_paths(filename);
  require_fresh_named_periodic_output(paths);
  std::ostringstream expected;
  mesh.Save(expected);
  if (!expected) throw std::runtime_error("Named periodic serialization failed");
  mesh.Save(paths.temporary);
  if (!std::filesystem::is_regular_file(paths.temporary))
    throw std::runtime_error("Native writer did not create the exact periodic temporary path");
  // Text input reverses the platform newline translation of the native text
  // writer, without changing the on-disk bytes. A short/silent write must fail.
  std::ifstream input(paths.temporary);
  input.exceptions(std::ios::badbit);
  if (!input) throw std::runtime_error("Cannot read periodic temporary output");
  const std::string actual((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  input.close();
  if (actual != expected.str())
    throw std::runtime_error("Native periodic temporary serialization is incomplete or differs");
  return paths;
}

inline void save_named_periodic_vol(const netgen::Mesh &mesh, const std::filesystem::path &filename)
{
  publish_named_periodic_output(stage_named_periodic_vol(mesh, filename));
}

// The final VOL is the publication marker. A failed mandatory record must
// never leave that marker visible. Owned partials remain explicit diagnostics.
template <class Writer>
inline void save_named_periodic_bundle(const netgen::Mesh &mesh,
    const std::filesystem::path &filename, Writer writer)
{
  const auto paths = stage_named_periodic_vol(mesh, filename);
  const auto record = write_named_periodic_sidecar(paths.sidecar, writer);
  publish_named_periodic_bundle(paths, record);
}

} // namespace cubit_mesh_export
#endif
