// Guards: geometric-edge identity survives all native callback entry points,
// including two curves with the same faces and a smooth spherical seam.
// Assumes: source contract with explicit geometry callbacks, no Cubit license;
// real RefEdge/ACIS integration remains the licensed export acceptance lane.
#include "callbackgeom.hpp"
#include "CurveIdentity.hpp"
#include "NetgenCurver.hpp"
#include "named_periodic_contract.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <string>

using namespace netgen;
static void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
template<class F> static void rejected(F action) {
  bool failed = false;
  try { action(); } catch (const std::invalid_argument&) { failed = true; }
  require(failed, "invalid identity did not fail loudly");
}

// Explicit qualification mode only: no input generation, Curve, refit or defaults.
static int q3_input_contract(int argc, char **argv) {
  const std::set<std::string> required{
      "--q3-input", "--q3-output-dir", "--name", "--master", "--slave", "--unit",
      "--tx", "--ty", "--tz", "--matching-distance", "--patches"};
  std::map<std::string, std::string> args;
  require(argc == 1 + 2 * static_cast<int>(required.size()), "All q3 input options are required");
  for (int i = 1; i < argc; i += 2) {
    require(required.count(argv[i]) != 0, "Unknown q3 input option");
    require(*argv[i+1] != '\0' && args.emplace(argv[i], argv[i+1]).second,
            "Empty or duplicate q3 input option");
  }
  auto number = [](const std::string &text) {
    size_t used = 0;
    double value = std::stod(text, &used);
    require(used == text.size() && std::isfinite(value), "Invalid explicit floating-point value");
    return value;
  };
  auto integer = [](const std::string &text) {
    require(!text.empty() && text.find_first_not_of("0123456789") == std::string::npos,
            "Invalid explicit boundary number");
    size_t used = 0;
    int value = std::stoi(text, &used);
    require(used == text.size() && value > 0, "Invalid explicit boundary number");
    return value;
  };
  cubit_mesh_export::NamedPeriodicContract contract;
  contract.name = args.at("--name");
  contract.master_label = args.at("--master");
  contract.slave_label = args.at("--slave");
  contract.coordinate_unit = args.at("--unit");
  contract.translation = {number(args.at("--tx")), number(args.at("--ty")), number(args.at("--tz"))};
  contract.matching_distance = number(args.at("--matching-distance"));
  const auto &patches = args.at("--patches");
  size_t start = 0;
  while (start < patches.size()) {
    const auto end = patches.find(',', start);
    const auto row = patches.substr(start, end == std::string::npos ? end : end - start);
    const auto colon = row.find(':');
    require(colon != std::string::npos && row.find(':', colon + 1) == std::string::npos,
            "Expected explicit masterBC:slaveBC patch map");
    contract.boundary_pairs.emplace_back(integer(row.substr(0, colon)), integer(row.substr(colon + 1)));
    if (end == std::string::npos) break;
    start = end + 1;
    require(start < patches.size(), "Trailing empty patch map entry");
  }
  const std::filesystem::path input(args.at("--q3-input")), output(args.at("--q3-output-dir"));
  require(input.is_absolute() && output.is_absolute(), "Qualification paths must be absolute");
  require(std::filesystem::is_regular_file(input), "q3 input must be an existing file");
  require(!std::filesystem::exists(output), "Qualification output directory must be fresh");
  require(std::filesystem::is_directory(output.parent_path()), "Output parent must exist");
  Mesh mesh;
  mesh.Load(input);
  require(mesh.GetNE() > 0 && mesh.GetCurvedElements().GetOrder() == 3 &&
          mesh.GetCurvedElements().IsHighOrder(), "Existing curved q3 volume mesh required; no Curve is allowed");
  auto serialize = [&] { std::ostringstream stream; mesh.Save(stream); return stream.str(); };
  const auto before = serialize();
  const auto geometry_before = cubit_mesh_export::periodic_geometry_serialization(mesh);
  require(std::filesystem::create_directory(output), "Cannot create fresh output directory");
  auto write = [&](const char *name, const std::string &data) {
    std::ofstream stream(output / name, std::ios::binary | std::ios::out);
    require(static_cast<bool>(stream), "Cannot open qualification output");
    stream.write(data.data(), static_cast<std::streamsize>(data.size()));
    stream.close();
    require(static_cast<bool>(stream), "Cannot write qualification output");
  };
  write("before.vol", before);
  write("before.geometry", geometry_before);
  const auto result = cubit_mesh_export::apply_named_periodic_identification(mesh, contract);
  const auto after = serialize();
  const auto geometry_after = cubit_mesh_export::periodic_geometry_serialization(mesh);
  write("after.vol", after);
  write("after.geometry", geometry_after);
  require(geometry_before == geometry_after, "q3 input geometry serialization changed");
  std::ostringstream table;
  table << "fd\tbc\tlabel\n";
  for (const auto &entry : result.descriptor_to_boundary) {
    const auto label = cubit_mesh_export::periodic_boundary_label(mesh, entry.second);
    require(label.find_first_of("\t\r\n") == std::string::npos, "Unsupported control character in label");
    table << entry.first << '\t' << entry.second << '\t' << label << '\n';
  }
  write("fd-bc-labels.tsv", table.str());
  std::ostringstream report;
  report << std::setprecision(17)
         << "{\"status\":\"passed\",\"curve_order\":3,\"geometry_serialization_equal\":true,"
         << "\"point_count\":" << mesh.GetNP() << ",\"volume_element_count\":" << mesh.GetNE()
         << ",\"pair_count\":" << result.point_pairs.size()
         << ",\"maximum_distance\":" << result.maximum_distance << ",\"pairs\":[";
  bool first = true;
  for (const auto &pair : result.point_pairs) {
    if (!first) report << ',';
    report << '[' << pair.first << ',' << pair.second << ",1]";
    first = false;
  }
  report << "]}\n";
  write("result.json", report.str());
  std::cout << report.str();
  return 0;
}

int main(int argc, char **argv) {
  try {
    if (argc != 1) return q3_input_contract(argc, argv);
    named_periodic_contract_test();
    Mesh mesh;
    const int first = mesh.AddEdgeDescriptor(EdgeDescriptor(31, 3, 4)).Nr1();
    const int second = mesh.AddEdgeDescriptor(EdgeDescriptor(32, 3, 4)).Nr1();
    auto identities = NetgenCurver::curve_identities(mesh, {{101, first}, {202, second}});
    require(identities.lookup(31, 3, 4).curve == 101, "production mapping lost first curve");
    require(identities.lookup(32, 3, 4).curve == 202, "production mapping lost second curve");
    rejected([&] { identities.add(31, {9, 303, 3, 4}); });
    rejected([&] { identities.lookup(-1, 3, 4); });
    rejected([&] { identities.lookup(31, 3, 5); });
    rejected([&] { NetgenCurver::curve_identities(mesh, {{101, 0}}); });
    rejected([&] { NetgenCurver::curve_identities(mesh, {{101, first}, {202, first}}); });
    auto surface = [](int, double x, double y, double z, double, double, bool) {
      return std::make_tuple(x, y, z, 0.0, 0.0);
    };
    auto normal = [](int, double x, double y, double z) {
      return std::make_tuple(x, y, z); // same smooth sphere on BOTH patches
    };
    auto project = [&](int edge, int s1, int s2, double, double, double z) {
      const auto & b = identities.lookup(edge, s1, s2);
      double x = std::sqrt(1.0-z*z) * (b.curve == 101 ? 1.0 : -1.0);
      return std::make_tuple(x, 0.0, z);
    };
    auto between = [&](int edge, int s1, int s2, double, double, double z1,
                       double, double, double z2, double t) {
      return project(edge, s1, s2, 0, 0, z1+t*(z2-z1));
    };
    auto tangent = [&](int edge, int s1, int s2, double x, double, double z) {
      identities.lookup(edge, s1, s2);
      return std::make_tuple(z, 0.0, -x);
    };
    CallbackGeometry geometry(surface, normal, 4, tangent, project, between);
    CallbackGeometry projected_midpoint(surface, normal, 4, tangent, project);
    int checked = 0;
    for (int edge : {31, 32}) {
      const double sign = edge == 31 ? 1.0 : -1.0;
      Point<3> p(0, 0, 0), midpoint;
      EdgePointGeomInfo info, next;
      geometry.ProjectPointEdge(3, 4, p, &info, edge);
      require(p[0] == sign && p[1] == 0 && p[2] == 0, "projection lost curve identity");
      auto t = geometry.GetTangent(p, 4, 3, info, edge);
      require(t.Length() == 1.0 && t[2] == -sign, "smooth seam tangent was lost");
      geometry.PointBetweenEdge(Point<3>(0,0,-1), Point<3>(0,0,1),
                                0.5, 3, 4, info, info, midpoint, next, edge);
      require(midpoint[0] == sign, "interpolation lost curve identity");
      projected_midpoint.PointBetweenEdge(Point<3>(0,0,-1), Point<3>(0,0,1),
                                0.5, 3, 4, info, info, midpoint, next, edge);
      require(midpoint[0] == sign, "projected midpoint lost curve identity");
      ++checked;
    }
    rejected([&] { Point<3> p(0,0,0); geometry.ProjectPointEdge(3,4,p,nullptr,99); });
    std::cout << "{\"status\":\"passed\",\"distinct_curves_same_faces\":" << checked
              << ",\"smooth_seam\":true,\"unknown_identity_rejected\":true}\n";
    return 0;
  } catch (const std::exception & e) {
    std::cerr << "native curve identity contract: " << e.what() << '\n';
    return 1;
  }
}
