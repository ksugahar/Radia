// Guards: geometric-edge identity survives all native callback entry points,
// including two curves with the same faces and a smooth spherical seam.
// Assumes: source contract with explicit geometry callbacks, no Cubit license;
// real RefEdge/ACIS integration remains the licensed export acceptance lane.
#include "callbackgeom.hpp"
#include "CurveIdentity.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace netgen;
static void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
template<class F> static void rejected(F action) {
  bool failed = false;
  try { action(); } catch (const std::invalid_argument&) { failed = true; }
  require(failed, "invalid identity did not fail loudly");
}

int main() {
  try {
    cubit_mesh_export::CurveIdentity identities;
    identities.add(31, {7, 101, 3, 4});
    identities.add(32, {8, 202, 3, 4});
    rejected([&] { identities.add(31, {9, 303, 3, 4}); });
    rejected([&] { identities.lookup(-1, 3, 4); });
    rejected([&] { identities.lookup(31, 3, 5); });
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
