// NetgenCurver::curve_identities: the Netgen edge -> Cubit curve mapping.
// Kept free of Cubit headers so the license-free contract tests link the same
// production code as radia_cubit_export.dll.
#include "NetgenCurver.hpp"
#include "netgen_meshing.hpp"
#include "CurveIdentity.hpp"

radia_cubit_export::CurveIdentity NetgenCurver::curve_identities(
    const netgen::Mesh & mesh, const std::unordered_map<int, int> & curve_descriptors)
{
  radia_cubit_export::CurveIdentity bindings;
  for (const auto & entry : curve_descriptors) {
    const int descriptor = entry.second;
    if (descriptor < 1 || descriptor > mesh.EdgeDescriptors().Size())
      throw std::invalid_argument("Cubit curve " + std::to_string(entry.first) +
          " has invalid Netgen descriptor " + std::to_string(descriptor));
    const auto & ed = mesh.GetEdgeDescriptor(descriptor);
    bindings.add(ed.EdgeNr(), {descriptor, entry.first, ed.SurfNr(0), ed.SurfNr(1)});
  }
  return bindings;
}
