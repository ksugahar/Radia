#ifndef CUBIT_MESH_EXPORT_CURVE_IDENTITY_HPP
#define CUBIT_MESH_EXPORT_CURVE_IDENTITY_HPP

#include <map>
#include <stdexcept>
#include <string>

namespace cubit_mesh_export {
// The geometric EdgeNr is distinct from the 1-based descriptor table index.
// Surface pairs are deliberately NOT keys: split smooth faces share >1 curve.
struct CurveBinding {
  int descriptor, curve, face1, face2;
};

class CurveIdentity {
  std::map<int, CurveBinding> bindings;
public:
  void add(int edge, CurveBinding binding) {
    if (edge < 1 || binding.descriptor < 1 || binding.curve < 1 ||
        binding.face1 < 1 || binding.face2 < 1)
      throw std::invalid_argument("Invalid Netgen/Cubit curve identity: edge=" + std::to_string(edge) +
          " descriptor=" + std::to_string(binding.descriptor) + " curve=" + std::to_string(binding.curve) +
          " faces=" + std::to_string(binding.face1) + "," + std::to_string(binding.face2));
    if (!bindings.emplace(edge, binding).second)
      throw std::invalid_argument("Ambiguous Netgen geometric edge " + std::to_string(edge));
  }
  const CurveBinding & lookup(int edge, int face1, int face2) const {
    auto it = bindings.find(edge);
    if (it == bindings.end())
      throw std::invalid_argument("Unknown Netgen geometric edge " + std::to_string(edge));
    const auto & b = it->second;
    if (!((b.face1 == face1 && b.face2 == face2) ||
          (b.face1 == face2 && b.face2 == face1)))
      throw std::invalid_argument("Surface mismatch for Netgen geometric edge " + std::to_string(edge));
    return b;
  }
};
}
#endif
