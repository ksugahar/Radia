#pragma once
// Copy curved edge/face coefficients across periodic seams.
// Moved out of the Radia Netgen fork (CurvedElements::CopyPeriodicCoefficients).
#include <array>
#include <vector>

namespace netgen { class Mesh; }

namespace radia_netgen_curving {

struct PeriodicTransfer
{
  bool face;
  int master, slave, mode, sign;   // 0-based topology edge/face numbers
  std::array<double, 3> original, expected;
};

// master[i] and slave[i] are triangles given by 1-based point numbers, with
// slave[i][k] the periodic image of master[i][k].  Copies every slave edge and
// face coefficient from its master (odd edge modes follow the orientation)
// and returns the applied transfers.  Polynomial orders 1..3 only.  On any
// failure the coefficients are restored before the exception propagates.
std::vector<PeriodicTransfer> copy_periodic_coefficients(
    netgen::Mesh & mesh,
    const std::vector<std::array<int, 3>> & master,
    const std::vector<std::array<int, 3>> & slave);

} // namespace radia_netgen_curving
