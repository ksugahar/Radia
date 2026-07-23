#pragma once

#include <cstddef>
#include <vector>

namespace radia { namespace ih {

// Transport a scalar field stored in workpiece coordinates by a cyclic
// angular shift.  The operation is conservative with respect to the supplied
// positive cell weights; it is used by the Thermal S-Function between two
// mechanical angles.
void transport_periodic(const std::vector<double>& previous,
                        const std::vector<double>& weights,
                        double delta_angle_rad,
                        std::vector<double>& current);

bool eddy_needs_update(double current_now_A, double current_prev_A,
                       double relative_tolerance = 0.0);

}}  // namespace radia::ih
