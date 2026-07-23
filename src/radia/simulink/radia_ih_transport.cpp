#include "radia_ih_transport.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace radia { namespace ih {

void transport_periodic(const std::vector<double>& previous,
                        const std::vector<double>& weights,
                        double delta_angle_rad,
                        std::vector<double>& current) {
    const std::size_t n = previous.size();
    if (n == 0 || weights.size() != n)
        throw std::invalid_argument("IH transport requires equal non-empty field and weight vectors");
    double total = 0.0;
    for (double w : weights) {
        if (!(w > 0.0) || !std::isfinite(w))
            throw std::invalid_argument("IH transport weights must be finite and positive");
        total += w;
    }
    if (!(total > 0.0) || !std::isfinite(delta_angle_rad))
        throw std::invalid_argument("IH transport received an invalid angle or weight sum");

    // The configuration contract uses equally spaced periodic workpiece
    // samples.  Linear interpolation is conservative after the weighted
    // correction below and remains continuous for large wrapped rotations.
    const double period = 2.0 * std::acos(-1.0);
    double x = std::fmod(delta_angle_rad / period * static_cast<double>(n),
                         static_cast<double>(n));
    if (x < 0.0) x += static_cast<double>(n);
    current.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double source = static_cast<double>(i) - x;
        double wrapped = std::fmod(source, static_cast<double>(n));
        if (wrapped < 0.0) wrapped += static_cast<double>(n);
        const std::size_t j0 = static_cast<std::size_t>(wrapped);
        const std::size_t j1 = (j0 + 1) % n;
        const double a = wrapped - static_cast<double>(j0);
        current[i] = (1.0 - a) * previous[j0] + a * previous[j1];
    }
    double before = 0.0;
    double after = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        before += weights[i] * previous[i];
        after += weights[i] * current[i];
    }
    if (std::abs(after) > 0.0) {
        const double scale = before / after;
        for (double& value : current) value *= scale;
    }
}

bool eddy_needs_update(double current_now_A, double current_prev_A,
                       double relative_tolerance) {
    if (!std::isfinite(current_now_A) || !std::isfinite(current_prev_A) ||
        relative_tolerance < 0.0 || !std::isfinite(relative_tolerance))
        throw std::invalid_argument("invalid current or tolerance");
    const double scale = std::max({1.0, std::abs(current_now_A),
                                   std::abs(current_prev_A)});
    return std::abs(current_now_A - current_prev_A) > relative_tolerance * scale;
}

}}  // namespace radia::ih
