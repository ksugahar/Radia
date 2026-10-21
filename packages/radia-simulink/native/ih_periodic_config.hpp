// Typed periodic-drive fragment for CoupledConfig; review candidate, not full ABI.
#pragma once
#include "ih_coupled_domain.hpp"
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace radia_ih {
struct PeriodicDriveConfig {
    double frequency_Hz;
    double thermal_macro_dt_s;
    std::uint32_t micro_steps_per_cycle;
    std::uint32_t max_settling_cycles;
    double periodic_relative_limit;
    double energy_relative_limit;
    double integer_cycle_absolute_tolerance;
    // Fixed physical mode: reversible isotropic odd single-valued BH(T), no
    // hysteresis; previous accepted T held through all trial magnetic cycles.
    // Fixed integration: nonlinear implicit midpoint, phase-zero endpoints,
    // final ONE settled cycle averaged, re-settle on amplitude/T changes.
    std::uint32_t checked_periods_per_macro() const {
        if (!std::isfinite(frequency_Hz) || frequency_Hz <= 0. ||
            !std::isfinite(thermal_macro_dt_s) || thermal_macro_dt_s <= 0. ||
            micro_steps_per_cycle < 8 || micro_steps_per_cycle > 65536 ||
            max_settling_cycles < 2 || max_settling_cycles > 10000 ||
            !std::isfinite(periodic_relative_limit) ||
            periodic_relative_limit <= 0. || periodic_relative_limit > 1e-6 ||
            !std::isfinite(energy_relative_limit) ||
            energy_relative_limit <= 0. || energy_relative_limit > 1e-6 ||
            !std::isfinite(integer_cycle_absolute_tolerance) ||
            integer_cycle_absolute_tolerance <= 0. ||
            integer_cycle_absolute_tolerance > 1e-10)
            throw std::invalid_argument("Invalid IH periodic drive configuration");
        const double count = frequency_Hz * thermal_macro_dt_s;
        if (!std::isfinite(count) || count < 1. || count > 1000000.)
            throw std::invalid_argument("IH macro period count outside declared bounds");
        const double nearest = std::round(count);
        if (std::abs(count-nearest) > integer_cycle_absolute_tolerance)
            throw std::invalid_argument("IH macro dt must span integer drive periods");
        // No physical timestep is rounded/rewritten; integer is metadata only.
        return static_cast<std::uint32_t>(nearest);
    }
};
inline void validate_controls(const Controls& value) {
    if (!std::isfinite(value.current_A) || !std::isfinite(value.ambient_K) ||
        value.ambient_K <= 0.)
        throw std::invalid_argument("IH controls must be finite with positive Kelvin ambient");
    switch (value.drive_mode) {
    case DriveMode::sinusoidal_periodic:
        if (value.current_A <= 0.)
            throw std::invalid_argument("IH sinusoidal drive requires positive peak current");
        break;
    case DriveMode::zero_drive_cooling:
        if (value.current_A != 0.)
            throw std::invalid_argument("IH tagged cooling mode requires zero current");
        break;
    default: throw std::invalid_argument("Unknown IH drive mode");
    }
}
} // namespace radia_ih
