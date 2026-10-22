// Host-only adapter-policy contract; no domain construction or GPU execution.
#include "../../native/ih_coupled_runtime.hpp"
#include <limits>
#include <stdexcept>

int main() {
    using radia_ih::Controls;
    using radia_ih::DriveMode;
    using Policy = radia_ih_runtime::CoupledPolicy;
    const Controls valid{2., 293.15, DriveMode::sinusoidal_periodic};
    Policy::validate(valid);
    Policy::validate(Controls{0., 293.15, DriveMode::zero_drive_cooling});
    const Controls zero_sine{0., 293.15, DriveMode::sinusoidal_periodic};
    const Controls cooling{0., 293.15, DriveMode::zero_drive_cooling};
    if (!Policy::equal(valid, valid) ||
        Policy::equal(valid, Controls{3., 293.15, DriveMode::sinusoidal_periodic}) ||
        Policy::equal(valid, Controls{2., 300., DriveMode::sinusoidal_periodic}) ||
        Policy::equal(zero_sine, cooling))
        throw std::runtime_error("IH sample identity omitted a control");
    const Controls invalid[] = {
        {0., 293.15, DriveMode::sinusoidal_periodic},
        {-1., 293.15, DriveMode::sinusoidal_periodic},
        {1., 293.15, DriveMode::zero_drive_cooling},
        {1., 0., DriveMode::sinusoidal_periodic},
        {1., -1., DriveMode::sinusoidal_periodic},
        {1., 293.15, static_cast<DriveMode>(99)},
        {std::numeric_limits<double>::quiet_NaN(), 293.15, DriveMode::sinusoidal_periodic},
        {std::numeric_limits<double>::infinity(), 293.15, DriveMode::sinusoidal_periodic},
        {1., std::numeric_limits<double>::quiet_NaN(), DriveMode::sinusoidal_periodic},
        {1., std::numeric_limits<double>::infinity(), DriveMode::sinusoidal_periodic}
    };
    for (const auto& controls : invalid) {
        bool rejected = false;
        try { Policy::validate(controls); }
        catch (const std::invalid_argument&) { rejected = true; }
        if (!rejected) throw std::runtime_error("Invalid IH controls accepted");
    }
}
