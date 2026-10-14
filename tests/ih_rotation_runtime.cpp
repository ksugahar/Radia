#include "radia_ih_runtime.h"
#include "radia_ih_transport.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace radia::ih;

// Independent material-frame oracle: a stationary source passes each material
// cell once per turn. With M=I, K=0 every cell gains dt*sum(q)/n per step.
int main() {
    std::vector<double> moved;
    bool rejected = false;
    try {
        transport_periodic({1, 1e-16, 0, -2}, {1, 2, 1, 2},
                           std::acos(-1.0)/2, moved);
    } catch (const std::invalid_argument&) { rejected = true; }
    if (!rejected) throw std::runtime_error("cancelling integral must fail fast");
    for (int subdivisions : {1, 2}) {
        const int n = 8, steps = n * subdivisions;
        EddyConfig e;
        e.n_unknown = 1; e.n_heat = n; e.n_temperature = n;
        e.matrix = {1.0}; e.rhs_per_amp = {1.0};
        e.heat_projection.assign(n, 0.0); e.heat_projection[0] = 1.0;
        e.heat_weights.assign(n, 1.0); e.temperature_weights.assign(n, 1.0);
        e.periodic_rotation = true;
        ThermalConfig t;
        t.n_heat = n; t.mass.n = n; t.stiffness.n = n;
        for (int i = 0; i <= n; ++i) t.mass.row_ptr.push_back(i);
        for (int i = 0; i < n; ++i) t.mass.col.push_back(i);
        t.mass.value.assign(n, 1.0); t.stiffness = t.mass;
        t.stiffness.value.assign(n, 0.0);
        t.initial_temperature_K.assign(n, 300.0); t.weights.assign(n, 1.0);
        t.heat_to_temperature.assign(n*n, 0.0);
        for (int i = 0; i < n; ++i) t.heat_to_temperature[i*n+i] = 1.0;
        t.options.dt_s = 1.0 / subdivisions; t.periodic_rotation = true;
        EddyRuntime eddy(e); ThermalRuntime thermal(t);
        for (int k = 0; k < steps; ++k) {
            const double angle = 2 * std::acos(-1.0) * k / steps;
            thermal.update(eddy.output(1.0, angle, thermal.output()), 300, angle);
        }
        for (double value : thermal.output()) {
            std::cout << value << ' ';
            if (std::abs(value - 301.0) > 1e-9)
                throw std::runtime_error("one turn must heat all material cells equally");
        }
        std::cout << '\n';
        // Rotation alone cannot move a temperature field in material coordinates.
        t.initial_temperature_K[0] = 310;
        ThermalRuntime passive(t);
        passive.update(std::vector<double>(n, 0.0), 300, 0.37);
        if (passive.output() != t.initial_temperature_K)
            throw std::runtime_error("passive material temperature moved");
        eddy.reset();
        const auto heat = eddy.output(1, 0, t.initial_temperature_K);
        if (heat != e.heat_projection) throw std::runtime_error("reset failed");
    }
}
