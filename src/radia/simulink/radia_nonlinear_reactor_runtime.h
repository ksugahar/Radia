#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace radia { namespace reactor {

struct Config {
    int n_modes = 0;
    int n_samples = 0;
    std::vector<double> demag;
    std::vector<double> magnetization_modes;
    std::vector<double> sample_weights;
    std::vector<double> excitation_per_amp;
    std::vector<double> magnetization_table_A_per_m;
    std::vector<double> field_table_A_per_m;
    double air_inductance_H = 0.0;
    double winding_resistance_Ohm = 0.0;
    double sample_time_s = 0.0;
    double initial_current_A = 0.0;
    double residual_tolerance = 1.0e-9;
    int max_iterations = 40;
    double line_search_minimum = 0x1p-20;
};

struct Output {
    double voltage_V = 0.0;
    double flux_linkage_Wb_turn = 0.0;
    double differential_inductance_H = 0.0;
    double peak_flux_density_T = 0.0;
    double magnetic_energy_J = 0.0;
    double residual_relative_norm = 0.0;
    int nonlinear_iterations = 0;
    std::vector<double> flux_density_T;
    std::vector<double> magnetization_coefficients;
};

struct Snapshot {
    double previous_current_A = 0.0;
    double previous_flux_linkage_Wb_turn = 0.0;
    std::uint64_t accepted_steps = 0;
    std::vector<double> magnetization_coefficients;
};

class Runtime {
public:
    explicit Runtime(Config config);

    const Output& output(double current_A);
    void update(double current_A);
    void reset();
    Snapshot snapshot() const;
    void restore(const Snapshot& state);

    int mode_count() const { return config_.n_modes; }
    int sample_count() const { return config_.n_samples; }
    std::uint64_t accepted_steps() const { return accepted_steps_; }

private:
    Output solve(double current_A, const std::vector<double>& initial) const;
    void validate() const;

    Config config_;
    std::vector<double> previous_magnetization_;
    double previous_current_A_ = 0.0;
    double previous_flux_linkage_ = 0.0;
    std::uint64_t accepted_steps_ = 0;
    bool cache_valid_ = false;
    double cached_current_A_ = 0.0;
    Output cached_output_;
};

}}  // namespace radia::reactor
