#pragma once

#include "radia_ih_thermal.h"

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace radia { namespace ih {

struct EddyConfig {
    int n_unknown = 0;
    int n_heat = 0;
    int n_temperature = 0;
    std::vector<std::complex<double>> matrix;
    std::vector<std::complex<double>> matrix_temperature_slope;
    double reference_temperature_K = 293.15;
    std::vector<std::complex<double>> rhs_per_amp;
    std::vector<double> heat_projection;
    std::vector<double> heat_weights;
    std::vector<double> temperature_weights;
    bool periodic_rotation = false;
    double angle_origin_rad = 0.0;
};

class EddyRuntime {
public:
    explicit EddyRuntime(EddyConfig config);
    std::vector<double> output(double current_A, double angle_rad,
                               const std::vector<double>& temperature_K);
    int heat_size() const { return config_.n_heat; }
    int temperature_size() const {
        return static_cast<int>(previous_temperature_.size());
    }
private:
    EddyConfig config_;
    std::vector<double> cached_heat_;
    std::vector<double> reference_heat_;
    std::vector<double> previous_temperature_;
    bool have_cache_ = false;
};

struct ThermalConfig {
    int n_heat = 0;
    CSRMatrix mass;
    CSRMatrix stiffness;
    CSRMatrix convection;
    bool has_convection = false;
    std::vector<double> initial_temperature_K;
    std::vector<double> weights;
    std::vector<double> heat_to_temperature;
    ThermalStepOptions options;
    bool periodic_rotation = false;
    double angle_origin_rad = 0.0;
};

class ThermalRuntime {
public:
    explicit ThermalRuntime(ThermalConfig config);
    const std::vector<double>& output() const { return state_.temperature_K; }
    void update(const std::vector<double>& heat_W_per_m3,
                double ambient_temperature_K, double angle_rad);
    void reset();
    std::size_t step_count() const { return step_count_; }
private:
    ThermalConfig config_;
    ThermalState state_;
    std::size_t step_count_ = 0;
};

}}  // namespace radia::ih
