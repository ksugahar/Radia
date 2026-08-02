#include "radia_ih_runtime.h"
#include "radia_ih_transport.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace radia { namespace ih {
namespace {

using Complex = std::complex<double>;

bool finite(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

bool equivalent_angle(double angle, double origin) {
    const double two_pi = 2.0 * std::acos(-1.0);
    return std::abs(std::remainder(angle - origin, two_pi)) <= 1.0e-12;
}

void require_finite(const std::vector<double>& values, const char* name) {
    if (!std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); }))
        throw std::invalid_argument(std::string(name) +
                                    " must contain finite values");
}

void require_positive(const std::vector<double>& values, const char* name) {
    if (!std::all_of(values.begin(), values.end(), [](double value) {
            return value > 0.0 && std::isfinite(value);
        }))
        throw std::invalid_argument(std::string(name) +
                                    " must contain finite positive values");
}

void validate_csr(const CSRMatrix& matrix, int n, const char* name) {
    if (matrix.n != n || matrix.row_ptr.size() !=
            static_cast<std::size_t>(n + 1) || matrix.row_ptr.empty() ||
        matrix.row_ptr.front() != 0 || matrix.col.size() != matrix.value.size() ||
        matrix.row_ptr.back() != static_cast<int>(matrix.col.size()))
        throw std::invalid_argument(std::string("invalid ") + name +
                                    " CSR matrix");
    for (int row = 0; row < n; ++row) {
        if (matrix.row_ptr[static_cast<std::size_t>(row)] >
            matrix.row_ptr[static_cast<std::size_t>(row + 1)])
            throw std::invalid_argument(std::string("non-monotone ") + name +
                                        " CSR rows");
    }
    for (std::size_t index = 0; index < matrix.col.size(); ++index) {
        if (matrix.col[index] < 0 || matrix.col[index] >= n ||
            !std::isfinite(matrix.value[index]))
            throw std::invalid_argument(std::string("invalid ") + name +
                                        " CSR entry");
    }
}

std::vector<Complex> solve(std::vector<Complex> matrix,
                           std::vector<Complex> rhs, int n) {
    for (int column = 0; column < n; ++column) {
        int pivot = column;
        double best = std::abs(matrix[static_cast<std::size_t>(column * n + column)]);
        for (int row = column + 1; row < n; ++row) {
            const double value =
                std::abs(matrix[static_cast<std::size_t>(row * n + column)]);
            if (value > best) {
                best = value;
                pivot = row;
            }
        }
        if (!(best > 0.0) || !std::isfinite(best))
            throw std::runtime_error("IH Eddy operator is singular");
        if (pivot != column) {
            for (int entry = column; entry < n; ++entry)
                std::swap(
                    matrix[static_cast<std::size_t>(column * n + entry)],
                    matrix[static_cast<std::size_t>(pivot * n + entry)]);
            std::swap(rhs[static_cast<std::size_t>(column)],
                      rhs[static_cast<std::size_t>(pivot)]);
        }
        for (int row = column + 1; row < n; ++row) {
            const Complex factor =
                matrix[static_cast<std::size_t>(row * n + column)] /
                matrix[static_cast<std::size_t>(column * n + column)];
            for (int entry = column + 1; entry < n; ++entry)
                matrix[static_cast<std::size_t>(row * n + entry)] -=
                    factor * matrix[static_cast<std::size_t>(column * n + entry)];
            rhs[static_cast<std::size_t>(row)] -=
                factor * rhs[static_cast<std::size_t>(column)];
        }
    }

    std::vector<Complex> solution(static_cast<std::size_t>(n));
    for (int row = n - 1; row >= 0; --row) {
        Complex value = rhs[static_cast<std::size_t>(row)];
        for (int column = row + 1; column < n; ++column)
            value -= matrix[static_cast<std::size_t>(row * n + column)] *
                     solution[static_cast<std::size_t>(column)];
        solution[static_cast<std::size_t>(row)] =
            value / matrix[static_cast<std::size_t>(row * n + row)];
    }
    return solution;
}

}  // namespace

EddyRuntime::EddyRuntime(EddyConfig config) : config_(std::move(config)) {
    const std::size_t matrix_size =
        static_cast<std::size_t>(config_.n_unknown) * config_.n_unknown;
    const std::size_t temperature_count =
        static_cast<std::size_t>(config_.n_temperature);
    if (config_.n_unknown <= 0 || config_.n_heat <= 0 ||
        config_.n_temperature <= 0 || config_.matrix.size() != matrix_size ||
        config_.rhs_per_amp.size() != static_cast<std::size_t>(config_.n_unknown) ||
        config_.heat_projection.size() !=
            static_cast<std::size_t>(config_.n_heat) * config_.n_unknown ||
        config_.heat_weights.size() != static_cast<std::size_t>(config_.n_heat) ||
        config_.temperature_weights.size() != temperature_count ||
        (!config_.matrix_temperature_slope.empty() &&
         config_.matrix_temperature_slope.size() !=
             temperature_count * matrix_size) ||
        !std::isfinite(config_.reference_temperature_K) ||
        !std::isfinite(config_.angle_origin_rad))
        throw std::invalid_argument("invalid IH Eddy configuration");
    if (!std::all_of(config_.matrix.begin(), config_.matrix.end(), finite) ||
        !std::all_of(config_.rhs_per_amp.begin(), config_.rhs_per_amp.end(), finite) ||
        !std::all_of(config_.matrix_temperature_slope.begin(),
                     config_.matrix_temperature_slope.end(), finite))
        throw std::invalid_argument("IH Eddy operators must contain finite values");
    require_finite(config_.heat_projection, "heat_projection");
    require_positive(config_.heat_weights, "heat_cell_weights");
    require_positive(config_.temperature_weights, "temperature_cell_weights");

    previous_temperature_.assign(temperature_count, 0.0);
    cached_heat_.assign(static_cast<std::size_t>(config_.n_heat), 0.0);
    reference_heat_.assign(static_cast<std::size_t>(config_.n_heat), 0.0);
}

std::vector<double> EddyRuntime::output(
        double current, double angle,
        const std::vector<double>& temperature) {
    if (!std::isfinite(current) || !std::isfinite(angle) ||
        temperature.size() != previous_temperature_.size())
        throw std::invalid_argument("invalid IH Eddy input");
    require_finite(temperature, "temperature_K");

    std::vector<double> local_temperature;
    if (config_.periodic_rotation) {
        transport_periodic(temperature, config_.temperature_weights,
                           angle - config_.angle_origin_rad,
                           local_temperature);
    } else {
        if (!equivalent_angle(angle, config_.angle_origin_rad))
            throw std::invalid_argument(
                "changing angle requires periodic rotation");
        local_temperature = temperature;
    }

    bool temperature_changed = !have_cache_;
    for (std::size_t index = 0; index < local_temperature.size(); ++index)
        temperature_changed = temperature_changed ||
            local_temperature[index] != previous_temperature_[index];

    if (!have_cache_ ||
        (!config_.matrix_temperature_slope.empty() && temperature_changed)) {
        std::vector<Complex> matrix = config_.matrix;
        const std::size_t matrix_size =
            static_cast<std::size_t>(config_.n_unknown) * config_.n_unknown;
        if (!config_.matrix_temperature_slope.empty()) {
            for (std::size_t temperature_index = 0;
                 temperature_index < local_temperature.size();
                 ++temperature_index) {
                const double delta_temperature =
                    local_temperature[temperature_index] -
                    config_.reference_temperature_K;
                for (std::size_t entry = 0; entry < matrix_size; ++entry)
                    matrix[entry] += delta_temperature *
                        config_.matrix_temperature_slope[
                            temperature_index * matrix_size + entry];
            }
        }
        const std::vector<Complex> solution =
            solve(std::move(matrix), config_.rhs_per_amp, config_.n_unknown);
        for (int heat_index = 0; heat_index < config_.n_heat; ++heat_index) {
            double heat = 0.0;
            for (int unknown = 0; unknown < config_.n_unknown; ++unknown)
                heat += config_.heat_projection[static_cast<std::size_t>(
                            heat_index * config_.n_unknown + unknown)] *
                        std::norm(solution[static_cast<std::size_t>(unknown)]);
            reference_heat_[static_cast<std::size_t>(heat_index)] = heat;
        }
        previous_temperature_ = local_temperature;
        have_cache_ = true;
    }

    std::vector<double> heat_at_origin = reference_heat_;
    for (double& heat : heat_at_origin) heat *= current * current;
    if (config_.periodic_rotation) {
        transport_periodic(heat_at_origin, config_.heat_weights,
                           -(angle - config_.angle_origin_rad), cached_heat_);
    } else {
        cached_heat_ = std::move(heat_at_origin);
    }
    return cached_heat_;
}

ThermalRuntime::ThermalRuntime(ThermalConfig config)
    : config_(std::move(config)) {
    const int n = static_cast<int>(config_.initial_temperature_K.size());
    if (config_.n_heat <= 0 || n <= 0 ||
        config_.weights.size() != static_cast<std::size_t>(n) ||
        config_.heat_to_temperature.size() !=
            static_cast<std::size_t>(n) * config_.n_heat ||
        !(config_.options.dt_s > 0.0) ||
        !std::isfinite(config_.options.dt_s) ||
        !(config_.options.tolerance > 0.0) ||
        !std::isfinite(config_.options.tolerance) ||
        config_.options.max_iterations <= 0 ||
        config_.options.convection_W_per_m2K < 0.0 ||
        !std::isfinite(config_.options.convection_W_per_m2K) ||
        !std::isfinite(config_.angle_origin_rad))
        throw std::invalid_argument("invalid IH Thermal configuration");
    require_positive(config_.initial_temperature_K, "initial_temperature_K");
    require_positive(config_.weights, "temperature_cell_weights");
    require_finite(config_.heat_to_temperature,
                   "heat_to_temperature_projection");
    validate_csr(config_.mass, n, "mass");
    validate_csr(config_.stiffness, n, "stiffness");
    if (config_.mass.row_ptr != config_.stiffness.row_ptr ||
        config_.mass.col != config_.stiffness.col)
        throw std::invalid_argument(
            "mass and stiffness CSR sparsity must match");
    if (config_.has_convection) {
        validate_csr(config_.convection, n, "convection");
        if (config_.mass.row_ptr != config_.convection.row_ptr ||
            config_.mass.col != config_.convection.col)
            throw std::invalid_argument(
                "mass and convection CSR sparsity must match");
    }
    reset();
}

void ThermalRuntime::reset() {
    state_.temperature_K = config_.initial_temperature_K;
    state_.time_s = 0.0;
    state_.previous_angle_rad = config_.angle_origin_rad;
    step_count_ = 0;
}

void ThermalRuntime::update(const std::vector<double>& heat, double ambient,
                            double angle) {
    if (heat.size() != static_cast<std::size_t>(config_.n_heat) ||
        !std::isfinite(ambient) || !std::isfinite(angle))
        throw std::invalid_argument("invalid IH Thermal input");
    require_finite(heat, "heat_density_W_per_m3");

    std::vector<double> source(state_.temperature_K.size(), 0.0);
    for (std::size_t temperature_index = 0;
         temperature_index < source.size(); ++temperature_index)
        for (int heat_index = 0; heat_index < config_.n_heat; ++heat_index)
            source[temperature_index] +=
                config_.heat_to_temperature[
                    temperature_index * config_.n_heat + heat_index] *
                heat[static_cast<std::size_t>(heat_index)];

    if (config_.periodic_rotation) {
        std::vector<double> moved;
        transport_periodic(state_.temperature_K, config_.weights,
                           angle - state_.previous_angle_rad, moved);
        state_.temperature_K = std::move(moved);
    } else if (!equivalent_angle(angle, state_.previous_angle_rad)) {
        throw std::invalid_argument(
            "changing angle requires periodic rotation");
    }

    config_.options.ambient_temperature_K = ambient;
    advance_thermal(config_.mass, config_.stiffness,
                    config_.has_convection ? &config_.convection : nullptr,
                    source, config_.weights, angle, config_.options, state_);
    ++step_count_;
}

}}  // namespace radia::ih
